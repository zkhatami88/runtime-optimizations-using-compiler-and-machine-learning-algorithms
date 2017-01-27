#include <string>

#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
#include "llvm/Support/raw_ostream.h"

using namespace clang;
using namespace clang::ast_matchers;
using namespace clang::driver;
using namespace clang::tooling;

struct statistics
{
    // Below static information related to the values of features are going to passed
    // through initilizer_list as the first argument in hpx::parallel::for_each:
    //
    // hpx::parallel::for_each(hpx::parallel::feature_container<double>({}), //:features
    //                            hpx::parallel::weights_container<double>({w0, w1, w2, w3, w4, w5, w6}), //:weights
    //                            policy, range.begin(), range.end(), lambda_fnc);
    // 
    // , which are going to be determined and assigned with this class

    unsigned num_threads;               //f0 = number of threads: will be assigned with runtime
    unsigned num_ops;                   //f1
    unsigned num_float_ops;             //f2
    unsigned num_comparison_ops;        //f3
    unsigned num_lambda_iterations;     //f4 : will be re-assigned with runtime
    unsigned deepest_loop_level;        //f5
    unsigned num_int_variables;
    unsigned num_float_variables;
    unsigned num_if_stmts;
    unsigned num_if_stmts_in_loop;
    unsigned num_func_calls;
    unsigned num_func_calls_in_loop;
    
    
    statistics() : num_threads(0), num_ops(0), num_float_ops(0), num_comparison_ops(0),
                    num_lambda_iterations(0), deepest_loop_level(0),
                    num_int_variables(0), num_float_variables(0),
                    num_if_stmts(0), num_if_stmts_in_loop(0),
                    num_func_calls(0), num_func_calls_in_loop(0) {}
    
    friend raw_ostream& operator<< (raw_ostream &out, statistics const& s)
    {
        return out  << s.num_threads << " " << s.num_ops << " " << s.num_float_ops
                    << " " << s.num_comparison_ops << " " << s.num_lambda_iterations
                    << " " << s.deepest_loop_level; /* << " " << s.num_int_variables
                    << " " << s.num_float_variables << " " << s.num_if_stmts
                    << " " << s.num_if_stmts_in_loop << " " << s.num_func_calls
                    << " " << s.num_func_calls_in_loop; */                                     
    }
};

namespace detail
{
void analyze_binary_operator(BinaryOperator* b, unsigned nested_for_loop_factor,
                                statistics& stats)
{
    stats.num_ops += nested_for_loop_factor;
                                
    Expr* lhs = b->getLHS();
    Expr* rhs = b->getRHS();
    
    if (lhs->getType()->isFloatingType()
            || rhs->getType()->isFloatingType())
        stats.num_float_ops += nested_for_loop_factor;
    
    if (b->isComparisonOp())
        stats.num_comparison_ops += nested_for_loop_factor;
}

void analyze_unary_operator(UnaryOperator* u, unsigned nested_for_loop_factor,
                                statistics& stats)
{
    stats.num_ops += nested_for_loop_factor;
}

void analyze_decl_statement(DeclStmt* s, unsigned nested_for_loop_factor,
                                statistics& stats)
{
    if (s->isSingleDecl())
    {
        // maybe multiply by nested_for_loop_factor at some point
        if (isa<IntegerLiteral>(*s->child_begin()))
            ++stats.num_int_variables;
        else if (isa<FloatingLiteral>(*s->child_begin()))
            ++stats.num_float_variables;
    }
}

unsigned analyze_for_loop(ForStmt* for_stmt, statistics& stats)
{    
    ++stats.num_comparison_ops;
    stats.num_ops += 2;
    
    // get the initial value
    const Stmt* init_stmt = for_stmt->getInit(); 
    const DeclStmt* init_decl = cast<DeclStmt>(init_stmt);
    const IntegerLiteral* init_literal =
        cast<IntegerLiteral>(*(init_decl->child_begin()));
    llvm::APInt initial_value = init_literal->getValue();

    // get the condition value
    const Expr* cond_expr = for_stmt->getCond();
    const BinaryOperator* cond_op = cast<BinaryOperator>(cond_expr);

    const IntegerLiteral* cond_literal =
        cast<IntegerLiteral>(cond_op->getRHS());
    const llvm::APInt cond_value = cond_literal->getValue();
    
    // condition code
    const StringRef cond_opcode = cond_op->getOpcodeStr();
    StringRef inc_opcode;
    llvm::APInt inc_value;

    const Expr* inc_expr = for_stmt->getInc();
    // handle cases with non unit stride i += n
    if (isa<BinaryOperator>(inc_expr))
    {
        const BinaryOperator* b = cast<BinaryOperator>(inc_expr);

        inc_opcode = b->getOpcodeStr();
        const IntegerLiteral* inc_literal = cast<IntegerLiteral>(b->getRHS());
        inc_value = inc_literal->getValue();
    }
    // handle unit stride i++
    else
    {
        const UnaryOperator* u = cast<UnaryOperator>(inc_expr);
        if (u->isIncrementOp())
            inc_opcode = "++";
        else
            inc_opcode = "--";
    }

    unsigned num_iter = 0;

    while(true)
    {
        if (cond_opcode == "<" && initial_value.uge(cond_value))
            break;
        else if (cond_opcode == "<=" && initial_value.ugt(cond_value))
            break;                
        else if (cond_opcode == ">" && !initial_value.ule(cond_value))
            break;     
        else if (cond_opcode == ">=" && !initial_value.ugt(cond_value))
            break;                
        else if (cond_opcode == "!=" && initial_value == cond_value)
            break;
                
        if (inc_opcode == "++")
            initial_value++;
        else if (inc_opcode == "--")
            initial_value--;
        else if (inc_opcode == "+=")
            initial_value += inc_value;
        else if (inc_opcode == "-=")
            initial_value -= inc_value;
        else if (inc_opcode == "*=")
            initial_value *= inc_value;

        num_iter++;
    }

    return num_iter;
}
}

static llvm::cl::OptionCategory MatcherSampleCategory("Matcher Sample");

class ForEachCallHandler : public MatchFinder::MatchCallback {
public:
  ForEachCallHandler(Rewriter &Rewrite) : Rewrite(Rewrite) {}

  virtual void run(const MatchFinder::MatchResult &Result) {

    if (const CallExpr *call =
            Result.Nodes.getNodeAs<clang::CallExpr>("functionCall"))
    {           
        statistics stats;
        const SourceManager &Sources = *Result.SourceManager;
      
        const Expr* num_iters_expr = call->getArg(4);
        
        if (isa<IntegerLiteral>(num_iters_expr))
            stats.num_lambda_iterations =
                cast<IntegerLiteral>(num_iters_expr)->getValue().getZExtValue();
      
        const Expr* lambda_expr = call->getArg(5);
     
        const CXXRecordDecl* lambda_record =
            lambda_expr->getBestDynamicClassType();
        
        if (!lambda_record->isLambda())
            return;
            
        const CXXMethodDecl* lambda_callop =
            lambda_record->getLambdaCallOperator();
        
        Stmt* lambda_body = lambda_callop->getBody();
        
        analyze_statement(lambda_body, stats);

        //printing out the extracted data
        llvm::outs() << stats;
       
        // Replacing first argument (arg(0)) to include these informations:
        // num_threads, num_ops, num_float_ops, num_comparison_ops, num_lambda_iterations and deepest_loop_level:
        const Expr* feature_container = call->getArg(0);
        const Expr* next_arg = call->getArg(1);
        including_stats_info_in_arg0(Sources, feature_container, next_arg, stats);        

        Rewrite.overwriteChangedFiles();   
    }
  }

private:
  Rewriter &Rewrite;
    
  void analyze_statement(Stmt* s, statistics& stats,
                            unsigned nested_for_loop_factor = 1,
                            unsigned loop_level = 0)
  {
    for (Stmt* child_stmt : s->children())
    {
        unsigned tmp_loop_factor = nested_for_loop_factor;
        
        if (!child_stmt)
            continue;
                    
        if (isa<BinaryOperator>(child_stmt))
            detail::analyze_binary_operator(cast<BinaryOperator>(child_stmt),
                                        tmp_loop_factor, stats); 
        else if (isa<UnaryOperator>(child_stmt))
            detail::analyze_unary_operator(cast<UnaryOperator>(child_stmt),
                                        tmp_loop_factor, stats);
        else if (isa<ForStmt>(child_stmt))
        {
            tmp_loop_factor *=
                detail::analyze_for_loop(cast<ForStmt>(child_stmt), stats);
            ++loop_level;
            
            stats.deepest_loop_level =
                std::max(stats.deepest_loop_level, loop_level);            
        }
        else if (isa<DeclStmt>(child_stmt))
        {
            detail::analyze_decl_statement(cast<DeclStmt>(child_stmt),
                                tmp_loop_factor, stats);
        }        
        else if (isa<IfStmt>(child_stmt))
        {
            stats.num_if_stmts += tmp_loop_factor;    
            if (tmp_loop_factor > 1)
                stats.num_if_stmts_in_loop += tmp_loop_factor;
        }
        else if(isa<CallExpr>(child_stmt))
        {
            stats.num_func_calls += tmp_loop_factor;
            if (tmp_loop_factor > 1)
                stats.num_func_calls_in_loop += tmp_loop_factor;
        }
        
        
        analyze_statement(child_stmt, stats, tmp_loop_factor,
                            loop_level);
    } 
  }

    void including_stats_info_in_arg0(const SourceManager &Sources, const Expr* feature_container, 
                                        const Expr* next_arg, statistics& stats) 
    {
        //setting arg0 to include static information:
        std::string data = "hpx::parallel::features_container<double>({" + std::to_string(stats.num_threads) + 
                            ", " + std::to_string(stats.num_ops) + 
                            ", " + std::to_string(stats.num_float_ops) + 
                            ", " + std::to_string(stats.num_comparison_ops) + 
                            ", " + std::to_string(stats.num_lambda_iterations) + 
                            ", " + std::to_string(stats.deepest_loop_level) + "}), ";

        StringRef features = StringRef(data);
        SourceLocation arg0_location = feature_container->getLocStart();
        SourceLocation arg1_location = next_arg->getLocStart();

        unsigned length = Sources.getSpellingColumnNumber(arg1_location) - 
            Sources.getSpellingColumnNumber(arg0_location);

        Rewrite.ReplaceText(arg0_location, length, features);
    }
};

// Implementation of the ASTConsumer interface for reading an AST produced
// by the Clang parser. It registers a couple of matchers and runs them on
// the AST.
class MyASTConsumer : public ASTConsumer {
public:
  MyASTConsumer(Rewriter &R) : HandlerForEach(R) {
    // Add a simple matcher for finding adaptive_for_each //'for_each' and 'for_each_n' statements.
    Matcher.addMatcher(callExpr(
          callee(functionDecl(hasName("hpx::parallel::adaptive_for_each")))
          
        ).bind("functionCall"), &HandlerForEach);
    
    /*
    Matcher.addMatcher(callExpr(
          callee(functionDecl(hasName("hpx::parallel::for_each_n")))
          
        ).bind("functionCall"), &HandlerForEach);*/
  }

  void HandleTranslationUnit(ASTContext &Context) override {
    // Run the matchers when we have the whole TU parsed.
    Matcher.matchAST(Context);
  }

private:
  ForEachCallHandler HandlerForEach;
  MatchFinder Matcher;
};

// For each source file provided to the tool, a new FrontendAction is created.
class MyFrontendAction : public ASTFrontendAction {
public:
  MyFrontendAction() {}
  void EndSourceFileAction() override {
   /* TheRewriter.getEditBuffer(TheRewriter.getSourceMgr().getMainFileID())
        .write(llvm::outs());*/
  }

  std::unique_ptr<ASTConsumer> CreateASTConsumer(CompilerInstance &CI,
                                                 StringRef file) override {
    TheRewriter.setSourceMgr(CI.getSourceManager(), CI.getLangOpts());
    return llvm::make_unique<MyASTConsumer>(TheRewriter);
  }

private:
  Rewriter TheRewriter;
};

int main(int argc, const char **argv) {
  CommonOptionsParser op(argc, argv, MatcherSampleCategory);
  ClangTool Tool(op.getCompilations(), op.getSourcePathList());

  return Tool.run(newFrontendActionFactory<MyFrontendAction>().get());
}
