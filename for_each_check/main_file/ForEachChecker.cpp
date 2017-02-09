//  Copyright (c) 2017 Lukas Troska and Zahra Khatami

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
#include "clang/Lex/Lexer.h"
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

    unsigned num_threads;               //f0 = will be assigned with runtime
    unsigned num_lambda_iterations;     //f1 = will be re-assigned with runtime
    unsigned num_ops;                   //f2
    unsigned num_float_ops;             //f3
    unsigned num_comparison_ops;        //f4    
    unsigned deepest_loop_level;        //f5
    unsigned num_int_variables;
    unsigned num_float_variables;
    unsigned num_if_stmts;
    unsigned num_if_stmts_in_loop;
    unsigned num_func_calls;
    unsigned num_func_calls_in_loop;
    
    
    statistics() : num_threads(0), num_lambda_iterations(0), 
                    num_ops(0), num_float_ops(0), 
                    num_comparison_ops(0), deepest_loop_level(0),
                    num_int_variables(0), num_float_variables(0),
                    num_if_stmts(0), num_if_stmts_in_loop(0),
                    num_func_calls(0), num_func_calls_in_loop(0) {}
    
    friend raw_ostream& operator<< (raw_ostream &out, statistics const& s)
    {
        return out  << /* "\n" << s.num_threads
                    << " " << s.num_lambda_iterations 
                    <<*/ " " << s.num_ops 
                    << " " << s.num_float_ops
                    << " " << s.num_comparison_ops 
                    << " " << s.deepest_loop_level;/* << "\n"; /* << " " << s.num_int_variables
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
        if (call->getNumArgs() == 0)
            return;
        const Expr* exec_policy_expr = call->getArg(0);
                
        if (!isa<DeclRefExpr>(exec_policy_expr))
            return;

        statistics stats;
        const SourceManager *SM = Result.SourceManager;

        //llvm::outs() << call->getLocStart().printToString(Result.Context->getSourceManager()) << "\n num_arg: " << call->getNumArgs() <<"\n";

        if (call->getNumArgs() >= 4 ) {
            /*
            //for tests
            for (auto it = call->arg_begin(); it != call->arg_end(); ++it)
            {   
                if (isa<DeclRefExpr>(*it))
                {
                    const Expr* lambda_expr = *it;
                                        
                    const CXXRecordDecl* lambda_record =
                        lambda_expr->getBestDynamicClassType();
                
                    if (!lambda_record->isLambda())
                        continue;            
                                                            
                    const CXXMethodDecl* lambda_callop =
                        lambda_record->getLambdaCallOperator();
                    
                    Stmt* lambda_body = lambda_callop->getBody();                                        
                    analyze_statement(lambda_body, stats);

                    // Printing out the extracted data
                    llvm::outs() << stats;
                }
            }
            */
            ///////////////////////////////////////////////////////////////////
            
            // Determining if policy is par_if
            SourceRange par_if_policy(call->getArg(0)->getExprLoc(), call->getArg(1)->getExprLoc().getLocWithOffset(-2));
                                      
            std::string par_if_string =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(par_if_policy), *SM,
                    LangOptions()
                ).str(); 

            //llvm::outs() << "\n first arg is : " << par_if_string << "\n"; //with param or exec : policy

            if (par_if_string.find("par_if") != std::string::npos)
            {
      
                llvm::outs() << "\nFound par_if exec policy in call at line "
                << call->getLocStart().printToString(
                                            Result.Context->getSourceManager())
                << "\n-------------------------\n";
                policy_determination(call, SM, stats);  
            }
            ///////////////////////////////////////////////////////////////////

            // Determining if policy is which_chunk
            SourceRange chunk_policy(call->getArg(0)->getExprLoc(), call->getArg(1)->getExprLoc().getLocWithOffset(-2));
                                      
            std::string chunk_policy_string =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(chunk_policy), *SM,
                    LangOptions()
                ).str(); 

            if(chunk_policy_string.find("which_chunk") != std::string::npos)
            {      
                llvm::outs() << "\nFound which_chunk exec policy in call at line "
                << call->getLocStart().printToString(
                                            Result.Context->getSourceManager())
                << "\n-------------------------\n";
                chunk_size_determination(call, SM, stats);  
            }
            ///////////////////////////////////////////////////////////////////

            // Determining if policy is which_prefetcher 
        }
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

  void chunk_size_determination(const CallExpr *call, const SourceManager *SM, statistics& stats) {
    for (auto it = call->arg_begin(); it != call->arg_end(); ++it)
    {                
        if (isa<DeclRefExpr>(*it))
        {
            const Expr* lambda_expr = *it;
                                        
            const CXXRecordDecl* lambda_record =
                lambda_expr->getBestDynamicClassType();
                
            if (!lambda_record->isLambda())
                continue;            
                                                            
            const CXXMethodDecl* lambda_callop =
                lambda_record->getLambdaCallOperator();
                    
            Stmt* lambda_body = lambda_callop->getBody();                                        
            analyze_statement(lambda_body, stats);

            // Printing out the extracted data
            //llvm::outs() << stats;     

            // Get the source range and manager.
            SourceRange range1 = call->getSourceRange();
            range1.setEnd(call->getArg(0)->getExprLoc());
            SourceRange range2 = call->getSourceRange();
            range2.setBegin(call->getArg(1)->getExprLoc().getLocWithOffset(-1));
            range2.setEnd(range2.getEnd().getLocWithOffset(2));
                    
            SourceRange iter_first_range(call->getArg(1)->getExprLoc(), call->getArg(2)->getExprLoc().getLocWithOffset(-2));
            SourceRange iter_last_range(call->getArg(2)->getExprLoc(), call->getArg(3)->getExprLoc().getLocWithOffset(-2));
                                      
            std::string first_iter =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(iter_first_range), *SM,
                    LangOptions()
                ).str();                
                        
            std::string last_iter =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(iter_last_range), *SM,
                    LangOptions()
                ).str();                    

            // Use LLVM's lexer to get source text.
            llvm::StringRef ref1 =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(range1), *SM,
                    LangOptions()
                );
            llvm::StringRef ref2 =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(range2), *SM,
                    LangOptions()
                );

            // Passing static info extracted at compile time into runtime
            std::string data = "{hpx::get_os_thread_count(), " + 
                                std::to_string(stats.num_ops) + 
                                ", " + std::to_string(stats.num_float_ops) + 
                                ", " + std::to_string(stats.num_comparison_ops) + 
                                ", std::size_t(std::distance(" + first_iter + ", " + last_iter + "))" +
                                ", " + std::to_string(stats.deepest_loop_level) + "}";

            std::string new_call;

            // Examining code if current policy has any executor to be reattached
            SourceRange policy_range(call->getArg(0)->getExprLoc(), call->getArg(1)->getExprLoc().getLocWithOffset(-2));
            std::string policy_range_string =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(policy_range), *SM,
                    LangOptions()
                ).str();

            std::size_t pos = policy_range_string.find(".on");
            std::string param_exec;
            if (pos != std::string::npos) {
                std::size_t next_pos = policy_range_string.find(")", (pos + 1));
                param_exec = policy_range_string.substr(pos, (next_pos - pos + 1));
            }
            else {
                param_exec = "";
            }

            if (param_exec == "") {                
            
                new_call = "//DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:"
                                    "\n\t" + ref1.str() +
                                    "hpx::parallel::par.with(hpx::parallel::chunk_size_determination(" + 
                                    data + ")), " + 
                                    ref2.str();
            }
            else {

                new_call = "//DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:"
                                    "\n\t" + ref1.str() +
                                    "hpx::parallel::par.with(hpx::parallel::chunk_size_determination(" + 
                                    data + "))" + param_exec + 
                                    ", " + ref2.str();
            }

            Rewrite.ReplaceText(SourceRange(range1.getBegin(), range2.getEnd()), new_call);
            Rewrite.overwriteChangedFiles();
        }            
    }
  }

  void policy_determination(const CallExpr *call, const SourceManager *SM, statistics& stats) {
    for (auto it = call->arg_begin(); it != call->arg_end(); ++it)
    {                
        if (isa<DeclRefExpr>(*it))
        {
            const Expr* lambda_expr = *it;
                                        
            const CXXRecordDecl* lambda_record =
                lambda_expr->getBestDynamicClassType();
                
            if (!lambda_record->isLambda())
                continue;            
                                                            
            const CXXMethodDecl* lambda_callop =
                lambda_record->getLambdaCallOperator();
                    
            Stmt* lambda_body = lambda_callop->getBody();                                        
            analyze_statement(lambda_body, stats);

            //printing out the extracted data
            //llvm::outs() << stats;     

            //Get the source range and manager.
            SourceRange range1 = call->getSourceRange();
            range1.setEnd(call->getArg(0)->getExprLoc());
            SourceRange range2 = call->getSourceRange();
            range2.setBegin(call->getArg(1)->getExprLoc().getLocWithOffset(-1));
            range2.setEnd(range2.getEnd().getLocWithOffset(2));
                    
            SourceRange iter_first_range(call->getArg(1)->getExprLoc(), call->getArg(2)->getExprLoc().getLocWithOffset(-2));
            SourceRange iter_last_range(call->getArg(2)->getExprLoc(), call->getArg(3)->getExprLoc().getLocWithOffset(-2));
                                      
            std::string first_iter =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(iter_first_range), *SM,
                    LangOptions()
                ).str();                
                        
            std::string last_iter =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(iter_last_range), *SM,
                    LangOptions()
                ).str();                    

            //Use LLVM's lexer to get source text.
            llvm::StringRef ref1 =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(range1), *SM,
                    LangOptions()
                );
            llvm::StringRef ref2 =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(range2), *SM,
                    LangOptions()
                );

            // Passing static info extracted at compile time into runtime
            std::string data = "{hpx::get_os_thread_count(), " + 
                            std::to_string(stats.num_ops) + 
                            ", " + std::to_string(stats.num_float_ops) + 
                            ", " + std::to_string(stats.num_comparison_ops) + 
                            ", std::size_t(std::distance(" + first_iter + ", " + last_iter + "))" +
                            ", " + std::to_string(stats.deepest_loop_level) + "}";

            std::string new_call;

            // Examining code if current policy has any executor or parameters to be reattached
            SourceRange policy_range(call->getArg(0)->getExprLoc(), call->getArg(1)->getExprLoc().getLocWithOffset(-2));
            std::string policy_range_string =
                Lexer::getSourceText(
                    CharSourceRange::getCharRange(policy_range), *SM,
                    LangOptions()
                ).str();

            std::size_t pos = policy_range_string.find(".");
            std::string param_exec;
            if (pos != std::string::npos) {
                param_exec = policy_range_string.substr(pos);
            }
            else {
                param_exec = "";
            }            

            if (param_exec == "") {
                    
                new_call = //adding exec and param to be attached
                    "//DETERMING CHUNK SIZE BASED ON STATIC AND DYNAMIC FEATURES:\n \tif (hpx::parallel::seq_or_par(" + 
                    data + ")) \n \t \t" + ref1.str() +
                    "hpx::parallel::seq," + ref2.str() +
                    "\n \telse \n \t \t"  + ref1.str() +
                    "hpx::parallel::par," + ref2.str();
            }
            else {

                new_call = //adding exec and param to be attached
                    "//DETERMING CHUNK SIZE BASED ON STATIC AND DYNAMIC FEATURES:\n \tif (hpx::parallel::seq_or_par(" + 
                    data + ")) \n \t \t" + ref1.str() +
                    "hpx::parallel::seq." + param_exec + "," + ref2.str() +
                    "\n \telse \n \t \t"  + ref1.str() +
                    "hpx::parallel::par." + param_exec + "," + ref2.str();
            }
                    
            Rewrite.ReplaceText(SourceRange(range1.getBegin(), range2.getEnd()), new_call);
            Rewrite.overwriteChangedFiles();
        }            
    }
  }

};

// Implementation of the ASTConsumer interface for reading an AST produced
// by the Clang parser. It registers a couple of matchers and runs them on
// the AST.
class MyASTConsumer : public ASTConsumer {
public:
  MyASTConsumer(Rewriter &R) : HandlerForEach(R) {
    // Add a simple matcher for finding hpx::parallel::* function calls
    Matcher.addMatcher(callExpr(
          callee(functionDecl(matchesName("hpx::parallel::[^:]*")))
          
        ).bind("functionCall"), &HandlerForEach);   
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
