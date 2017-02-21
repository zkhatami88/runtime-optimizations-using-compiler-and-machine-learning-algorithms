#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendAction.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "llvm/Support/CommandLine.h"

using namespace llvm;
using namespace clang;
using namespace clang::tooling;

static QualType CleanUpType(QualType type) {
    type = type.getNonReferenceType();
    type.removeLocalConst();
    type.removeLocalRestrict();
    type.removeLocalVolatile();
    return type;
}

static void PrintType(const std::string& print_name, const std::string& name, QualType type, int depth = 0) {
    std::string prefix;
    prefix.append(depth * 4, ' ');
    if(depth == 0) {
        prefix += print_name + ":" + name;
    }
    if(isa<const TemplateSpecializationType>(type.getTypePtr())) {
        const auto *template_specialization_type = dyn_cast<const TemplateSpecializationType>(type.getTypePtr());
        for (TemplateSpecializationType::iterator i = template_specialization_type->begin(); i != template_specialization_type->end(); ++i) {
            QualType subType = CleanUpType(i->getAsType());
            printf("%s - Template %s\n", prefix.c_str(), template_specialization_type->getTemplateName().getAsTemplateDecl()->getNameAsString().c_str());
            PrintType(print_name, name, subType, depth + 1);
        }
    } else {
        printf("%s - %s\n", prefix.c_str(), CleanUpType(type).getAsString().c_str());
    }
}

class FindNamedClassVisitor : public RecursiveASTVisitor<FindNamedClassVisitor> {
public:
  explicit FindNamedClassVisitor(ASTContext *Context) : Context(Context) {}

  bool VisitCXXRecordDecl(CXXRecordDecl *Declaration) {
    std::string s = Declaration->getQualifiedNameAsString();
    if(s.find("TestClass") == -1) {
        return true;
    }
    printf("\n\nCLASS: %s\n", Declaration->getNameAsString().c_str());
    
    auto ns = Declaration->getEnclosingNamespaceContext();
    while(ns->isNamespace()) {
        NamespaceDecl* nsd = dyn_cast_or_null<NamespaceDecl>(ns);
        printf("Namespace: %s\n", nsd->getNameAsString().c_str());
        ns = ns->getParent();
    }
    for(auto field = Declaration->field_begin();field!=Declaration->field_end();++field) {
        std::string param_name = (*field)->getName();
        PrintType("Field", param_name, CleanUpType((*field)->getTypeSourceInfo()->getType()));
    }
    printf("\n");
    return true;
  }
  
  bool VisitCXXMethodDecl(CXXMethodDecl* Declaration) {
    std::string s = Declaration->getQualifiedNameAsString();
    if(s.find("TestClass") == -1) {
        return true;
    }
    bool isConstructorOrDestructor = isa<CXXConstructorDecl>(Declaration) || isa<CXXDestructorDecl>(Declaration);
    if(isConstructorOrDestructor || Declaration->isMoveAssignmentOperator() || Declaration->isCopyAssignmentOperator()) {
        return true;
    }
    
    std::string name = Declaration->getAsFunction()->getNameAsString();
    if(Declaration->isVirtual()) {
        printf("Virtual function: %s\n", name.c_str());
    } else if(Declaration->isStatic()) {
        printf("Static function: %s\n", name.c_str());
    } else if(Declaration->isCXXClassMember()) {
        printf("Member function: %s\n", name.c_str());
    } else {
        printf("Other function: %s\n", name.c_str());
    }
    
    
    for(auto param = Declaration->param_begin();param!=Declaration->param_end();++param) {
        std::string param_name = (*param)->getName();
        PrintType("Parameter", param_name, CleanUpType((*param)->getTypeSourceInfo()->getType()));
    }
    printf("\n");
    return true;
  }

private:
  ASTContext *Context;
};

class FindNamedClassConsumer : public clang::ASTConsumer {
public:
  explicit FindNamedClassConsumer(ASTContext *Context) : Visitor(Context) {}

  virtual void HandleTranslationUnit(clang::ASTContext &Context) {
    Visitor.TraverseDecl(Context.getTranslationUnitDecl());
  }
private:
  FindNamedClassVisitor Visitor;
};

class FindNamedClassAction : public clang::ASTFrontendAction {
public:
  virtual std::unique_ptr<clang::ASTConsumer> CreateASTConsumer(clang::CompilerInstance &Compiler, llvm::StringRef InFile) {
    return std::unique_ptr<clang::ASTConsumer>(new FindNamedClassConsumer(&Compiler.getASTContext()));
  }
};

// Apply a custom category to all command-line options so that they are the
// only ones displayed.
static llvm::cl::OptionCategory MyToolCategory("my-tool options");

// CommonOptionsParser declares HelpMessage with a description of the common
// command-line options related to the compilation database and input files.
// It's nice to have this help message in all tools.
static cl::extrahelp CommonHelp(CommonOptionsParser::HelpMessage);

// A help message for this specific tool can be added afterwards.
static cl::extrahelp MoreHelp("\nMore help text...");

int main(int argc, const char **argv) {
  CommonOptionsParser OptionsParser(argc, argv, MyToolCategory);
  ClangTool Tool(OptionsParser.getCompilations(),
                 OptionsParser.getSourcePathList());
  return Tool.run(newFrontendActionFactory<FindNamedClassAction>().get());
}