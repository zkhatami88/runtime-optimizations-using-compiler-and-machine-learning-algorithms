set(LLVM_LINK_COMPONENTS support)
add_compile_options(-std=c++14 -Wall -Wextra -Wno-unused-parameter)

add_clang_executable(for-each-checker
  ForEachChecker.cpp
  )
target_link_libraries(for-each-checker
  clangTooling
  clangBasic
  clangASTMatchers
  )
