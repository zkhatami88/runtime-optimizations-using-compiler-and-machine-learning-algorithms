CHECKER_CMD="/home/zahra/clang-llvm/build/bin/loop-convert /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Real_applications/matrix_mult/main.cpp -- -std=c++11 -I/home/zahra/Projects/HPX/repo -I/home/zahra/Projects/HPX/build -I/home/zahra/Projects/HPX/repo/tests -I/home/zahra/Projects/boost_1_63_0 "

cd build

$CHECKER_CMD
make -j12
./main -t1
