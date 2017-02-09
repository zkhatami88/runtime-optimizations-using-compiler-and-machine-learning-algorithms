CHECKER_CMD="/home/zahra/clang-llvm/build/bin/loop-convert /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/other_tests/seq_par/test2/main.cpp -- -std=c++11 -I/home/zahra/Projects/HPX/repo -I/home/zahra/Projects/HPX/build -I/home/zahra/Projects/HPX/repo/tests -I/home/zahra/Projects/boost_1_63_0 "
FILE="../computed_features/features.txt"

cd build

#rm $FILE
echo "num_ops num_float_ops num_comp_ops deepest_loop_level" >> $FILE

file="my_test_cases/test1.hpp" 
echo $file
$CHECKER_CMD  >> $FILE
echo -n " " >> $FILE
make -j12
./main -t1 >> $FILE
