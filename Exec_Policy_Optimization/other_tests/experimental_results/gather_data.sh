CHECKER_CMD="/home/zahra/clang-llvm/build/bin/loop-convert /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/other_tests/experimental_results/main.cpp -- -std=c++11 -I/home/zahra/Projects/HPX/repo -I/home/zahra/Projects/HPX/build -I/home/zahra/Projects/HPX/repo/tests -I/home/zahra/Projects/boost_1_63_0"
FILE="../execution_times/par_seq.dat"

cd build

#rm $FILE
echo "filename num_threads num_lambda_iter num_ops num_float_ops num_comp_ops deepest_loop_level time_par(s)" >> $FILE

file="functions/mat_mat.hpp" 
echo $file
$CHECKER_CMD  >> $FILE
echo -n " " >> $FILE
make -j12
./main -t8 >> $FILE

