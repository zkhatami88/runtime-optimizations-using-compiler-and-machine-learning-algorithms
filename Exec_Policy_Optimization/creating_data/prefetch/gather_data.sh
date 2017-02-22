CHECKER_CMD="/home/zahra/clang-llvm/build/bin/loop-convert /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/prefetch/main.cpp -- -std=c++11 -I/home/zahra/Projects/HPX/repo -I/home/zahra/Projects/HPX/build -I/home/zahra/Projects/HPX/repo/tests -I/home/zahra/Projects/boost_1_63_0 "
FILE="../input/data_new.dat"
MAX_NUM_THREADS=32

cd build

#rm $FILE
echo "num_ops num_float_ops num_comp_ops deepest_loop_level num_threads num_lambda_iter t(1) t(5) t(10) t(500) t(100) t(500)" >> $FILE

for file in ../functions/*
do
	if [ "$file" != "../functions/func.hpp" ]
	then
		echo $file
		cp $file "../functions/func.hpp"
		echo -n "$file" >> $FILE
		$CHECKER_CMD  >> $FILE
		echo -n " " >> $FILE
		make -j12
		./main -t1 >> $FILE

		for (( COUNTER=2; COUNTER<=MAX_NUM_THREADS; COUNTER+=2 )); do
			echo -n "$file" >> $FILE
			$CHECKER_CMD  >> $FILE
			echo -n " " >> $FILE
			./main -t$COUNTER >> $FILE
		done

	fi
done