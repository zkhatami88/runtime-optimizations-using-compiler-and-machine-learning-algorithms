CHECKER_CMD="/home/lukas/compiler/build/bin/for-each-checker /home/lukas/code/hpx/analyze.cpp -- -std=c++14"
FILE="../data.dat"
MAX_NUM_THREADS=32

cd build

#rm $FILE
echo "filename num_threads num_ops num_float_ops num_comp_ops num_lambda_iter deepest_loop_level num_int_var num_float_var num_if_stmts num_if_stmts_in_loop time_seq(s) time_par(s)" >> $FILE

for file in ../functions/*
do
	if [ "$file" != "../functions/func.hpp" ]
	then
		echo $file
		cp $file "../functions/func.hpp"
		echo -n "$file 1 " >> $FILE
		$CHECKER_CMD  >> $FILE
		echo -n " " >> $FILE
		make -j12
		./main -t1 >> $FILE

		for (( COUNTER=2; COUNTER<=MAX_NUM_THREADS; COUNTER+=2 )); do
			echo -n "$file $COUNTER " >> $FILE
			$CHECKER_CMD  >> $FILE
			echo -n " " >> $FILE
			./main -t$COUNTER >> $FILE
		done

	fi
done
