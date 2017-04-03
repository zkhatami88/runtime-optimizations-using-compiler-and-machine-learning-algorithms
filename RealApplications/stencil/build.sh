CHECKER_CMD="/home/zahra/clang-llvm/build/bin/loop-convert ../main.cpp -- -std=c++11 -I/home/zahra/Projects/HPX/repo -I/home/zahra/Projects/HPX/build -I/home/zahra/Projects/HPX/repo/tests -I/home/zahra/Projects/boost_1_63_0"


	
cd /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/RealApplications/stencil
mkdir build && cd build
cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ -DHPX_DIR=/home/zahra/Projects/HPX/build/lib/cmake/HPX -DCMAKE_BUILD_TYPE=Release -std=c++11 ..
	
$CHECKER_CMD

make -j
	
for (( COUNTER=1; COUNTER<=16; COUNTER+=2 )); do
	./main -t$COUNTER
done
