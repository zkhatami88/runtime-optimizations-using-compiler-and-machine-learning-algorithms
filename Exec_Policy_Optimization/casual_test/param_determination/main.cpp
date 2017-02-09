#include <iostream>

#include <hpx/hpx_init.hpp>
//#include <hpx/parallel/algorithms/for_each.hpp>
//#include <hpx/util/high_resolution_timer.hpp>
#include <boost/range.hpp>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <typeinfo>
#include <iterator>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>
#include <hpx/parallel/seq_or_par.hpp>
#include <hpx/parallel/chunk_size_determination.hpp>
//#include "seq_or_par.hpp"
//#include "chunk_size_determination.hpp"

int hpx_main(int argc, char* argv[]) {

	std::vector<std::size_t> f_sp1 = {10, 30, 20, 50, 10, 30};
	std::vector<std::size_t> f_sp2 = {1, 3, 2, 5, 1, 3};

	std::vector<std::size_t> f_pd1 = {1, 3, 2, 5};
	std::vector<std::size_t> f_pd2 = {10, 30, 20, 50};

	bool res1 = hpx::parallel::seq_or_par(f_sp1);
	bool res2 = hpx::parallel::seq_or_par(f_sp2);


	hpx::parallel::dynamic_chunk_size dcs1 = hpx::parallel::chunk_size_determination(f_pd1);
	hpx::parallel::dynamic_chunk_size dcs2 = hpx::parallel::chunk_size_determination(f_pd2);

	std::cout << res1 << ", " << res2 << std::endl;

	return hpx::finalize();
}


int main(int argc, char* argv[]) {
 
    return hpx::init(argc, argv);
}