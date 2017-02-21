//  Copyright (c) 2017 Zahra Khatami 
//  Copyright (c) 2016 David Pfander
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <stdlib.h> 
#include <vector>

#include <hpx/hpx_init.hpp>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/util/high_resolution_timer.hpp>
#include <boost/range.hpp>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <typeinfo>
#include <iterator>
#include <hpx/parallel/executors/sequential_execution_wrapper.hpp>
#include <hpx/parallel/seq_or_par.hpp>
#include <hpx/parallel/chunk_size_determination.hpp>
#include <hpx/parallel/prefetching_distance_determination.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>
#include <hpx/parallel/executors/adaptive_chunk_size.hpp>


#define vector_size 5

template<typename T>
void comparing_perfromances(std::vector<T>& A, std::vector<T>& B, std::vector<T>& C) {

    auto time_range = boost::irange(0, 5);
    
    auto f = [&](int i) {
        // i += 4
        if(i % 4 == 0) {
            for (int j = 0; j < vector_size; j += 4) {

                T result1 = 0.0;
                T result2 = 0.0;
                T result3 = 0.0;
                T result4 = 0.0;

                T result5 = 0.0;
                T result6 = 0.0;
                T result7 = 0.0;
                T result8 = 0.0;

                T result9 = 0.0;
                T result10 = 0.0;
                T result11 = 0.0;
                T result12 = 0.0;

                T result13 = 0.0;
                T result14 = 0.0;
                T result15 = 0.0;
                T result16 = 0.0;

                for (int k = 0; k < vector_size; k++) {                  

                    result1 += A[i * vector_size + k] * B[j * vector_size + k];
                    result2 += A[i * vector_size + k] * B[(j + 1) * vector_size + k];
                    result3 += A[i * vector_size + k] * B[(j + 2) * vector_size + k];
                    result4 += A[i * vector_size + k] * B[(j + 3) * vector_size + k];

                    result5 += A[(i + 1) * vector_size + k] * B[j * vector_size + k];
                    result6 += A[(i + 1) * vector_size + k] * B[(j + 1) * vector_size + k];
                    result7 += A[(i + 1) * vector_size + k] * B[(j + 2) * vector_size + k];
                    result8 += A[(i + 1) * vector_size + k] * B[(j + 3) * vector_size + k];

                    result9 += A[(i + 2) * vector_size + k] * B[j * vector_size + k];
                    result10 += A[(i + 2) * vector_size + k] * B[(j + 1) * vector_size + k];
                    result11 += A[(i + 2) * vector_size + k] * B[(j + 2) * vector_size + k];
                    result12 += A[(i + 2) * vector_size + k] * B[(j + 3) * vector_size + k];

                    result13 += A[(i + 3) * vector_size + k] * B[j * vector_size + k];
                    result14 += A[(i + 3) * vector_size + k] * B[(j + 1) * vector_size + k];
                    result15 += A[(i + 3) * vector_size + k] * B[(j + 2) * vector_size + k];
                    result16 += A[(i + 3) * vector_size + k] * B[(j + 3) * vector_size + k];
                }
                
                C[i * vector_size + j] = result1;
                C[i * vector_size + (j + 1)] = result2;
                C[i * vector_size + (j + 2)] = result3;
                C[i * vector_size + (j + 3)] = result4;

                C[(i + 1) * vector_size + j] = result5;
                C[(i + 1) * vector_size + (j + 1)] = result6;
                C[(i + 1) * vector_size + (j + 2)] = result7;
                C[(i + 1) * vector_size + (j + 3)] = result8;

                C[(i + 2) * vector_size + j] = result9;
                C[(i + 2) * vector_size + (j + 1)] = result10;
                C[(i + 2) * vector_size + (j + 2)] = result11;
                C[(i + 2) * vector_size + (j + 3)] = result12;

                C[(i + 3) * vector_size + j] = result13;
                C[(i + 3) * vector_size + (j + 1)] = result14;
                C[(i + 3) * vector_size + (j + 2)] = result15;
                C[(i + 3) * vector_size + (j + 3)] = result16;
            }
        }
    };

    ////////////////////////////////////////////////////////////////////////
    // [1] Original code implemantion with HPX
    std::size_t t_origin = hpx::util::high_resolution_clock::now();

    hpx::parallel::for_each(hpx::parallel::par, time_range.begin(), time_range.end(), f);

    std::size_t elapsed_origin = hpx::util::high_resolution_clock::now() - t_origin;

    ////////////////////////////////////////////////////////////////////////
    // [2] Efficient chunk size

    std::size_t t_chunk = hpx::util::high_resolution_clock::now();

    hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::adaptive_chunk_size()), 
        time_range.begin(), time_range.end(), f);

    std::size_t elapsed_chunk = hpx::util::high_resolution_clock::now() - t_chunk;

    ////////////////////////////////////////////////////////////////////////
    // [3] Prefetching:
    std::size_t pref_dist_fac = 2;
    auto policy = hpx::parallel::par;

    std::size_t t_prefetch = hpx::util::high_resolution_clock::now();

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, pref_dist_fac, A, B, C), 
        time_range.begin(), time_range.end(), f);

    std::size_t elapsed_prefetch = hpx::util::high_resolution_clock::now() - t_prefetch;

    ////////////////////////////////////////////////////////////////////////
    // Printing results
    std::cout << std::left << std::setw(6) << std::setfill(' ') << "time_origin = ";
    std::cout << std::left << std::setw(6) << std::setfill(' ') << elapsed_origin;
    std::cout << std::left << std::setw(6) << std::setfill(' ') << "\ntime_chunk = ";
    std::cout << std::left << std::setw(6) << std::setfill(' ') << elapsed_chunk;
    std::cout << std::left << std::setw(6) << std::setfill(' ') << "\ntime_prefetching = ";
    std::cout << std::left << std::setw(6) << std::setfill(' ') << elapsed_prefetch << std::endl;
}

int hpx_main(int argc, char* argv[])
{
    // Initialization
    std::size_t size_of_mat = vector_size * vector_size;
    std::vector<double> A(size_of_mat, 0.0);
    std::vector<double> B(size_of_mat, 0.0);
    std::vector<double> C(size_of_mat, 0.0);

    for(std::size_t r = 0; r < size_of_mat; r++) {
        srand(r);
        A[r] = (r + 1) * (r + 10) * (rand() % 10 + 1);
        B[r] = (r + 1) * (r + 10) * (rand() % 10 + 1);
        C[r] = (r + 1) * (r + 10) * (rand() % 10 + 1);
    }

    comparing_perfromances(A, B, C);
  
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
