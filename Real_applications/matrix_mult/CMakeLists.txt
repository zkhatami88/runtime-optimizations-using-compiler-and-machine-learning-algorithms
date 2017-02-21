#################################
#Compile nast_hpx
#@date 2016-04
#@author troska@ins.uni-bonn.de
#################################

#HPX requires this version
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#C++-based project
project(sga CXX)

set(HPX_ROOT "/home/zahra/Projects/HPX/repo/hpx")
set(HPX_DIR "/home/zahra/Projects/HPX/build/lib/cmake/HPX")

find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

add_compile_options(-std=c++11)

add_hpx_executable(
    main
    ESSENTIAL
    SOURCES main.cpp
)

