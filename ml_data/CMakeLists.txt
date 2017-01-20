#################################
#Compile nast_hpx
#@date 2016-04
#@author troska@ins.uni-bonn.de
#################################

#HPX requires this version
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#C++-based project
project(sga CXX)

find_package(HPX)

add_compile_options(-std=c++14)

add_hpx_executable(
    main
    ESSENTIAL
    SOURCES main.cpp
)

