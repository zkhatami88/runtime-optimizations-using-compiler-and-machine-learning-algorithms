cmake_minimum_required(VERSION 2.8.12)
project(hpxloop CXX)

set(HPX_ROOT "/home/zahra/Projects/HPX/build")
set(HPX_DIR "${HPX_ROOT}/lib/cmake/HPX")
set(BOOST_ROOT "/home/zahra/Projects/boost_1_63_0")
#set(CMAKE_C_COMPILER "/usr/local/bin/clang")
#set(CMAKE_CXX_COMPILER "/usr/local/bin/clang++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(HPX)

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

set(SOURCE_FILES hpxloop.cpp)
add_hpx_executable(hpxloop
    ESSENTIAL
    SOURCES ${SOURCE_FILES}
    COMPONENT_DEPENDENCIES iostreams)
