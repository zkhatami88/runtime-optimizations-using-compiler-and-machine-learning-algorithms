#HPX requires this version
cmake_minimum_required(VERSION 2.8.4)
project(main CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(HPX_ROOT "/home/zahra/Projects/HPX/repo/hpx")
set(HPX_DIR "/home/zahra/Projects/HPX/build/lib/cmake/HPX")

find_package(HPX)

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

set(SOURCE_FILES main.cpp)
add_hpx_executable(
    main
    ESSENTIAL
    SOURCES ${SOURCE_FILES}
    COMPONENT_DEPENDENCIES iostreams)
