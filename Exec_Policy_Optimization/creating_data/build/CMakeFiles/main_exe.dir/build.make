# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/zahra/.local/bin/cmake

# The command to remove a file.
RM = /home/zahra/.local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build

# Include any dependencies generated for this target.
include CMakeFiles/main_exe.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main_exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main_exe.dir/flags.make

CMakeFiles/main_exe.dir/main.cpp.o: CMakeFiles/main_exe.dir/flags.make
CMakeFiles/main_exe.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main_exe.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main_exe.dir/main.cpp.o -c /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/main.cpp

CMakeFiles/main_exe.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main_exe.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/main.cpp > CMakeFiles/main_exe.dir/main.cpp.i

CMakeFiles/main_exe.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main_exe.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/main.cpp -o CMakeFiles/main_exe.dir/main.cpp.s

CMakeFiles/main_exe.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/main_exe.dir/main.cpp.o.requires

CMakeFiles/main_exe.dir/main.cpp.o.provides: CMakeFiles/main_exe.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/main_exe.dir/build.make CMakeFiles/main_exe.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/main_exe.dir/main.cpp.o.provides

CMakeFiles/main_exe.dir/main.cpp.o.provides.build: CMakeFiles/main_exe.dir/main.cpp.o


# Object files for target main_exe
main_exe_OBJECTS = \
"CMakeFiles/main_exe.dir/main.cpp.o"

# External object files for target main_exe
main_exe_EXTERNAL_OBJECTS =

main: CMakeFiles/main_exe.dir/main.cpp.o
main: CMakeFiles/main_exe.dir/build.make
main: /home/zahra/Projects/HPX/build/lib/libhpx_init.a
main: /home/zahra/Projects/HPX/build/lib/libhpx.so.1.0.0
main: /home/zahra/Projects/HPX/build/lib/libhpx_init.a
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_chrono.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_date_time.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_filesystem.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_program_options.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_regex.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_system.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_thread.so
main: /usr/lib/x86_64-linux-gnu/libpthread.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_context.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_random.so
main: /home/zahra/Projects/boost_1_63_0/stage/lib/libboost_atomic.so
main: /usr/lib/libtcmalloc_minimal.so
main: /usr/lib/x86_64-linux-gnu/libhwloc.so
main: CMakeFiles/main_exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main_exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main_exe.dir/build: main

.PHONY : CMakeFiles/main_exe.dir/build

CMakeFiles/main_exe.dir/requires: CMakeFiles/main_exe.dir/main.cpp.o.requires

.PHONY : CMakeFiles/main_exe.dir/requires

CMakeFiles/main_exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main_exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main_exe.dir/clean

CMakeFiles/main_exe.dir/depend:
	cd /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build /home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/build/CMakeFiles/main_exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main_exe.dir/depend

