# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/erik/downloads/rand/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/erik/downloads/rand/lemon-1.3.1/build

# Include any dependencies generated for this target.
include test/CMakeFiles/error_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/error_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/error_test.dir/flags.make

test/CMakeFiles/error_test.dir/error_test.cc.o: test/CMakeFiles/error_test.dir/flags.make
test/CMakeFiles/error_test.dir/error_test.cc.o: ../test/error_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/erik/downloads/rand/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/error_test.dir/error_test.cc.o"
	cd /home/erik/downloads/rand/lemon-1.3.1/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/error_test.dir/error_test.cc.o -c /home/erik/downloads/rand/lemon-1.3.1/test/error_test.cc

test/CMakeFiles/error_test.dir/error_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/error_test.dir/error_test.cc.i"
	cd /home/erik/downloads/rand/lemon-1.3.1/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erik/downloads/rand/lemon-1.3.1/test/error_test.cc > CMakeFiles/error_test.dir/error_test.cc.i

test/CMakeFiles/error_test.dir/error_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/error_test.dir/error_test.cc.s"
	cd /home/erik/downloads/rand/lemon-1.3.1/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erik/downloads/rand/lemon-1.3.1/test/error_test.cc -o CMakeFiles/error_test.dir/error_test.cc.s

test/CMakeFiles/error_test.dir/error_test.cc.o.requires:

.PHONY : test/CMakeFiles/error_test.dir/error_test.cc.o.requires

test/CMakeFiles/error_test.dir/error_test.cc.o.provides: test/CMakeFiles/error_test.dir/error_test.cc.o.requires
	$(MAKE) -f test/CMakeFiles/error_test.dir/build.make test/CMakeFiles/error_test.dir/error_test.cc.o.provides.build
.PHONY : test/CMakeFiles/error_test.dir/error_test.cc.o.provides

test/CMakeFiles/error_test.dir/error_test.cc.o.provides.build: test/CMakeFiles/error_test.dir/error_test.cc.o


# Object files for target error_test
error_test_OBJECTS = \
"CMakeFiles/error_test.dir/error_test.cc.o"

# External object files for target error_test
error_test_EXTERNAL_OBJECTS =

test/error_test: test/CMakeFiles/error_test.dir/error_test.cc.o
test/error_test: test/CMakeFiles/error_test.dir/build.make
test/error_test: lemon/libemon.a
test/error_test: /lib/libglpk.so
test/error_test: test/CMakeFiles/error_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erik/downloads/rand/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable error_test"
	cd /home/erik/downloads/rand/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/error_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/error_test.dir/build: test/error_test

.PHONY : test/CMakeFiles/error_test.dir/build

test/CMakeFiles/error_test.dir/requires: test/CMakeFiles/error_test.dir/error_test.cc.o.requires

.PHONY : test/CMakeFiles/error_test.dir/requires

test/CMakeFiles/error_test.dir/clean:
	cd /home/erik/downloads/rand/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -P CMakeFiles/error_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/error_test.dir/clean

test/CMakeFiles/error_test.dir/depend:
	cd /home/erik/downloads/rand/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erik/downloads/rand/lemon-1.3.1 /home/erik/downloads/rand/lemon-1.3.1/test /home/erik/downloads/rand/lemon-1.3.1/build /home/erik/downloads/rand/lemon-1.3.1/build/test /home/erik/downloads/rand/lemon-1.3.1/build/test/CMakeFiles/error_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/error_test.dir/depend

