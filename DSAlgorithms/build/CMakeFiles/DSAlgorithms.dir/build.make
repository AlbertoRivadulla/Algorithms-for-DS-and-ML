# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build

# Include any dependencies generated for this target.
include CMakeFiles/DSAlgorithms.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/DSAlgorithms.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/DSAlgorithms.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DSAlgorithms.dir/flags.make

CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o: CMakeFiles/DSAlgorithms.dir/flags.make
CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o: ../src/linalg.cpp
CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o: CMakeFiles/DSAlgorithms.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o -MF CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o.d -o CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o -c /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/linalg.cpp

CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/linalg.cpp > CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.i

CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/linalg.cpp -o CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.s

CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o: CMakeFiles/DSAlgorithms.dir/flags.make
CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o: ../src/clustering.cpp
CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o: CMakeFiles/DSAlgorithms.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o -MF CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o.d -o CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o -c /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/clustering.cpp

CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/clustering.cpp > CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.i

CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/clustering.cpp -o CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.s

CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o: CMakeFiles/DSAlgorithms.dir/flags.make
CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o: ../src/dimensionalReduction.cpp
CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o: CMakeFiles/DSAlgorithms.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o -MF CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o.d -o CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o -c /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/dimensionalReduction.cpp

CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/dimensionalReduction.cpp > CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.i

CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/dimensionalReduction.cpp -o CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.s

# Object files for target DSAlgorithms
DSAlgorithms_OBJECTS = \
"CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o" \
"CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o" \
"CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o"

# External object files for target DSAlgorithms
DSAlgorithms_EXTERNAL_OBJECTS =

libDSAlgorithms.so.1.0.1: CMakeFiles/DSAlgorithms.dir/src/linalg.cpp.o
libDSAlgorithms.so.1.0.1: CMakeFiles/DSAlgorithms.dir/src/clustering.cpp.o
libDSAlgorithms.so.1.0.1: CMakeFiles/DSAlgorithms.dir/src/dimensionalReduction.cpp.o
libDSAlgorithms.so.1.0.1: CMakeFiles/DSAlgorithms.dir/build.make
libDSAlgorithms.so.1.0.1: CMakeFiles/DSAlgorithms.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library libDSAlgorithms.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DSAlgorithms.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library libDSAlgorithms.so.1.0.1 libDSAlgorithms.so.1.0.1 libDSAlgorithms.so

libDSAlgorithms.so: libDSAlgorithms.so.1.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate libDSAlgorithms.so

# Rule to build all files generated by this target.
CMakeFiles/DSAlgorithms.dir/build: libDSAlgorithms.so
.PHONY : CMakeFiles/DSAlgorithms.dir/build

CMakeFiles/DSAlgorithms.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DSAlgorithms.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DSAlgorithms.dir/clean

CMakeFiles/DSAlgorithms.dir/depend:
	cd /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/CMakeFiles/DSAlgorithms.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DSAlgorithms.dir/depend

