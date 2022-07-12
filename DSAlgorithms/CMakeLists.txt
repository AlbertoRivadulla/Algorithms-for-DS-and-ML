################################################################################
#
# Compile the C++ source to a shared library
# Based on: https://stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake
#
################################################################################
#
# After running 
#     cmake ..
# inside the directory /build, run also
#     make 
#     make install
# This will copy the files .so and the PUBLIC_HEADER to the directory specified in
# the install statement (in this case PROJECT_SOURCE_DIR)
#
################################################################################

cmake_minimum_required(VERSION 3.15)

# Obtain a file compile_commands.json used by ccls (through the plugin coc.nvim)
# to provide code completion in Neovim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

# Name of the project
project("DSAlgorithms" 
    VERSION 1.0.1 
    DESCRIPTION "Library with common algorithms for Machine Learning and Data Science.")

set(CMAKE_BUILD_TYPE Release)

# Declare a new library target, with the source files
add_library(DSAlgorithms SHARED
    ${PROJECT_SOURCE_DIR}/src/linalg.cpp
    ${PROJECT_SOURCE_DIR}/src/clustering.cpp
    ${PROJECT_SOURCE_DIR}/src/dimensionalReduction.cpp
    )

# Set the VERSION property
# Optional, but a good practice
set_target_properties( DSAlgorithms PROPERTIES VERSION ${PROJECT_VERSION} )

# Also possible to do the following, so libDSAlgorithms.so.1 is a symlink to libDSAlgorithms.so.1.0.1
# set_target_properties( DSAlgorithms PROPERTIES SOVERSION 1 )

# Declare a public API header of the library, 
# It can be more than one header, separating them with ;
set_target_properties( DSAlgorithms PROPERTIES PUBLIC_HEADER 
    "${PROJECT_SOURCE_DIR}/src/DSAlgorithms.h;${PROJECT_SOURCE_DIR}/src/clustering.h"
    )

# Add the include directories for the target
target_include_directories( DSAlgorithms PRIVATE 
    ${PROJECT_SOURCE_DIR}/src
    )

# The following copies the .so files and the PUBLIC_HEADER to PROJECT_SOURCE_DIR.
install(TARGETS DSAlgorithms 
        DESTINATION ${PROJECT_SOURCE_DIR}
       )

# The following copies the .so files and the public headers in different directories.
# 'include' is the system's include directory, by default (in Linux)
#     /usr/local/include
# install(TARGETS DSAlgorithms 
#         LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}
#         # PUBLIC_HEADER DESTINATION include
#        )

# # Get rid of the cmake_install.cmake file created
# set(CMAKE_SKIP_INSTALL_RULES True)

# Move the file compile_commands.json to the directory of the project.
execute_process(COMMAND cp compile_commands.json ..) 
