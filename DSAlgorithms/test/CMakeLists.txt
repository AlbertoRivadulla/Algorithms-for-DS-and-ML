################################################################################
#
# At the end of this file there is another one that should be used if I wanted
# to import the library as it should be distributed (this is, only using its public
# header).
#
################################################################################
#
# The following code allows me to access all the components (headers) of the library,
# without need to install it.
#
################################################################################

cmake_minimum_required(VERSION 3.7)

# Obtain a file compile_commands.json used by ccls (through the plugin coc.nvim)
# to provide code completion in Neovim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

# Name of the project
project(test)

# Your-external "mylib", add GLOBAL if the imported library is located in directories above the current.
# Add the library DSAlgorithms
add_library( DSAlgorithms SHARED IMPORTED GLOBAL )
# # set_target_properties( DSAlgorithms PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR} )
# Directory of the library
set_target_properties( DSAlgorithms PROPERTIES 
    IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/../build/libDSAlgorithms.so 
    )

# Create a variable with all the include directories
set(INCLUDE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src

    # Include directory of the linked library
    ${PROJECT_SOURCE_DIR}/../src
   )
include_directories(${INCLUDE})

# Create a variable with a link to all cpp files to compile
set(SOURCES
    main.cpp
    # ${PROJECT_SOURCE_DIR}/src/rungeKutta.cpp
   )

add_executable(main ${SOURCES})

# Link to the library
target_link_libraries(main DSAlgorithms)

# Get rid of the cmake_install.cmake file created
set(CMAKE_SKIP_INSTALL_RULES True)

# Move the file compile_commands.json to the directory of the project.
execute_process(COMMAND cp compile_commands.json ..) 

################################################################################
#
# The next code is an example of how to import the library as it is distributed
# (only using its public header).
# Inside the code I need to include the public header of the library
#
################################################################################

# cmake_minimum_required(VERSION 3.7)
#
# # Obtain a file compile_commands.json used by ccls (through the plugin coc.nvim)
# # to provide code completion in Neovim
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
#
# # Name of the project
# project(test)
#
# # Directory of the library
# set( DSALGORITHMS_DIR "${PROJECT_SOURCE_DIR}/.." )
#
# # Find the library.
# # From now on, I refer to it as DSALGORITHMS_LIB
# find_library(DSALGORITHMS_LIB
#     NAMES DSAlgorithms
#     # The directory of the .so and public header files
#     HINTS ${DSALGORITHMS_DIR}
# )
# # Print the directory of the library
# message(STATUS ${DSALGORITHMS_LIB})
#
# # Create a variable with all the include directories
# set(INCLUDE
#     ${PROJECT_SOURCE_DIR}
#     ${PROJECT_SOURCE_DIR}/src
#
#     # Include directory of the linked library
#     ${DSALGORITHMS_DIR}
#    )
# include_directories(${INCLUDE})
#
# # Create a variable with a link to all cpp files to compile
# set(SOURCES
#     main.cpp
#     # ${PROJECT_SOURCE_DIR}/src/rungeKutta.cpp
#    )
# add_executable(main ${SOURCES})
#
# # Link to the library
# target_link_libraries(main PRIVATE ${DSALGORITHMS_LIB})
#
# # Get rid of the cmake_install.cmake file created
# set(CMAKE_SKIP_INSTALL_RULES True)
#
# # Move the file compile_commands.json to the directory of the project.
# execute_process(COMMAND cp compile_commands.json ..) 
