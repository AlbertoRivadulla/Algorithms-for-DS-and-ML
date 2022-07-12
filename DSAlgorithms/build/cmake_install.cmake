# Install script for directory: /home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms" TYPE SHARED_LIBRARY FILES "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/libDSAlgorithms.so.1.0.1")
  if(EXISTS "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so.1.0.1")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms" TYPE SHARED_LIBRARY FILES "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/libDSAlgorithms.so")
  if(EXISTS "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/libDSAlgorithms.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/DSAlgorithms.h;/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/clustering.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms" TYPE FILE FILES
    "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/DSAlgorithms.h"
    "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/src/clustering.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/albertors/Datos/Data_science_and_machine_learning/Projects/001_Algorithms_in_data_science/DSAlgorithms/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
