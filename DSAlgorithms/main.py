#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ctypes
# from ctypes import cdll

import os
import numpy as np

# Load the library
# The three methods below work for me
# lib = cdll.LoadLibrary( os.getcwd() + "/build/libmylib.so.1.0.1" )
# lib = cdll.LoadLibrary( "./build/libmylib.so.1.0.1" )
lib = ctypes.cdll.LoadLibrary( "./build/libmylib.so" )

###################################################
# Test the methods that call the class Foo
###################################################

foo_instance = lib.Foo_new()
lib.Foo_bar( foo_instance )
print()

###################################################
# Test the sum_array function
###################################################

# Tell Python the arguments and result types for the function sum_array
# The second argument of sum_array is a pointer to an array of integers
lib.sum_array.restype  = ctypes.c_longlong
lib.sum_array.argtypes = [ ctypes.c_int, np.ctypeslib.ndpointer( dtype=np.int32 ) ]

# Create an array
array = np.arange( 0, 10, 1, dtype=np.int32 )

# Sum the elements of the array in C++
array_sum = lib.sum_array( len(array), array )
print( "The sum of the elements of the array is: {}".format( array_sum ) )
print()

###################################################
# Test the first_integers function
###################################################

# Number of integers that I want
nr = 10

# Tell Python the arguments and result types for the function first_integers
lib.first_integers.restype  = ctypes.POINTER( ctypes.c_int * nr )
lib.first_integers.argtypes = [ ctypes.c_int ]

# Get the first integers, and print them
integers = lib.first_integers( nr )
# Convert this to a list:
integers_list = list(integers.contents)
integers_list = [ i for i in integers.contents ]
print("First {} integers:".format(nr))
print(integers_list)
print()
