#ifndef MYTESTLIB_H
#define MYTESTLIB_H

#include <iostream>

/* 
   Based on:
    https://nesi.github.io/perf-training/python-scatter/ctypes
    https://stackoverflow.com/questions/145270/calling-c-c-from-python
*/

// Class with simple methods
class Foo 
{
    public:
        void bar();
};

void testPrint();

// Wrappers for the class above, to be used by ctypes
extern "C"
{
    //////////////////////////////////////////////////////
    // Wrappers without pointers

    // // Function to create a new instance of Foo
    // Foo Foo_new();
    //
    // // Function to run the bar method of Foo
    // void Foo_bar( Foo foo );

    //////////////////////////////////////////////////////
    // Wrappers with a pointer to the class above

    // Function to create a new instance of Foo
    Foo* Foo_new();

    // Function to run the bar method of Foo
    void Foo_bar( Foo* foo );

    //////////////////////////////////////////////////////
    // Other test functions

    // Function that sums all the elements in an array of integers
    // The return type is a 64-bit integer
    long long sum_array( int n, int* my_array );
}

#endif
