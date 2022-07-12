#include "mytestlib.h"

/* 
   Based on:
    https://nesi.github.io/perf-training/python-scatter/ctypes
    https://stackoverflow.com/questions/145270/calling-c-c-from-python
*/

void testPrint()
{
    std::cout << "Mensaxe de proba\n";
}

//////////////////////////////////////////////////////
// Method of the class Foo
void Foo::bar()
{
    std::cout << "Called bar method of Foo class in C++.\n";
}

//////////////////////////////////////////////////////
// Wrappers without pointers for the class above

// // Function to create a new instance of Foo
// Foo Foo_new()
// {
//     std::cout << "Created a new instance of Foo in C++.\n";
//     return Foo();
// }
//
// // Function to run the bar method of Foo
// void Foo_bar( Foo foo )
// {
//     foo.bar();
// }

//////////////////////////////////////////////////////
// Wrappers with a pointer to the class above

// Function to create a new instance of Foo
Foo* Foo_new()
{
    std::cout << "Created a new instance of Foo in C++.\n";
    return new Foo();
}

// Function to run the bar method of Foo
void Foo_bar( Foo* foo )
{
    foo->bar();
}

//////////////////////////////////////////////////////
// Other test functions

// Function that sums all the elements in an array of integers
// The return type is a 64-bit integer
long long sum_array( int n, int* my_array )
{
    std::cout << "Summing the array in C++...\n";
    long long sum = 0;
    for ( int i = 0; i < n; ++i )
        sum += my_array[i];
    return sum;
}

