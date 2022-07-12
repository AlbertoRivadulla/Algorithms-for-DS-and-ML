#include "mytestlib_2.h"

//////////////////////////////////////////////////////
// Other test functions

// Function that returns an array with the first n integers
int* first_integers( int n )
{
    std::cout << "computing the first " << n << "integers...\n";
    int* integers = new int[n];
    for (int i = 0; i < n; ++i)
    {
        integers[i] = i;
    }
    return integers;
}
