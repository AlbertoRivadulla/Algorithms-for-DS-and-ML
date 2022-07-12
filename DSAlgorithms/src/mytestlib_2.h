#ifndef MYTESTLIB_2_H
#define MYTESTLIB_2_H

#include "mytestlib.h"

// Wrappers for the class above, to be used by ctypes
extern "C"
{
    //////////////////////////////////////////////////////
    // Other test functions

    // Function that returns an array with the first n integers
    int* first_integers( int n );
}

#endif
