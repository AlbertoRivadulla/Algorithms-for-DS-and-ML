#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <iostream>
#include <vector>

// Function to find the indices of the largest n elements in a vector
std::vector<int> indicesOfLargestVals( std::vector<double> vec, int nIndices );

// Map the indices (i, j) of a square matrix to indices of an array that represents
// sequentially the elements of the upper triangle of the original one
inline int mapIndicesUpperTriang( const int& i, const int& j, const int& n )
{
    return j - 1 - (i*i + 3*i) / 2 + i*n;
}

// Map the indices (i, j) of a square matrix to indices of an array that represents
// sequentially the elements of the upper triangle of the original one.
// This also works for the lower triangle of a symmetric matrix.
inline int mapIndicesUpperTriangSym( const int& i, const int& j, const int& n )
{
    return (i < j) ? j - 1 - (i*i + 3*i) / 2 + i*n
                   : i - 1 - (j*j + 3*j) / 2 + j*n ;
}

#endif
