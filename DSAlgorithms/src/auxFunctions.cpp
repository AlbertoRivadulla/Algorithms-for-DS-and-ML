#include "auxFunctions.h"

// Function to find the indices of the largest n elements in a vector
std::vector<int> indicesOfLargestVals( std::vector<double> vec, int nIndices )
{
    // Lambda function to check if a number is present in a list
    auto nrInList = [] ( const int& in, const std::vector<int>& list ) -> bool
    {
        for ( auto el : list )
            if ( el == in )
                return true;
        return false;
    };

    // Initialize the vector with the indices
    std::vector<int> indicesLargest ( nIndices, -1 );

    double maxVal;
    for (int i = 0; i < nIndices; ++i)
    {
        maxVal = 0.;
        for (int j = 0; j < vec.size(); ++j)
        {
            // for ( auto in : indicesLargest )
            // {
            //     if ( j == in )
            //         continue;
            // }
            if ( nrInList( j, indicesLargest ) )
                continue;
            if ( vec[j] > maxVal )
            {
                maxVal = vec[j];
                indicesLargest[i] = j;
            }
        }
    }

    return indicesLargest;
}


