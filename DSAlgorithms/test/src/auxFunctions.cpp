#include "auxFunctions.h"

void printEigen( std::vector<double> eigenvalues, 
                 std::vector<std::vector<double>> eigenvectors, int& n)
{
    std::cout << "Eigenvalues: ( ";
    for (int i = 0; i < n-1; ++i)
    {
        std::cout << eigenvalues[i] << ", ";
    }
    std::cout << eigenvalues[n-1] << " )\n\n";

    std::cout << "Eigenvectors:\n";
    for (int i = 0; i < n; ++i)
    {
        std::cout << "( ";
        for (int j = 0; j < n-1; ++j)
        {
            std::cout << eigenvectors[i][j] << ", ";
        }
        std::cout << eigenvectors[i][n-1] << " )\n";
    }
}

