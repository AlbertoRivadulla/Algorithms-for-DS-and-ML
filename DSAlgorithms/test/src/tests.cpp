#include "tests.h"

void testComputeEigenvectors()
{
    std::cout << "\n\n==========================================\n";
    std::cout << "Tests of the methods that compute eigenfunctions and eigenvalues\n";
    std::cout << "==========================================\n\n";

    // Input matrix
    double matrix[] { 2.92, 0.86, -1.15,
                      0.86, 6.51, 3.32,
                      -1.15, 3.32, 4.57 };
    // Size of the matrix
    int n = 3;
    // Vectors for the output
    std::vector<double> eigenvalues(n);
    std::vector<std::vector<double>> eigenvectors(n, std::vector<double>(n));

    // Compute the eigenvalues and eigenvectors
    computeEigenPower( matrix, n, eigenvalues, eigenvectors );
    std::cout << "\n\n\t============\n\n";
    std::cout << "Eigenvalues computed with the Power method:\n";
    std::cout << "(this only produces the first eigenvalue/vector)\n\n";
    printEigen(eigenvalues, eigenvectors, n);

    std::cout << "\n\n\t============\n\n";
    computeEigenQR( matrix, n, eigenvalues, eigenvectors );
    std::cout << "Eigenvalues computed with the QR method:\n\n";
    printEigen(eigenvalues, eigenvectors, n);

    std::cout << "\n\n\t============\n\n";
}
