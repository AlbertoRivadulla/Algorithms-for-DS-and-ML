#ifndef LINALG_H
#define LINALG_H

#include <iostream>
#include <vector>
#include <random>
#include <ctime>

//=============================================================================
// Auxiliary functions
//=============================================================================

// Print a vector
template <typename T>
void printVector( const T* v, const int& n );
template <typename T>
void printVector( const std::vector<T>& v );
void printVector( const std::vector<double>& v );

// Print a square matrix
template <typename T>
void printSqMatrix( const T* M, const int& n );
void printSqMatrix( const double* M, const int& n );

// Print the eigenvalues and eigenvectors
void printEigen( const std::vector<double> eigenvalues, 
                 const std::vector<std::vector<double>> eigenvectors, const int& n);

// Check if all components of two vectors are close, within a given tolerance
bool allCloseVectors( const std::vector<double>& v1, 
                      const std::vector<double>& v2, const double& tolerance );

// Check if all components of two matrices are close, within a given tolerance
template <typename T>
bool allCloseSqMatrices( const T* M1, const T* M2, const int& n, 
                         const double& tolerance); 

// Norm of a vector
double normVector( const std::vector<double>& v );

// Multiply two vectors
template <typename T>
T dotProductVectors( const std::vector<T>& v1, const std::vector<T>& v2 );
template <typename T>
T dotProductVectors( const T* v1, const T* v2, const int& n );

// Normalize a vector
void normalizeVector( std::vector<double>& v );

// Multiply a square matrix (of dimension nxn) by a vector
template <typename T>
std::vector<T> multiplySqMatrixVector( const T* M, const std::vector<T>& V,
                                       const int& n );

// Multiply two square matrices
template <typename T>
void multiplySqMatrices( const T* M1, const T* M2, T* R, const int& n );
// void multiplySqMatrices( const double* M1, const double* M2, double* R, const int& n );

// QR decomposition of a matrix using the Gram-Schmidt process
// Implemented following https://en.wikipedia.org/wiki/QR_decomposition
// Parameters:
//    M: matrix to decompose (a square matrix)
//    R, Q: output
//    n: size of the side of the matrix
template <typename T>
void computeQRDecomposition( const T* A, T* Q, T* R, const int& n );


//=============================================================================
// Functions for computing eigenvalues and eigenvectors
//=============================================================================

// Function to compute the eigenvalues and eigenvectors of a hermitian matrix M,
// using the Power Iteration algorithm
// This tends to produce always the eigenvector with the largest eigenvalue
void computeEigenPower( const double* matrix, const int& n, std::vector<double>& eigenVals,
                        std::vector<std::vector<double>>& eigenVectors, double tolerance=0.000001 );

// Function to compute the eigenvalues and eigenvectors of a hermitian and positive
// definite matrix M, using the QR algorithm
void computeEigenQR( const double* matrix, const int& n, std::vector<double>& eigenVals,
                     std::vector<std::vector<double>>& eigenVectors, double tolerance=0.0001 );

#endif
