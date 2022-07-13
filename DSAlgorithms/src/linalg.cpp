#include "linalg.h"

//=============================================================================
// Auxiliary functions
//=============================================================================

// Print a vector
template <typename T>
void printVector( const T* v, const int& n )
{
    std::cout << "( ";
    for (int i = 0; i < n; ++i)
        std::cout << v[i] << ' ';
    std::cout << ")\n";
}
template <typename T>
void printVector( const std::vector<T>& v )
{
    std::cout << "( ";
    for (int i = 0; i < v.size(); ++i)
        std::cout << v[i] << ' ';
    std::cout << ")\n";
}

// Print a square matrix
template <typename T>
void printSqMatrix( const T* M, const int& n )
{
    for (int j = 0; j < n; ++j)
    {
        std::cout << "[ ";
        for (int i = 0; i < n; ++i)
            std::cout << M[j*n + i] << ' ';
        std::cout << "]\n";
    }
}

// Check if all components of two vectors are close, within a given tolerance
bool allCloseVectors( const std::vector<double>& v1, 
                      const std::vector<double>& v2, const double& tolerance )
{
    for (int i = 0; i < v1.size(); ++i)
        if ( std::fabs( v1[i] - v2[i] ) > tolerance )
            return false;
    return true;
}

// Check if all components of two matrices are close, within a given tolerance
template <typename T>
bool allCloseSqMatrices( const T* M1, const T* M2, const int& n, 
                         const double& tolerance)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if ( std::fabs( M1[i*n+j] - M2[i*n+j] ) > tolerance )
                return false;
    return true;
}

// Norm of a vector
double normVector( const std::vector<double>& v )
{
    double norm = 0;
    for (int i = 0; i < v.size(); ++i)
        norm += v[i]*v[i];
    return std::sqrt( norm );
}

// Normalize a vector
void normalizeVector( std::vector<double>& v )
{
    // Compute the norm
    double norm = normVector( v );
    // Normalize all the components of the vector
    for (int i = 0; i < v.size(); ++i)
        v[i] /= norm;
}

// Multiply two vectors
template <typename T>
T dotProductVectors( const std::vector<T>& v1, const std::vector<T>& v2 )
{
    T result = 0;
    for (int i = 0; i < v1.size(); ++i)
        result += v1[i] * v2[i];
    return result;
}

template <typename T>
T dotProductVectors( const T* v1, const T* v2, const int& n )
{
    T result = 0;
    for (int i = 0; i < n; ++i)
        result += v1[i] * v2[i];
    return result;
}

// Multiply a square matrix (of dimension nxn) by a vector
template <typename T>
std::vector<T> multiplySqMatrixVector( const T* M, const std::vector<T>& V,
                                       const int& n )
{
    // Initialize the result vector to zeros
    std::vector<T> result(n, 0);
    // Compute the product
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            result[i] += M[i*n + j] * V[j];
    }

    return result;
}

// Multiply two square matrices
template <typename T>
void multiplySqMatrices( const T* M1, const T* M2, T* R, const int& n )
{
    // Compute the product
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            R[ i*n + j ] = 0;
            for (int k = 0; k < n; ++k)
                R[ i*n + j ] += M1[ i*n + k ] * M2[ k*n + j ];
        }
    }
}

// QR decomposition of a matrix using the Gram-Schmidt process
// Implemented following:
//        https://en.wikipedia.org/wiki/QR_decomposition
// Parameters:
//    M: matrix to decompose (a square matrix)
//    R, Q: output
//    n: size of the side of the matrix
template <typename T>
void computeQRDecomposition( const T* A, T* Q, T* R, const int& n )
{
    // Initialize the matrix U
    T* U = new T[ n*n ];

    // Compute the components of the matrices U and Q
    for (int k = 0; k < n; ++k)
    {
        // Compute the vector a_k
        T* ak = new T[ n ];
        for (int i = 0; i < n; ++i)
            ak[i] = A[i*n + k];
        // Compute the elements U_{ik} (or the elements of the vector u_k)
        for (int i = 0; i < n; ++i)
        {
            U[i*n + k] = ak[i];
            // The following does not apply for the first column
            if (k > 0)
            {
                for (int j = 0; j < k; ++j)
                {
                    // Compute the vector u_j
                    T* uj = new T[ n ];
                    for (int l = 0; l < n; ++l)
                        uj[l] = U[ l*n + j ];
                    // Add to this component of the matrix U
                    U[i*n + k] -= dotProductVectors( uj, ak, n ) * uj[i] / dotProductVectors( uj, uj, n );
                }
            }
        }

        // Compute the vector of the matrix Q
        double ukNorm = 0;
        for (int i = 0; i < n; ++i)
            ukNorm += U[ i*n + k ] * U[ i*n + k ];
        ukNorm = std::sqrt(ukNorm);
        for (int i = 0; i < n; ++i)
            Q[i*n + k] = U[i*n + k] / ukNorm;
    }

    // Compute the components of R
    // This is given by
    //      R = Q^T A
    // Since Q has to be transposed, I compute it explicitly here
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            R[ i*n + j ] = 0.;
            for (int k = 0; k < n; ++k)
                // The transposition of Q is manifest in the location of the indices here,
                // since
                //      R_{ij} = Q^T_{ik} A_{kj} = Q_{ki} A_{kj}
                R[ i*n + j ] += Q[ k*n + i ] * A[ k*n + j ];
        }
    }
}

//=============================================================================
// Functions for computing eigenvalues and eigenvectors
//=============================================================================

// Function to compute the eigenvalues and eigenvectors of a hermitian matrix M,
// using the Power Iteration algorithm
// This tends to produce always the eigenvector with the largest eigenvalue
void computeEigenPower( const double* matrix, const int& n, std::vector<double>& eigenVals,
                        std::vector<std::vector<double>>& eigenVectors, double tolerance)
{
    // Set the seed for the random number generator
    std::srand(std::time(0));

    // Random initial vector
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i)
        v[i] = (double)std::rand() / RAND_MAX;
    // Normalize it
    normalizeVector( v );

    // Vector to remember the previous values of v
    std::vector<double> vPrevious = v;

    // Iterate until the defined precision is met
    do
    {
        // Copy the elements of v
        vPrevious = v;
        // Update v and normalize it
        v = multiplySqMatrixVector( matrix, v, n );
        normalizeVector( v );
    }
    while ( !allCloseVectors( v, vPrevious, tolerance ) );

    // Multiply the eigenvector by the matrix to obtain the eigenvalue
    std::vector<double> Mtimesv = multiplySqMatrixVector( matrix, v, n );
    
    // Save the eigenvector and eigenvalue
    eigenVectors[0] = v;
    eigenVals[0] = normVector( Mtimesv );
}

// Function to compute the eigenvalues and eigenvectors of a hermitian matrix M,
// using the QR algorithm
void computeEigenQR( const double* matrix, const int& n, std::vector<double>& eigenVals,
                     std::vector<std::vector<double>>& eigenVectors, double tolerance )
{
    // Matrix to hold the value of the products of the Q matrices obtained
    double* QProduct = new double[ n*n ];
    double* QProductCopy = new double[ n*n ];

    // Identity matrix, needed to check the convergence
    double* identity = new double[ n*n ];
    for (int i = 0; i < n*n; ++i)
        identity[ i ] = 0.;
    for (int i = 0; i < n; ++i)
        identity[ i*n + i ] = 1.;

    // Dummy matrix to compute the products of the original matrix and Q 
    // computed in each iteration
    double* X = new double[ n*n ];

    // QR decomposition of the matrix
    double* Q = new double[ n*n ];
    double* R = new double[ n*n ];
    computeQRDecomposition( matrix, Q, R, n );

    // Initialize the product of Q matrices with the first one
    for (int i = 0; i < n*n; ++i)
        QProduct[ i ] = Q[ i ];

    // Iterate until the desired accuracy is reached
    int count = 0;
    int countMax = 500;
    do
    {
        // Copy the previous components of the matrix of the product
        for (int i = 0; i < n*n; ++i)
            QProductCopy[ i ] = QProduct[ i ];

        // Update the dummy matrix X 
        //      X = R * Q
        multiplySqMatrices( R, Q, X, n );
        // QR-decompose this new matrix
        //      X = Q * R
        computeQRDecomposition( X, Q, R, n );

        // Update the product of Q matrices
        multiplySqMatrices( QProductCopy, Q, QProduct, n );

        // Count the iterations
        ++count;
    }
    while ( !allCloseSqMatrices( identity, Q, n, tolerance ) && count < countMax );

    // If the iteration stopped because we reached the maximum amount of iterations,
    // print a message
    if (count >= countMax)
        std::cout << "The computation of the eigenvectors did not converge!\n";

    // Get the eigenvectors from the columns of the product of matrices Q
    // From these, compute the eigenvalues
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            eigenVectors[i][j] = QProduct[ j*n + i ];

        // Compute the eigenvalue from this
        eigenVals[i] = normVector( multiplySqMatrixVector( matrix, eigenVectors[i], n ) );
    }
}


