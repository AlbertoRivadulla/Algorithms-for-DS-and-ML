#include "dimensionalReduction.h"

//------------------------------------------------------------------------------
// Principal component analysis

// Principal component analysis of a dataset
//      data: the dataset
//      nData: the number of points
//      dim: the number of features
//      reducedData: dataset with a reduced number of features
//      dimReduced: the number of features after the reduction
double* computePCA( const double* data, const int nData, const int dim, 
                    double* reducedData, const int dimReduced )
{
    std::cout << "Performing PCA on the dataset...\n";

    // Normalize the data and compute the means
    double* dataNorm = new double [ nData * dim ];
    double* means = new double [ dim ];
    normalizeDataset( data, nData, dim, dataNorm, means );

    // Since the dataset is normalized, all the means are equal to zero now
    for (int i = 0; i < dim; ++i) 
        means[i] = 0.;

    std::cout << "Computing the correlation matrix...";
    // Compute the correlation matrix
    double* corrMatrix = new double [ dim * dim ];
    compCorrMatrixDataset( dataNorm, nData, dim, means, corrMatrix );
    std::cout << "\tDone.\n";

    std::cout << "Computing the eigenvectors of the correlation matrix...";
    // Find the eigenvalues and eigenvectors of the correlation matrix
    std::vector<double> eigenValues( dim );
    std::vector<std::vector<double>> eigenVectors( dim, std::vector<double>(dim));
    computeEigenQR( corrMatrix, dim, eigenValues, eigenVectors );
    std::cout << "\tDone.\n";

    // printEigen( eigenValues, eigenVectors, dim );

    // Get the indices of the largest eigenvalues
    std::vector<int> indicesLargestEigen = indicesOfLargestVals( eigenValues, dimReduced );

    std::cout << "Projecting the data on the directions of the largest eigenvalues...";
    // Project the data into the directions with the largest eigenvalues
    for (int i = 0; i < dimReduced; ++i)
    {
        for (int j = 0; j < nData; ++j)
        {
            reducedData[ j*dimReduced + i ] = 0.;
            for (int k = 0; k < dim; ++k)
                reducedData[ j*dimReduced + i ] += eigenVectors[indicesLargestEigen[i]][k] * dataNorm[ j*dim + k ];
        }
    }
    std::cout << "\tDone.\n";

    // Return the most significant eigenvectors
    double* mostSignificantEigenvectors = new double [ dimReduced * dim ];
    for (int i = 0; i < dimReduced; ++i)
    {
        for (int j = 0; j < dim; ++j)
            mostSignificantEigenvectors[ j + i*dim ] = eigenVectors[indicesLargestEigen[i]][j];
    }

    delete[] dataNorm;
    delete[] means;
    delete[] corrMatrix;

    return mostSignificantEigenvectors;
}


//------------------------------------------------------------------------------
// t-SNE

// Compute the components of the matrix p_ij for t-SNE, from an array of data
// points of a given dimension, and with a desired perplexity
// void computeP_ijFortSNE( double* p, const double* data, const int& nData,
//                          const int& dim, const double& perplexity)
void computeP_ijFortSNE( double* p, const double* data, const int& nData,
                         const int& dim, const double& perplexity)
{
    // Tolerance for the value of the local entropy
    double toleranceEntropy = 0.0000001;
    // Desired value of the local entropy.
    // I use natural logarithms instead of base 2 as in my notes.
    double HGoal = std::log( perplexity );

    // Compute the distances between the points in an upper-triangle matrix
    double* distSqX = new double [ (int)( nData * (nData - 1.) / 2. ) ];
    for ( int i = 0; i < nData; ++i )
    {
        for ( int j = i+1; j < nData; ++j )
        {
            int index = mapIndicesUpperTriang( i, j, nData );
            distSqX[ index ] = 0.;
            for ( int k = 0; k < dim; ++k )
                distSqX[ index ] += ( data[ i*dim + k ] - data[ j*dim + k ] ) *
                                    ( data[ i*dim + k ] - data[ j*dim + k ] );
        }
    }

    // Find the maximum value of the squared distance, to initialize the betas
    double maxDistSqX = 0.;
    for ( int i = 0; i < (int)(nData*(nData-1.)/2.); ++i )
        if ( distSqX[ i ] > maxDistSqX )
            maxDistSqX = distSqX[ i ];

    // Initialize the vector for the bandwidth parameters.
    // Instead of the variables sigma_i, I will compute beta_i defined as:
    //      beta_i = 1/2 sigma_i^2
    double* betas = new double [ nData ];

    // Compute the bandwidth parameter for each point x_i
    for ( int i = 0; i < nData; ++i )
    {
        // Initialize the bandwidth parameter
        betas[ i ] = 1. / maxDistSqX;

        // Initialize parameters for the computation
        double H;
        // Initialize the limits for the window of values of beta, to find the
        // root of H - HGoal
        double betaLower = 0.;
        double betaUpper = std::numeric_limits<double>::max();

        // Search for the values of sigma that reproduce the desired perplexity
        int iteration = 0;
        do
        {

            // Compute the denominator of p_{j|i}
            double denominatorP = 0.;
            for ( int k = 0; k < nData; ++k )
                if ( k != i )
                    denominatorP += std::exp( - betas[i] * distSqX[ mapIndicesUpperTriangSym(i,k,nData) ] );

            // Compute the components of the matrix p_{j|i} for different j
            // Use the same loop to compute the local entropy H(p_i)
            H = 0.;
            for ( int j = 0; j < nData; ++j )
            {
                if ( j == i ) 
                    p[ i*nData + j ] = 0.;
                else
                {
                    p[ i*nData + j ] = std::exp( - betas[i] * distSqX[ mapIndicesUpperTriangSym(i,j,nData) ] ) 
                                       / denominatorP;
                    H += - p[ i*nData + j ] * std::log( p[ i*nData + j ] );
                }
            }

            // Update the window of values of beta, to find the root of H - HGoal
            // In the following I use the fact that dH/dbeta < 0 always.
            if ( H - HGoal > 0. )
            {
                // In this case the root is above the current value of beta
                betaLower = betas[ i ];
                if ( betaUpper == std::numeric_limits<double>::max() )
                    betas[ i ] = betas[ i ] * 2.;
                else
                    betas[ i ] = ( betas[ i ] + betaUpper ) / 2.;
            }
            else
            {
                // In this case the root is below the current value of beta
                betaUpper = betas[ i ];
                betas[ i ] = ( betas[ i ] + betaLower ) / 2.;
            }

        }
        while ( std::fabs( H - HGoal ) > toleranceEntropy && iteration++ < 100000 );
    }

    // Delete the previously defined arrays
    delete[] distSqX;
    delete[] betas;
}

// Compute the distances squared between the variables U_i.
void computeDistSqU( const double* U, const int& nData, const int& dimReduced,
                     double* distSqU )
{
    // Iterate through the points.
    // Since the matrix of distances is symmetric and with null diagonal, it is
    // enough to compute the elements in the upper triangle.
    for ( int i = 0; i < nData; ++i )
    {
        for ( int j = i + 1; j < nData; ++j )
        {
            int index = mapIndicesUpperTriang( i, j, nData );
            distSqU[ index ] = 0.;
            for ( int k = 0; k < dimReduced; ++k )
                distSqU[ index ] += ( U[i*dimReduced + k] - U[j*dimReduced + k] ) *
                                    ( U[i*dimReduced + k] - U[j*dimReduced + k] );
        }
    }
}

// Compute the elements of the similar probability distribution q_{ij}
void computeSimilarDistrtSNE( double* q, const double* distSqU, const int& nData )
{
    // Iterate through the points
    for ( int i = 0; i < nData; ++i )
    {
        // Compute the denominator
        double denominator = 0.;
        for ( int k = 0; k < nData; ++k )
        {
            if ( k != i )
                denominator += 1. / ( 1. + distSqU[ mapIndicesUpperTriangSym(i,k,nData) ] );
        }

        // Normalization factor, such that the sum of all components of q_{ij} is 1
        denominator *= nData;

        for ( int j = i + 1; j < nData; ++j )
        {
            int index = mapIndicesUpperTriang( i, j, nData );
            q[ index ] = 1. / ( 1. + distSqU[ index ] ) / denominator;
        }
    }
}

// Compute the cost function as the Kullback-Leibler divergence of the probability 
// distributions q_{ij} and p_{ij}
double computeCosttSNE( const double* pSym, const double* q, const int& nData )
{
    double cost = 0.;

    // Since both distributions are symmetric and with null diagonal, it is 
    // enough to iterate over the upper triangle and multiply by two
    for ( int i = 0; i < (int)( nData * (nData - 1.) / 2. ); ++i )
        cost += 2. * pSym[ i ] * std::log( pSym[ i ] / q[ i ] );

    return cost;
}

// Compute the gradient of the cost function for t-SNE
void computeGradientCostFortSNE( const int& i, const double* pSym, const double* q,
                                 const double* U, const double* distSqU,
                                 const int& nData, const int& dimReduced,
                                 double* gradient )
{
    // Iterate over the dimensions of the vector
    for ( int k = 0; k < dimReduced; ++k )
    {
        // Set to zero the component of the gradient
        gradient[ k ] = 0.;

        // Sum over the other points
        for ( int j = 0; j < nData; ++j )
        {
            // If the point is the same, continue
            if ( j == i )
                continue;
            else
            {
                // Index for the upper-triangle matrices
                int index = mapIndicesUpperTriangSym( i, j, nData );

                if ( index > nData * (nData-1.)/2.)
                    std::cout << "Wrong index\n";

                // Sum the corresponding value to the component of the gradient
                gradient[ k ] += 4. * ( pSym[index] - q[index] ) *
                                 1. / ( 1. + distSqU[index] ) *
                                 ( U[ i*dimReduced + k ] - U[ j*dimReduced + k ] );
            }
        }
    }
}

// t-SNE (t-stochastic neighbor embedding)
//      data: the dataset
//      nData: the number of points
//      dim: the number of features
//      reducedData: dataset with a reduced number of features
//      dimReduced: the number of features after the reduction
//      perplexity: the perplexity parameter
//      gdRate: the rate parameter for gradient descent
//      gdMomentum: the momentum parameter for gradient descent
void computetSNE( const double* data, const int nData, const int dim, 
                  double* reducedData, const int dimReduced,
                  const double perplexity, const double gdRate,
                  const double gdMomentum )
{
    std::cout << "Applying t-SNE on the data...\n";

    // Setup the random seed
    std::srand( time(0) );

    // Tolerance for the gradient descent and maximum number of iterations
    double gdTolerance = 0.000001;
    int maxIterations = 5000;
    int itersForPrint = maxIterations / 20;
    int itersForEndExaggeration = itersForPrint * 2;

    std::cout << "Computing the matrix p_{j|i}... ";
    // Compute the components of p_{j|i}
    double* p = new double [ nData * nData ];
    computeP_ijFortSNE( p, data, nData, dim, perplexity );
    std::cout << "Done.\n\n";

    // From p_{j | i} above, compute the symmetrized version p_{ij}
    // Since p_{ii} = 0, we only need to compute the upper triangle
    double* pSym = new double [ (int)( nData * (nData - 1.) / 2. ) ];
    for ( int i = 0; i < nData; ++i )
    {
        for ( int j = i + 1; j < nData; ++j )
            pSym[ mapIndicesUpperTriang( i, j, nData ) ] = ( p[i*nData + j] + p[j*nData + i] ) / ( 2. * nData );
    }
    // Early exaggeration on the values of p
    for ( int i = 0; i < (int)( nData * (nData - 1.) / 2. ); ++i)
        pSym[i] *= 4.; 

    double sump = 0.;
    for ( int i = 0; i < nData; ++i )
    {
        for ( int j = 0; j < nData; ++j )
        {
            if ( i != j )
                sump += pSym[ mapIndicesUpperTriangSym( i, j, nData ) ];
        }
    }
    // Delete the array p now, since it won't be needed anymore.
    delete[] p;

    // Initialize the projection variables U_i to random values in the range [-1, 1]
    double* U = new double [ nData * dimReduced ];
    for (int i = 0; i < nData*dimReduced; ++i)
        U[i] = 2. * std::rand() / (double)RAND_MAX - 1.;

    // Initialize the variable for the distances between the U_i
    double* distSqU = new double [ (int)( nData * (nData - 1.) / 2. ) ];

    // Initialize the similar probability distribution q_{ij} with the values of 
    // U_i above.
    // This is also a symmetric matrix with the diagonal elements equal to zero.
    double* q = new double [ (int)( nData * (nData - 1.) / 2. ) ];

    // Compute the cost function
    double cost = 0.;
    double previousCost;

    // Vector for the velocities of the variables U_i, used in gradient descent
    double* velU = new double [ nData * dimReduced ];
    // Initialize these to zero
    for ( int i = 0; i < nData * dimReduced; ++i )
        velU[ i ] = 0.;

    // Vector to store the gradient computed for each point in each iteration
    double* gradient = new double [ dimReduced ];

    std::cout << "Minimizing the cost function for t-SNE...\n";
    // Perform gradient descent until a threshold condition is met
    int iteration = 0;
    do
    {
        // Compute the distances between the variables U_i
        computeDistSqU( U, nData, dimReduced, distSqU );
        // Compute the elements of the similar probability distribution q_{ij}
        computeSimilarDistrtSNE( q, distSqU, nData );
        // Store the previous value of the cost function
        previousCost = cost;
        // Compute the cost function
        cost = computeCosttSNE( pSym, q, nData );

        // Iterate through the different points U_i
        for ( int i = 0; i < nData; ++i )
        {
            // Compute the corresponding gradient of the cost function
            computeGradientCostFortSNE( i, pSym, q, U, distSqU, nData, 
                                        dimReduced, gradient );

            // Iterate through the components of U_i
            for ( int k = 0; k < dimReduced; ++k )
            {
                // Update the velocity, with the contribution of the previous value
                // (momentum) and the gradient
                velU[ i*dimReduced + k ] = gdMomentum * velU[ i*dimReduced + k ] + gdRate * gradient[ k ];
            }
        }

        // Move U along the new velocity
        for ( int i = 0; i < nData * dimReduced; ++i )
            U[ i ] -= velU[ i ];

        if ( iteration % itersForPrint == 0 )
        {
            std::cout << "Iteration " << iteration <<  " of " << maxIterations 
                      <<  ", cost function: " << cost 
                      << ". Previous value of the cost function: " << previousCost << '\n';
            if ( iteration == itersForEndExaggeration )
            {
                // Revert early exageration on the values of the matrix p
                for ( int i = 0; i < (int)( nData * (nData - 1.) / 2. ); ++i)
                    pSym[i] /= 4.; 
            }
        }
    }
    while ( std::fabs( cost - previousCost ) > gdTolerance && iteration++ < maxIterations );

    // Save the reduced data in the output variable
    for ( int i = 0; i < nData * dimReduced; ++i )
        reducedData[i] = U[i];

    // Delete the previously defined arrays
    // delete[] p;
    delete[] pSym;
    delete[] U;
    delete[] distSqU;
    delete[] q;
    delete[] velU;
    delete[] gradient;
}
