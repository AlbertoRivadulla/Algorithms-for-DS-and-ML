#include "dimensionalReduction.h"

// Principal component analysis of a dataset
//      data: the dataset
//      nData: the number of points
//      dim: the number of features
//      reducedData: dataset with a reduced number of features
//      dimReduced: the number of features after the reduction
double* computePCA( const double* data, const int nData, const int dim, 
                    double* reducedData, const int dimReduced )
// void computePCA( const double* data, const int nData, const int dim, 
//                  double* reducedData, const int dimReduced )
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

    return mostSignificantEigenvectors;
}
