#include "basicsStatistics.h"

// Means of all the features in a dataset
void compMeansDataset( const double* data, const int& nData, const int& dim,
                       double* means )
{
    // Iterate through the features
    for (int i = 0; i < dim; ++i)
    {
        means[i] = 0.;
        // Iterate through the data points
        for (int j = 0; j < nData; ++j)
            means[i] += data[j*dim + i];
        // Divide the mean by the number of points
        means[i] /= nData;
    }
}

// Standard deviations of all the features in a dataset
void compStdDevDataset( const double* data, const int& nData, const int& dim,
                        double* means, double* stdDevs )
{
    // Iterate through the features
    for (int i = 0; i < dim; ++i)
    {
        stdDevs[i] = 0.;
        // Iterate through the data points
        for (int j = 0; j < nData; ++j)
            stdDevs[i] += ( data[j*dim + i] - means[i] ) * (data[j*dim + i] - means[i] );
        // Divide the standard deviation by the number of points, and compute its
        // square root
        stdDevs[i] = std::sqrt( stdDevs[i] / nData );
    }
}

// Correlation matrix of a dataset
void compCorrMatrixDataset( const double* data, const int& nData, const int& dim,
                            double* means, double* corrMatrix )
{
    // Iterate through the features twice
    for (int i1 = 0; i1 < dim; ++i1)
    {
        for (int i2 = 0; i2 < dim; ++i2)
        {
            // Initialize to zero this component of the correlation matrix
            corrMatrix[i1*dim + i2] = 0.;
            // Iterate through the data points
            for (int j = 0; j < nData; ++j)
            {
                corrMatrix[i1*dim + i2] += ( data[j*dim+i1] - means[i1] ) * ( data[j*dim+i2] - means[i2] );
            }
            // Divide by the number of points
            corrMatrix[i1*dim + i2] /= nData;
        }
    }
}

// Normalize a dataset
void normalizeDataset( const double* data, const int& nData, const int& dim,
                       double* normDataset, double* means )
{
    // Compute the means
    compMeansDataset( data, nData, dim, means );

    // Compute the standard deviations
    double* stdDevs = new double [ dim ];
    compStdDevDataset( data, nData, dim, means, stdDevs );

    // Normalize the data
    for (int i = 0; i < nData; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            normDataset[i*dim + j] = ( data[i*dim + j] - means[j] ) / stdDevs[j];
        }
    }
}
