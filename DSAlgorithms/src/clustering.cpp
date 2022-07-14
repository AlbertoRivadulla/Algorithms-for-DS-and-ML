#include "clustering.h"

//------------------------------------------------------------------------------
// K-means clustering

// Initialize the means randomly, between the minimum and maximum values in each
// direction
void initializeMeans( const double* data, const int nData, const int dim, 
                      const int K, double* means, double& scale )
{
    // Iterate over the dimensions
    for (int i = 0; i < dim; ++i)
    {
        // Get the maximum and maximum values of this feature in the dataset
        float max = std::numeric_limits<float>::min();
        float min = std::numeric_limits<float>::max();
        for (int j = 0; j < nData; ++j)
        {
            if ( data[j*dim + i] > max )
                max = data[j*dim + i];
            if ( data[j*dim + i] < min )
                min = data[j*dim + i];
        }

        // Initialize this component of each mean
        for (int k = 0; k < K; ++k)
            means[ k*dim + i ] = (max - min) * std::rand() / (float)RAND_MAX + min;

        // If max-min is larger than the scale, update it
        double newScale = max - min;
        if ( newScale > scale )
            scale = newScale;
    }
}

// Distance between the point i and the mean k
float distSqMeanToPoint( const double* data, const int& dim, const double* means,
                         const int& i, const int& k)
{
    float distSq = 0.;
    for (int j = 0; j < dim; ++j)
        distSq += ( data[dim*i + j] - means[dim*k + j] ) * ( data[dim*i + j] - means[dim*k + j] );
    return distSq;
}

// Maximization step
// Get the cluster assignments so the distances are minimized
// Update also the cost function
void maximizationStep( const double* data, const int nData, const int dim,
                       const int K, const double* means, int* assignments, 
                       double& cost )
{
    // Initialize the cost function
    cost = 0.;

    // Iterate over the data points
    for (int i = 0; i < nData; ++i)
    {
        // Minimum distance
        double minDistSq = std::numeric_limits<float>::max();

        // Iterate through the cluster means
        for (int k = 0; k < K; ++k)
        {
            // Compute the distance between this cluster mean and the point i
            float distanceSq = distSqMeanToPoint( data, dim, means, i, k );
            // If the distance is smaller than the one found previously, store
            // the cluster index
            if ( distanceSq < minDistSq )
            {
                assignments[ i ] = k;
                minDistSq = distanceSq;
            }
        }
        // Update the cost function, adding to it the smallest distance squared
        cost += minDistSq;
    }
}

// Expectation step
// Compute the cluster means with the updated cluster assignments
void expectationStep( const double* data, const int nData, const int dim,
                      const int K, double* means, const int* assignments )
{
    // Counters of the number of points in each cluster
    int* counters = new int [ K ];
    for (int k = 0; k < K; ++k)
        counters[ k ] = 0;

    // Set all the components of the means to zero
    for (int i = 0; i < K*dim; ++i)
        means[ i ] = 0.;

    // Iterate through the points
    for (int i = 0; i < nData; ++i)
    {
        // Add the location of the current point to the corresponding mean
        for (int j = 0; j < dim; ++j)
            means[ assignments[i]*dim + j ] += data[ i*dim + j ];
        // Add one to the corresponding counter
        counters[ assignments[i] ] += 1;
    }

    // Divide all the means by the number of points in its cluster
    for (int k = 0; k < K; ++k)
    {
        for (int j = 0; j < dim; ++j)
            means[ k*dim + j ] /= (float)counters[ k ];
    }
}

// K-means clustering
//      data: the dataset
//      nData: the number of points
//      dim: the number of features
//      K: the number of clusters to form
//      clusterIndices: array to output the cluster index of each data point
// This function is declared in the header DSAlgorithms.h
void computeKMeansClusters( const double* data, const int nData, const int dim,
                            const int K, int* clusterIndices, double* means )
{
    // Setup the random seed
    std::srand( time(0) );

    // // Vector of cluster means
    // double* means = new double [ K * dim ];
    // // Vector of cluster assignments
    // int* assignments = new int [ nData ];
    // Initialize the cost (objective) function
    double cost = 0.;
    double previousCost = 0.;
    // An scale of the problem, to use when computing the end procedure for the loop
    double scale;

    // Initialize the means randomly in the plane
    initializeMeans( data, nData, dim, K, means, scale );

    // Iterate the algorithm until the change in the sum of distances is small
    int iteration = 0;
    do
    {
        // Copy the previous value of the cost function
        previousCost = cost;

        // Maximization step
        // Get the cluster assignments so the distances are minimized
        // Update also the cost function
        maximizationStep( data, nData, dim, K, means, clusterIndices, cost );

        // Expectation step
        // Compute the cluster means with the updated cluster assignments
        expectationStep( data, nData, dim, K, means, clusterIndices );

        iteration++;
    }
    while ( std::fabs( cost - previousCost ) > 0.00001 * scale && iteration < 300 );

    std::cout << "Data classified in " << K << " clusters, found after " << iteration << " iterations.\n";
}
