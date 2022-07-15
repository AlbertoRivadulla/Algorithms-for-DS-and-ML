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
//      means: array to output the positions of the means of the clusters
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



//------------------------------------------------------------------------------
// Agglomerative clustering

// Map the indices (i, j) of a square matrix to indices of an array that represents
// sequentially the elements of the upper triangle of the original one
inline int mapIndices( const int& i, const int& j, const int& n )
{
    return j - 1 - (i*i + 3*i) / 2 + i*n;
}

// Compute the distances between pairs of clusters using Ward linkage
void computeWardClusterDistances( const double* means, const int* nPointsClusters, 
                                  const int dim, double* distances, 
                                  const int& nPoints )
{
    // Iterate through the clusters
    for (int i = 0; i < nPoints; ++i)
    {
        for (int j = i+1; j < nPoints; ++j)
        {
            // Check if there are points in the two clusters
            if ( nPointsClusters[i] != 0 && nPointsClusters[j] != 0 )
            {
                // Compute the distances between the means
                double distSqMeans = 0.;
                for (int k = 0; k < dim; ++k)
                {
                    distSqMeans += ( means[i*dim+k] - means[j*dim+k] ) * 
                                   ( means[i*dim+k] - means[j*dim+k] );
                }
                // Store the value in the matrix of distances
                // The expression in the index maps (i,j) to sequential indices of
                // the upper half diagonal
                distances[ mapIndices( i, j, nPoints ) ] = distSqMeans * nPointsClusters[i] * 
                                                           nPointsClusters[j] / (nPointsClusters[i] + nPointsClusters[j]);
            }
            else
            {
                distances[ mapIndices( i, j, nPoints ) ] = 0.;
            }
        }
    }
}

// Agglomerative clustering, with Ward linkage for the distances
//      data: the dataset
//      nData: the number of points
//      dim: the number of features
//      nClusters: the number of clusters to form
//      clusterIndices: array to output the cluster index of each data point
//      means: array to output the positions of the means of the clusters
// This function is declared in the header "DSAlgorithms.h"
void computeAgglomerativeClusters( const double* data, const int nData, 
                                   const int dim, const int nClusters, 
                                   int* clusterIndices, double* means )
{
    std::cout << "Computing agglomerative clusters...\n";

    // Initialize the cluster means
    double* clusterMeans = new double [ dim * nData ];
    for (int i = 0; i < dim*nData; ++i)
        clusterMeans[i] = data[i];
    // Initialize the cluster counters
    int* nPointsClusters = new int [ nData ];
    for (int i = 0; i < nData; ++i)
        nPointsClusters[i] = 1;

    // Matrix of distances between clusters
    double* distances = new double[ (int)(nData * (nData - 1.) / 2.) ];
    // Number of clusters
    int nClustersCurrent = nData;
    // In the beginning, each point is its own cluster
    // Set this in the cluster labels
    for (int i = 0; i < nData; ++i)
        clusterIndices[i] = i;

    // Iterate until the desired amount of clusters is reached
    int iterations = 0;
    do
    {
        // Compute the cluster distances
        computeWardClusterDistances( clusterMeans, nPointsClusters, dim, distances, nData );

        // Find the clusters with the smallest distance
        int iMin1 = 0;
        int iMin2 = 1;
        double minDistSq = std::numeric_limits<double>::max();
        for (int i = 0; i < nData; ++i)
        {
            for (int j = i+1; j < nData; ++j)
            {
                if ( distances[ mapIndices( i, j, nData ) ] < minDistSq && nPointsClusters[i] != 0 && nPointsClusters[j] != 0 )
                {
                    minDistSq = distances[ mapIndices( i, j, nData ) ];
                    iMin1 = i;
                    iMin2 = j;
                }
            }
        }

        // Merge the clusters with the smallest distance
        // Update cluster labels
        for ( int i = 0; i < nData; ++i )
        {
            if ( clusterIndices[i] == iMin2 )
                clusterIndices[i] = iMin1;
        }

        // Update cluster means
        for ( int j = 0; j < dim; ++j )
        {
            clusterMeans[ iMin1*dim + j ] *= nPointsClusters[ iMin1 ];
            clusterMeans[ iMin2*dim + j ] *= nPointsClusters[ iMin2 ];
            clusterMeans[ iMin1*dim + j ] += clusterMeans[ iMin2*dim + j ];
            clusterMeans[ iMin1*dim + j ] /= ( nPointsClusters[ iMin2 ] + nPointsClusters[ iMin1 ] );
        }

        // Update cluster counts
        nPointsClusters[ iMin1 ] += nPointsClusters[ iMin2 ];
        nPointsClusters[ iMin2 ] = 0;

        // Subtract one to the number of clusters
        nClustersCurrent -= 1;

        ++iterations;
    }
    while ( nClustersCurrent > nClusters );

    // Put the cluster means in the output array
    for ( int i = 0; i < nClusters; ++i )
    {
        for ( int j = 0; j < dim; ++j )
            means[ i*dim+j ] = clusterMeans[ i*dim+j ];
    }

    // Map the cluster indices to numbers in ( 0, nClusters - 1 )
    // Get the numbers in this range that do not label any cluster
    std::vector<int> freeIndices;
    for (int i = 0; i < nClusters; ++i)
    {
        if ( nPointsClusters[i] == 0 )
            freeIndices.push_back( i );
    }
    // Iterate through the cluster labels
    if ( freeIndices.size() > 0 )
    {
        for ( int i = 0; i < nData; ++i )
        {
            if ( clusterIndices[i] >= nClusters )
            {
                int thisIndex = clusterIndices[i];
                // Replace the number of points in the clusters
                nPointsClusters[freeIndices[freeIndices.size()-1]] = nPointsClusters[clusterIndices[i]];
                nPointsClusters[clusterIndices[i]] = 0;

                // Replace the means
                for ( int j = 0; j < dim; ++j )
                    means[ freeIndices[freeIndices.size()-1]*dim+j ] = clusterMeans[ thisIndex*dim+j ];

                // Replace this index by the last in the vector of free indices
                for ( int j = i; j < nData; ++j )
                {
                    if ( clusterIndices[j] == thisIndex )
                        clusterIndices[j] = freeIndices[freeIndices.size()-1];
                }

                // Delete the index
                freeIndices.pop_back();
            }
        }
    }

    std::cout << "Data classified in " << nClusters << " clusters, found after " << iterations << " iterations.\n";
}


//------------------------------------------------------------------------------
// DBSCAN clustering

// Get the neighbors of a point index in the dataset, within a radius epsilon
std::queue<int> getNeighbors( const double* data, const int& dim, const int& nData, 
                              const int& index, const double& epsilonSq )
{
    // Initialize the empty queue
    std::queue<int> neighbors;

    // Iterate over all the points in the dataset
    for ( int j = 0; j < nData; ++j )
    {
        // Check that they are not the same point!
        if ( j == index )
            continue;

        // Compute the distance
        double distSq = 0;
        for ( int k = 0; k < dim; ++k )
            distSq += ( data[j*dim+k] - data[index*dim+k] ) * 
                      ( data[j*dim+k] - data[index*dim+k] );
        
        // If the distance is smaller than the threshold, add it to the list of 
        // neighbors
        if ( distSq < epsilonSq )
            neighbors.push( j );
    }

    return neighbors;
}

// DBSCAN clustering
//      data: the dataset
//      nData: the number of points
//      dim: the number of features
//      epsilon: positive number that sets the scale of the clusters
//      nrNeighborsCore: minimum nr of neighbors for a point to be a core point
//      clusterLabels: array to output the cluster index of each data point
// This function is declared in the header "DSAlgorithms.h"
void computeDBSCANClusters( const double* data, const int nData, 
                            const int dim, const double epsilon, 
                            const int nrNeighborsCore, int* clusterLabels )
{
    std::cout << "Computing clusters with the DBSCAN algorithm...\n";

    // Compute the value of epsilon squared, to speed up later computations
    double epsilonSq = epsilon * epsilon;

    // Cluster counter
    int cluster = 1;

    // Initialize all the labels to minus one
    // This means the point has not been visited yet
    for (int i = 0; i < nData; ++i)
        clusterLabels[i] = -1;

    // Process all points in the dataset in a sequential manner
    for (int i = 0; i < nData; ++i)
    {
        // If the point is marked as unvisited (-1), continue
        if ( clusterLabels[ i ] > 0 )
            continue;

        // Find the neighbors of the point within the distance epsilon
        std::queue<int> neighbors = getNeighbors( data, dim, nData, i, epsilonSq );

        // Check if the number of neighbors is bigger than the minimum in order
        // to be a core-point
        if ( neighbors.size() > nrNeighborsCore )
        {
            // If it is not, (mark it as noise) and continue
            clusterLabels[ i ] = 0;
            continue;
        }

        // Label the point as belonging to the current cluster
        clusterLabels[ i ] = cluster;

        // Iterate over the neighbors of the point
        while ( !neighbors.empty() )
        {
            // Get the current neighbor and remove it from the queue
            int currentNeigh = neighbors.front();
            neighbors.pop();

            // If the point already belongs to a cluster, continue
            if ( clusterLabels[ currentNeigh ] > 0 )
                continue;

            // Otherwise, set the label of this point to the current cluster
            clusterLabels[ currentNeigh ] = cluster;

            // Compute the neighbors of this neighbor 
            std::queue<int> neighsOfNeigh = getNeighbors( data, dim, nData, 
                                                          currentNeigh, epsilonSq );

            // Check if the neighbor is also a core-point
            if ( neighsOfNeigh.size() > nrNeighborsCore )
            {
                // If it is, add them to the other list of neighbors
                while ( !neighsOfNeigh.empty() )
                {
                    int currentNeighOfNeigh = neighsOfNeigh.front();
                    neighsOfNeigh.pop();
                    if ( clusterLabels[ currentNeighOfNeigh ] <= 0 )
                        neighbors.push( currentNeighOfNeigh );
                }
                // neighbors.insert( neighsOfNeigh.begin(), neighsOfNeigh.end() );
            }
        }

        // Add one to the cluster count
        ++cluster;
    }

    std::cout << "\nData classified in " << cluster << " clusters.\n\n";
}
