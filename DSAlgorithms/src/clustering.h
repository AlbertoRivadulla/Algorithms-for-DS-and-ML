#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <iostream>
#include <random>
#include <ctime>
#include <limits>
// #include <iterator>
// #include <set>
#include <queue>

#include "linalg.h"
#include "basicsStatistics.h"
#include "DSAlgorithms.h"
#include "auxFunctions.h"

//------------------------------------------------------------------------------
// K-means clustering

// Initialize the means randomly, between the minimum and maximum values in each
// direction
void initializeMeans( const double* data, const int nData, const int dim, 
                      const int K, double* means, double& scale );

// Distance between the point i and the mean k
float distSqMeanToPoint( const double* data, const int& dim, const double* means,
                         const int& i, const int& k);

// Maximization step
// Get the cluster assignments so the distances are minimized
// Update also the cost function
void maximizationStep( const double* data, const int nData, const int dim,
                       const int K, const double* means, int* assignments, 
                       double& cost );

// Expectation step
// Compute the cluster means with the updated cluster assignments
void expectationStep( const double* data, const int nData, const int dim,
                      const int K, double* means, const int* assignments );


//------------------------------------------------------------------------------
// Agglomerative clustering

// Map the indices (i, j) of a square matrix to indices of an array that represents
// sequentially the elements of the upper triangle of the original one
inline int mapIndices( const int& i, const int& j, const int& n );

// Compute the distances between pairs of clusters using Ward linkage
void computeWardClusterDistances( const double* means, const int* nPointsClusters, 
                                  const int dim, double* distances, 
                                  const int& nPoints );


//------------------------------------------------------------------------------
// DBSCAN clustering

// Get the neighbors of a point index in the dataset, within a radius epsilon
std::queue<int> getNeighbors( const double* data, const int& dim, const int& nData, 
                              const int& index, const double& epsilonSq );

#endif
