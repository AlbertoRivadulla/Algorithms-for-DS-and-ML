#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <iostream>
#include <random>
#include <ctime>
#include <limits>

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

#endif
