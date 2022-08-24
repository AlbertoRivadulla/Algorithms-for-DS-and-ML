#ifndef DSALGORITHMS_H
#define DSALGORITHMS_H

#include <iostream>

#include "basicsStatistics.h"
#include "linalg.h"
#include "clustering.h"
#include "dimensionalReduction.h"

// Write here the wrappers for the functions that will be called from Python
extern "C"
{
    //--------------------------------------------------------------
    // Methods implemented in dimensionalReduction.cpp

    // Principal component analysis of a dataset
    //      data: the dataset
    //      nData: the number of points
    //      dim: the number of features
    //      reducedData: dataset with a reduced number of features
    //      dimReduced: the number of features after the reduction
    double* computePCA( const double* data, const int nData, const int dim, 
                        double* reducedData, const int dimReduced=2 );

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
                  const double gdMomentum );


    //--------------------------------------------------------------
    // Methods implemented in clustering.cpp

    // K-means clustering
    //      data: the dataset
    //      nData: the number of points
    //      dim: the number of features
    //      K: the number of clusters to form
    //      clusterIndices: array to output the cluster index of each data point
    //      means: array to output the positions of the means of the clusters
    void computeKMeansClusters( const double* data, const int nData, const int dim,
                                const int K, int* clusterIndices, double* means );

    // Agglomerative clustering
    //      data: the dataset
    //      nData: the number of points
    //      dim: the number of features
    //      nClusters: the number of clusters to form
    //      clusterIndices: array to output the cluster index of each data point
    //      means: array to output the positions of the means of the clusters
    void computeAgglomerativeClusters( const double* data, const int nData, 
                                       const int dim, const int nClusters, 
                                       int* clusterIndices, double* means );

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
                                const int nrNeighborsCore, int* clusterLabels );
}

#endif
