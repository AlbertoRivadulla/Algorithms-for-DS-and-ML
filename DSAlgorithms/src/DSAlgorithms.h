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
    // void computePCA( const double* data, const int nData, const int dim, 
    //                  double* reducedData, const int dimReduced=2 );
    double* computePCA( const double* data, const int nData, const int dim, 
                        double* reducedData, const int dimReduced=2 );


    //--------------------------------------------------------------
    // Methods implemented in clustering.cpp

    // K-means clustering
    //      data: the dataset
    //      nData: the number of points
    //      dim: the number of features
    //      K: the number of clusters to form
    //      clusterIndices: array to output the cluster index of each data point
    void computeKMeansClusters( const double* data, const int nData, const int dim,
                                const int K, int* clusterIndices, double* means );
}

#endif
