#ifndef BASICSSTATISTICS_H
#define BASICSSTATISTICS_H

#include <iostream>
#include <vector>
#include <cmath>

// Means of all the features in a dataset
void compMeansDataset( const double* data, const int& nData, const int& dim,
                       double* means );

// Standard deviations of all the features in a dataset
void compStdDevDataset( const double* data, const int& nData, const int& dim,
                        double* means, double* stdDevs );

// Correlation matrix of a dataset
void compCorrMatrixDataset( const double* data, const int& nData, const int& dim,
                            double* means, double* corrMatrix );

// Normalize a dataset
void normalizeDataset( const double* data, const int& nData, const int& dim,
                       double* normDataset, double* means );

#endif
