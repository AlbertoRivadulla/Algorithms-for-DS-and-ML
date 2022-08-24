#ifndef DIMENSIONAL_REDUCTION_H
#define DIMENSIONAL_REDUCTION_H

#include <iostream>
#include <random>

#include "linalg.h"
#include "basicsStatistics.h"
#include "DSAlgorithms.h"
#include "auxFunctions.h"


//------------------------------------------------------------------------------
// t-SNE

// Compute the components of the matrix p_ij for t-SNE, from an array of data
// points of a given dimension, and with a desired perplexity
void computeP_ijFortSNE( double* p, const double* data, const int& nData,
                         const int& dim, const double& perplexity);

// Compute the distances squared between the variables U_i.
void computeDistSqU( const double* U, const int& nData, const int& dimReduced,
                     double* distSqU );

// Compute the elements of the similar probability distribution q_{ij}
void computeSimilarDistrtSNE( double* q, const double* distSqU, const int& nData );

// Compute the cost function as the Kullback-Leibler divergence of the probability 
// distributions q_{ij} and p_{ij}
double computeCosttSNE( const double* p, const double* q, const int& nData );

// Compute the gradient of the cost function for t-SNE
void computeGradientCostFortSNE( const int& i, const double* pSym, const double* q,
                                 const double* U, const double* distSqU,
                                 const int& nData, const int& dimReduced,
                                 double* gradient );

#endif
