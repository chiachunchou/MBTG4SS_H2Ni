// ==============================================================================
//
//  Eigen.h
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 8/4/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_EIGEN_H
#define QTR_EIGEN_H

#include "./libraries/Eigen/Dense"
#include "./libraries/Eigen/Eigenvalues"

using namespace Eigen;

typedef Eigen::Matrix<double,Eigen::Dynamic,3> AtomMatrix;
typedef Eigen::Matrix<int,Eigen::Dynamic,3> AtomMatrixI;

#endif /* QTR_EIGEN_H */
