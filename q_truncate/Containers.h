// ==============================================================================
//
//  Containers.h
//  QTR
//
//  Created by Albert Lu on 8/19/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/30/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_CONTAINERS_H
#define QTR_CONTAINERS_H

#include <complex>
#include <vector>
//#include "boost/multi_array.hpp"

/*typedef boost::multi_array<std::complex<double>, 2> MeshCX2D;
typedef boost::multi_array<double, 2> MeshD2D;
typedef boost::multi_array<std::complex<double>, 3> MeshCX3D;
typedef boost::multi_array<double, 3> MeshD3D;
typedef boost::multi_array<std::complex<double>, 4> MeshCX4D;
typedef boost::multi_array<double, 4> MeshD4D;
typedef boost::multi_array<bool, 2> Mask2D;
typedef boost::multi_array<bool, 3> Mask3D;
typedef boost::multi_array<bool, 4> Mask4D;*/
typedef std::vector<int> MeshIndex;
typedef std::vector<std::vector<int>> MeshIndexFull;
typedef std::vector<unsigned long int> MeshIndexLU;
typedef std::vector<unsigned int> MeshIndexU;
typedef std::vector<double> MeshValueD;
//typedef std::vector<std::vector<int>> MeshGrid;

#endif /* QTR_CONTAINERS_H */
