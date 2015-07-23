#ifndef CTVM_H
#define CTVM_H

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "ctvm_util.h"

typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
typedef boost::numeric::ublas::scalar_matrix<double> BoostScalarDoubleMatrix;
typedef boost::numeric::ublas::vector<double> BoostDoubleVector;
typedef boost::numeric::ublas::scalar_vector<double> BoostScalarDoubleVector;

BoostDoubleVector 2DGradient(BoostDoubleVector U, unsigned long pixel);

// BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles);

#endif