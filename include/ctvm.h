#ifndef CTVM_H
#define CTVM_H


#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "ctvm_util.h"

BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long pixel);
BoostDoubleMatrix Gradient2DMatrix(BoostDoubleVector U);
double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu);

// BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles);

#endif