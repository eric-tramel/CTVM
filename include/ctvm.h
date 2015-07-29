#ifndef CTVM_H
#define CTVM_H


#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "ctvm_util.h"

BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long pixel);
BoostDoubleMatrix Gradient2DMatrix(BoostDoubleVector U);
double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu);
BoostDoubleVector Shrike(BoostDoubleVector DiUk, BoostDoubleVector NUi, double beta);

double u_subfunction(BoostDoubleMatrix A, double u, BoostDoubleVector B, BoostDoubleVector Dl, BoostDoubleVector Wl, BoostDoubleVector NUl, BoostDoubleVector LAMBDA, double beta, double mu);
// double onestep_direction(BoostDoubleMatrix A, double U, BoostDoubleVector B, BoostDoubleVector Wl, BoostDoubleVector NUl, BoostDoubleVector LAMBDA, double beta, double mu);
// BoostDoubleVector alternating_minimisation(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu, unsigned long rank);

// BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles);

#endif