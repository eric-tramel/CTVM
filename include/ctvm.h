#ifndef CTVM_H
#define CTVM_H

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "ctvm_util.h"

BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long pixel);
BoostDoubleMatrix Gradient2DMatrix(BoostDoubleVector U);
BoostDoubleMatrix Unit_Gradient2DMatrix(BoostDoubleVector U, unsigned long pixel);

/* Shrinkage-like Operators */
// EWT: I have verified the operation of both of these Shrike operators
//      (the Anisotropic and Isotropic versions), so we can be assured in
//      their correct operation.
enum TVType {ISOTROPIC,ANISOTROPIC};
BoostDoubleVector ShrikeIsotropic(BoostDoubleVector W, BoostDoubleVector Nu, 
                                  double beta);
BoostDoubleVector ShrikeAnisotropic(BoostDoubleVector W, BoostDoubleVector Nu, 
                                    double beta);
BoostDoubleMatrix ApplyShrike(BoostDoubleMatrix AllW, BoostDoubleMatrix AllNu,
                              double beta, TVType ShrikeMode);

double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, 
                  BoostDoubleVector B, BoostDoubleMatrix W, 
                  BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                  double beta, double mu);

BoostDoubleVector Onestep_direction(BoostDoubleMatrix A, BoostDoubleVector Uk, 
                                    BoostDoubleVector B, BoostDoubleMatrix Wk, 
                                    BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                                    double beta, double mu);

double U_subfunction(BoostDoubleMatrix A, BoostDoubleVector U, 
                     BoostDoubleVector B, BoostDoubleMatrix Wk, 
                     BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                     double beta, double mu);

BoostDoubleMatrix alternating_minimisation(BoostDoubleMatrix A, BoostDoubleVector U, 
                                           BoostDoubleVector B, BoostDoubleMatrix W, 
                                           BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                                           double beta, double mu);

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles);

#endif