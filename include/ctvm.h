#ifndef CTVM_H
#define CTVM_H

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "ctvm_util.h"



/* Gradient Operations */
#define HORZ 0
#define VERT 1
BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long pixel);
BoostDoubleMatrix Gradient2DMatrix(BoostDoubleVector U);
BoostDoubleMatrix Unit_Gradient2DMatrix(BoostDoubleVector U, unsigned long pixel);

BoostDoubleVector PixelGradient(BoostDoubleVector X, unsigned long Index,
                                unsigned int SideLength);
BoostDoubleMatrix AllPixelGradients(BoostDoubleVector X, unsigned int SideLength);
BoostDoubleVector PixelGradientAdjointSum(BoostDoubleMatrix G, unsigned int Sidelength);

// TODO: Implement the following
// 2D Gradients...
// BoostDoubleVector PixelGradientAdjoint(BoostDoubleVector g, unsigned long index,
//                                        unsigned int SideLength);
// 3D Gradients...
BoostDoubleVector VoxelGradient(BoostDoubleVector g, unsigned long index,
                                unsigned int SideLength);
BoostDoubleVector VoxelGradientAdjoint(BoostDoubleVector X, unsigned long index,
                                       unsigned int SideLength);
BoostDoubleMatrix AllVoxelGradients(BoostDoubleVector X, unsigned int SideLength);
// ENDTODO

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

/* Optimization */
double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, 
                  BoostDoubleVector B, BoostDoubleMatrix W, 
                  BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                  double beta, double mu);

BoostDoubleVector Onestep_Direction(BoostDoubleMatrix A, BoostDoubleVector Uk, 
                                    BoostDoubleVector B, BoostDoubleMatrix Wk, 
                                    BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                                    double beta, double mu, unsigned long l);

double U_subfunction(BoostDoubleMatrix A, BoostDoubleVector U, 
                     BoostDoubleVector B, BoostDoubleMatrix Wk, 
                     BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                     double beta, double mu);

BoostDoubleMatrix alternating_minimisation(BoostDoubleMatrix A, BoostDoubleVector U, 
                                           BoostDoubleVector B, BoostDoubleMatrix W, 
                                           BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, 
                                           double beta, double mu, unsigned long l);

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles);

#endif