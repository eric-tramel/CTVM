#ifndef CTVM_H
#define CTVM_H

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>
#include "ctvm_util.h"



/* Gradient Operations */
#define HORZ 0
#define VERT 1
BoostDoubleVector PixelGradient(BoostDoubleVector X, unsigned long Index,
                                unsigned int SideLength);
BoostDoubleMatrix AllPixelGradients(BoostDoubleVector X, unsigned long SideLength);
BoostDoubleVector PixelGradientAdjointSum(BoostDoubleMatrix G, unsigned long SideLength);



// TODO: Implement the following
// Periodic 2D Gradients...
BoostDoubleVector PeriodicPixelGradient(BoostDoubleVector X, unsigned long Index,
                                unsigned int SideLength);
BoostDoubleMatrix AllPeriodicPixelGradients(BoostDoubleVector X, unsigned long SideLength);
BoostDoubleVector PeriodicPixelGradientAdjointSum(BoostDoubleMatrix G, unsigned long SideLength);

// 3D Gradients...
BoostDoubleVector VoxelGradient(BoostDoubleVector g, unsigned long index,
                                unsigned int SideLength);
BoostDoubleMatrix AllVoxelGradients(BoostDoubleVector X, unsigned int SideLength);
BoostDoubleVector VoxelGradientAdjointSum(BoostDoubleVector X, unsigned long index,
                                          unsigned int SideLength);
// ENDTODO


/* Shrinkage-like Operators */
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
                  BoostDoubleMatrix Nu, BoostDoubleVector Lambda, 
                  double beta, double mu,
                  unsigned long SideLength, TVType GradNorm);

BoostDoubleVector Onestep_Direction(BoostDoubleMatrix A, BoostDoubleVector Uk, 
                                    BoostDoubleVector B, BoostDoubleMatrix Wk, 
                                    BoostDoubleMatrix Nu, BoostDoubleVector Lambda, 
                                    double beta, double mu,
								    unsigned long SideLength);

double U_Subfunction(BoostDoubleMatrix A, BoostDoubleVector U, 
                     BoostDoubleVector B, BoostDoubleMatrix Wk, 
                     BoostDoubleMatrix Nu, BoostDoubleVector Lambda, 
                     double beta, double mu, unsigned long SideLength);

void Alternating_Minimisation(BoostDoubleMatrix A, BoostDoubleVector &U, 
                              BoostDoubleVector B, BoostDoubleMatrix &W, 
                              BoostDoubleMatrix Nu, BoostDoubleVector Lambda, 
                              double beta, double mu, 
							                unsigned long SideLength);

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix A, BoostDoubleVector y,
									   unsigned long SideLength);

#endif