#include "ctvm.h"

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, 
                                       BoostDoubleVector TiltAngles){
/*
    Function: tval3_reconstruction
    ------------------------------
    Perform two-dimensional TVAL3 tomographic reconstruction given a
    set of measurements (Sinogram) and the measured tilt angles.
*/

	/* <TODO: Everything!> */
	unsigned double l = Sinogram.size1(); /* Size of the sample (in pixels) */
	unsigned double o = Sinogram.size2(); /* Numbers of tilt angles */
	unsigned double m = l * o; /* Numbers of measurements */
	unsigned double n = l * l;
	unsigned double nu, lambda, beta, mu;
	unsigned double
	unsigned double
	unsigned double Lagrangian;
	float alpha = 1.05;
	float tol = 0.5;
	float delta, rho, eta;
	
	BoostDoubleMatrix W (rows, cols);
	BoostDoubleMatrix A (m, n);
	
	BoostDoubleVector Y (m);
	BoostDoubleVector X (n);
	BoostDoubleVector NU (/*?*/);
	
	A = CreateRandomMatrix(m, n);
	Y = MatrixToVector(Sinogram);
	
	do{
		
	
	} while (Lagrangian > tol)


    BoostDoubleMatrix RecoveredImage(32,32);  // Create a dummy matrix to return

return RecoveredImage;
}