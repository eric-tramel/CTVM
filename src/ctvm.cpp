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
	unsigned double nu = 1;
	unsigned double lambda = 1;
	unsigned double beta = 1;
	unsigned double mu = 1;
	unsigned double
	unsigned double Lagrangian; // ?
	float alpha = 1.05;
	float tol = 0.5;
	float delta, rho, eta;
	
	BoostDoubleMatrix W (rows, cols); // w(i)?
	BoostDoubleMatrix A (m, n);
	
	BoostDoubleVector Y (m);
	BoostDoubleVector X (n);
	BoostDoubleVector NU (); // size ?
	
	A = CreateRandomMatrix(m, n);
	Y = MatrixToVector(Sinogram); // u(i) ?
	
	do{
		for (unsigned int k = 0; k < m; ++k){
			w = w(k+1);
			u = u(k+1);
			for (unsigned int i = 0; i < m; ++i){
				L(w,u) = L(w(k),u(k)) + [norm_1(w(i)) - transpose(nu(i)) * (D*u - w(i)) + beta/2 * norm_2(D*u - w(i))] - [transpose(lambda) * (A*u - Y) + mu/2 * norm_2(A*u - Y)];
			}
		}
	} while (Lagrangian > tol)


    BoostDoubleMatrix RecoveredImage(32,32);  // Create a dummy matrix to return

return RecoveredImage;
}