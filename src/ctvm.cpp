#include "ctvm.h"

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, 
                                       BoostDoubleVector TiltAngles)
{
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
	unsigned double alpha =
	unsigned double Lagrangian; // ?
	float alpha = 1.05;
	float innerstop, outerstop;
	float innertol = 0.5, outertol = 0.5, tol = 0.5; // multiple tol
	float delta, rho, eta;
	
	BoostDoubleMatrix W (N, 2); // w(i) = 0 for all i
	BoostDoubleMatrix NU (N, 2);
	BoostDoubleMatrix A (m, n);
	
	BoostDoubleVector Y (m);
	BoostDoubleVector X (n);
	
	A = CreateRandomMatrix(m, n);
	Y = MatrixToVector(Sinogram); // u(0) ?
	
	for (unsigned int i = 0; i < 2*n; ++i)
	{
		do /********************************** TVAL3 SCHEME **********************************/
		{
			for (unsigned int k = 0; k < m; ++k)
			{
				w(k+1) = w(k);3
				u(k+1) = u(k);
				for (unsigned int j = 0; j < m; ++j)
				{
					delta = 0.5;
					rho = 0.5;
					eta = 0.5;
					w(0) = // ?
					u(0) = // ?
					
					do /****************** ALTERNATING MINIMISATION SCHEME ******************/
					{
						do /* "w sub-problem" */
						{
							alpha(j) = rho * alpha(j);
						} while (/*(2.33)*/);          // ...
						  /* "u sub-problem" */
						innerstop = norm_2(u(j+1) - u(j));
					} while (innerstop > tol);
				}
				// Lagrangian ?
				//		L(w(k+1),u(k+1)) = L(w(k),u(k)) + [norm_1(w(i)) - transpose(nu(i)) * (D*u - w(i)) + beta/2 * norm_2(D*u - w(i))] - [transpose(lambda) * (A*u - Y) + mu/2 * norm_2(A*u - Y)];
			}
			outerstop = norm_2(u(k+1) - u(k));
		} while (outerstop > tol);
	}

    BoostDoubleMatrix RecoveredImage(32,32);  // Create a dummy matrix to return

return RecoveredImage;
}

DoubleBoostVector Gradient2D(DoubleBoostVector U, unsigned double pixel) // pixel = actual pixel number -1 ??
{
	unsigned double n = U.size ();
	unsigned double l = sqrt(n);
	unsigned double rows = 0;
	unsigned double cols = 0;
	float q = 0;
	
	DoubleBoostVector Di (2);
	DoubleBoostMatrix X (l, l);
	
	X = VectorToMatrix(U, l, l);

	for (unsigned double i = 0; i < l; ++i)
	{
		if (!pixel)
		{
			rows = i;
			cols = 0;
			break;
		}
		else
		{
			for (unsigned double j = 0; j < l; ++j)
			{
				q = j / pixel;
				if (q)
				{
					rows = i;
					cols = j;
					break;
				}
			}
		pixel = pixel - l;
	}
	
	
return Di;
}