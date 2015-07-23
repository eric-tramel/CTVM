#include "ctvm.h"

BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long rank) // rank = actual rank number -1! (first rank = U(0))
{
	unsigned long n = U.size();
	unsigned long pixel = rank;
	unsigned long l = sqrt(n); // !
	unsigned long rows = 0;
	unsigned long cols = 0;
	float q = 0;

	BoostDoubleVector Di (2);
	BoostDoubleMatrix X (l, l);

	X = VectorToMatrix(U, l, l);

	for (unsigned long i = 0; i < l; ++i) /* Find pixel place in the matrix */ 
	{
		if (!pixel)
		{
			rows = i;
			cols = 0;
			break;
		}
		else
		{
			for (unsigned long j = 0; j < l; ++j)
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
	}
	if (cols == l-1 || rows == l-1)
	{
		Di(0) = X(rows, cols) - X(rows, cols);
		Di(1) = X(rows, cols) - X(rows, cols);
	}
	else if (cols == l-1)
	{
		Di(0) = X(rows, cols) - X(rows, cols);
		Di(1) = X(rows + 1, cols) - X(rows, cols);
	}
	else if (rows == l-1)
	{
		Di(0) = X(rows, cols + 1) - X(rows, cols);
		Di(1) = X(rows, cols) - X(rows, cols);
	}
	else
	{
		Di(0) = X(rows, cols + 1) - X(rows, cols);
		Di(1) = X(rows + 1, cols) - X(rows, cols);
	}

	for (unsigned long i = 0; i < l; ++i) /* Find the down gradient */
	{
		for (unsigned long i = 0; i < l; ++i) /* Find the right gradient */
		{
		}
	}
return Di;
}

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles)
{
	/*
	Function: tval3_reconstruction
	------------------------------
	Perform two-dimensional TVAL3 tomographic reconstruction given a
	set of measurements (Sinogram) and the measured tilt angles.
	*/

	/* <TODO: Everything!> */

	unsigned long l = Sinogram.size1(); /* Size of the sample (in pixels) */
	unsigned long o = Sinogram.size2(); /* Numbers of tilt angles */
	unsigned long m = l * o; /* Numbers of measurements */
	unsigned long n = l * l;
	double lambda = 1;
	double beta = 1;
	double mu = 1;
	double Lagrangian; // ?
	double alpha = 1.05;
	double innerstop, outerstop;
	double innertol = 0.5, outertol = 0.5, tol = 0.5; // multiple tol
	double delta, rho, eta;

	BoostDoubleMatrix W(n, 2); // w(i) = 0 for all i
	BoostDoubleMatrix NU(n, 2);
	BoostDoubleMatrix A(m, n);

	BoostDoubleVector ALPHA(); // ?
	BoostDoubleVector Y(m);
	BoostDoubleVector X(n);

	A = CreateRandomMatrix(m, n);
	Y = MatrixToVector(Sinogram); // u(0) ?

	for (unsigned int i = 0; i < 2 * n; ++i)
	{
		do /********************************** TVAL3 SCHEME **********************************/
		{
			for (unsigned int k = 0; k < m; ++k)
			{
				w(k + 1) = w(k);
				u(k + 1) = u(k);
				for (unsigned int j = 0; j < m; ++j)
				{
					delta = 0.5;
					rho = 0.5;
					eta = 0.5;
					w(0) = ; // ?
					u(0) = ; // ?

					do /****************** ALTERNATING MINIMISATION SCHEME ******************/
					{
						do /* "w sub-problem" */
						{
							ALPHA(j) = rho * ALPHA(j);
						} while (/*(2.33)*/);          // ...
													   /* "u sub-problem" */
						innerstop = norm_2(u(j + 1) - u(j));
					} while (innerstop > tol);
				}
				// Lagrangian ?
				//		L(w(k+1),u(k+1)) = L(w(k),u(k)) + [norm_1(w(i)) - transpose(nu(i)) * (D*u - w(i)) + beta/2 * norm_2(D*u - w(i))] - [transpose(lambda) * (A*u - Y) + mu/2 * norm_2(A*u - Y)];
			}
			outerstop = norm_2(u(k + 1) - u(k));
		} while (outerstop > tol);
	}

	BoostDoubleMatrix RecoveredImage(32, 32);  // Create a dummy matrix to return

	return RecoveredImage;
}