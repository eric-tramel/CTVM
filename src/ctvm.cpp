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
	
	if (pixel >= n)
	{
		std::cout << "WARNING: Gradient's rank bigger than specimen size" << std::endl;
		exit(EXIT_FAILURE);
	}

	for (unsigned long i = 0; i < l; ++i) /****** Find pixel place in the matrix *******/ 
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

	if (cols == l - 1 && rows == l - 1) /******* Find gradient at the rank i *******/
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

return Di;
}

BoostDoubleMatrix Gradient2DMatrix(BoostDoubleVector U)
{
	/* List all gradients for all i */

	unsigned long n = U.size();
	BoostDoubleVector Di (2);
	BoostDoubleMatrix Du (n,2);

	for (unsigned long i = 0; i < n; ++i)
	{
		Di = Gradient2D(U, i);
		for (int j = 0; j < 2; ++j)
		{
			Du(i, j) = Di(j);
		}
	}
return Du;
}

double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{
	double L = 0;
	unsigned long n = U.size();
	unsigned long m = B.size();
	BoostDoubleVector Wi (2);
	BoostDoubleVector Di (2);
	BoostDoubleVector NUi (2);
	BoostDoubleVector PROD (m);
	BoostDoubleVector DIFF(m);
	PROD = prod(A, U);
	DIFF = PROD - B;
	double norm_diff = norm_2(DIFF);

	for (unsigned long i = 0; i < n; ++i)
	{
		Di = Gradient2D(U, i);
		for (unsigned long j = 0; j < 2; ++j)
		{
			Wi(j) = W(i,j);
			NUi(j) = NU(i,j);
		}
		BoostDoubleVector tNUi = trans(NUi);
		BoostDoubleVector DIFFi = Di*U(i) - Wi;
		double norm_wi = norm_2(Wi);
		double norm_diffi = norm_2(DIFFi);

		L = L + norm_wi - prod(tNUi, DIFFi) + (beta/2)*norm_diffi; // VECOTRS PRODUCT ON UBLAS??
	}
	L = L - prod(LAMBDA, DIFF) + (mu/2)*norm_diff;
return L;
}

/*BoostDoubleVector alternating_minimisation(...)
{
	BoostDoubleVector MIN (2);

	do
	{
		do // "w sub-problem"
		{
			ALPHA(j) = rho * ALPHA(j);
		} while (*//*(2.33)*//*);          // ...
									   // "u sub-problem"
		innerstop = norm_2(u(j + 1) - u(j));
	} while (innerstop > tol);

return MIN;
} */

/*BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles)
{
	unsigned long l = Sinogram.size1(); // Size of the sample (in pixels)
	unsigned long o = Sinogram.size2(); // Numbers of tilt angles
	unsigned long m = l * o; // Numbers of measurements
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
	BoostDoubleMatrix Z(l, l);

	BoostDoubleVector ALPHA(); // ?
	BoostDoubleVector Y(m);
	BoostDoubleVector X(n);

	A = CreateRandomMatrix(m, n);
	Y = MatrixToVector(Sinogram); // u(0) ?

	for (unsigned int i = 0; i < 2 * n; ++i)
	{
		do
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
				}
			}
			outerstop = norm_2(u(k + 1) - u(k));
		} while (outerstop > tol);
	}
	return Z;
} */
