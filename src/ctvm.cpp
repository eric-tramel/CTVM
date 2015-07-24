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

/****************** Find pixel place in the matrix *******************/
	for (unsigned long i = 0; i < l; ++i) 
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

/******************** Find gradient at the rank i ********************/
	if (cols == l - 1 && rows == l - 1)
	{
		Di(0) = X(rows, cols) - X(rows, cols);
		Di(1) = X(rows, cols) - X(rows, cols);
	}
	else if (cols == l - 1)
	{
		Di(0) = X(rows, cols) - X(rows, cols);
		Di(1) = X(rows + 1, cols) - X(rows, cols);
	}
	else if (rows == l - 1)
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
/******************* List all gradients for all i *******************/
	unsigned long n = U.size();
	BoostDoubleVector Di(2);
	BoostDoubleMatrix Du(n, 2);

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

	for (unsigned long i = 0; i < n; ++i)
	{
		Di = Gradient2D(U, i);
		for (unsigned long j = 0; j < 2; ++j)
		{
			Wi(j) = W(i, j);
			NUi(j) = NU(i, j);
		}
		BoostDoubleVector DIFFi = Di*U(i) - Wi;
		double norm_wi = norm_2(Wi);
		double square_norm_diffi = norm_2(DIFFi)*norm_2(DIFFi);

		L = L + norm_wi - inner_prod(NUi, DIFFi) + (beta / 2)*square_norm_diffi;
	}
	BoostDoubleVector PROD(m);
	BoostDoubleVector DIFF(m);
	PROD = prod(A, U);
	DIFF = PROD - B;
	double square_norm_diff = norm_2(DIFF)*norm_2(DIFF);

	L = L - inner_prod(LAMBDA, DIFF) + (mu/2)*square_norm_diff;
return L;
}

/*BoostDoubleVector alternating_minimisation(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu, unsigned long rank)
{
unsigned long i = rank;
double n = U.size();
double delta = 0.5;
double rho = 0.5;
double eta = 0.5;
double tol = 0.01;
double innerstop;

BoostDoubleVector ALT_MIN (2);
BoostDoubleVector Q (n);
BoostDoubleVector C (n);
BoostDoubleVector D (n);
BoostDoubleVector ALPHA (n);

ALT_MIN(0) = 1; ALT_MIN(1) = 1; // ALT_MIN(0) = w(i,0), ALT_MIN(1) = u(0)
Q(0) = 1;
C(0) = Lagrangian(A, U, B, W, NU, LAMBDA, beta, mu); // NOT OK, Lagrangian function for all values and not the initialisation

double uj_1 = ALT_MIN(1);

do
{
//******************** "w sub-problem" ********************
for (unsigned long j = 0; j < n - 1; ++j)
{

ALT_MIN(0) = W(j+1);
do
{
ALPHA(j) = rho * ALPHA(j);
} while (/*(2.33)*//*);          // NOT OK

//******************** "u sub-problem" ********************
D(0) = onestep_gradient(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu, unsigned long rank, unsigned long j);
D(1) = ;
double dj_1 = D(0);
double dj = D(1);

U(j + 1) = U(j) - ALPHA(j) * d(j);

innerstop = norm_2(U(j + 1) - U(j));
uj_1 = ALT_MIN(1); // uj_1 = U(j-1)
ALT_MIN(1) = U(j + 1);
}
} while (innerstop > tol);
return ALT_MIN;
}*/

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
