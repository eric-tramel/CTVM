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
		for (int j = 0; j < 2; ++j)
		{
			Wi(j) = W(i, j);
			NUi(j) = NU(i, j);
		}
		BoostDoubleVector DIFFi = Di*U(i) - Wi;
		double norm_wi = norm_2(Wi); // norm_1 for anisotropic TV or norm_2 for isotropic TV
		double square_norm_diffi = norm_2(DIFFi)*norm_2(DIFFi);

		L = L + norm_wi - inner_prod(NUi, DIFFi) + (beta/2) * square_norm_diffi;
	}
	BoostDoubleVector PROD = prod(A, U);
	BoostDoubleVector DIFF = PROD - B;
	double square_norm_diff = norm_2(DIFF)*norm_2(DIFF);

	L = L - inner_prod(LAMBDA, DIFF) + (mu/2)*square_norm_diff;
return L;
}

BoostDoubleVector Shrike(BoostDoubleVector DiUk, BoostDoubleVector NUi, double beta)
{
	BoostDoubleVector W (2);

	BoostDoubleVector DIFF = DiUk - NUi/beta;
	double norm_diff = norm_2(DIFF);
	double x = norm_diff - 1 / beta;
	if (x < 0) x = 0;
	W = x * (DIFF / norm_diff);

return W;
}

/*double onestep_direction(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{

}*/

double u_subfunction(BoostDoubleMatrix A, double u, BoostDoubleVector B, BoostDoubleVector Dl, BoostDoubleVector Wl, BoostDoubleVector NUl, BoostDoubleVector LAMBDA, double beta, double mu)
{
	double Q;
	unsigned long n = U.size();
	unsigned long m = B.size();

	for (unsigned long i = 0; i < n; ++i)
	{
		BoostDoubleVector DIFFl = Dl*u - Wl;
		double square_norm_diffl = norm_2(DIFFl)*norm_2(DIFFl);

		Q = Q - inner_prod(NUl, DIFFl) + (beta / 2) * square_norm_diffl;
	}
	BoostDoubleVector PROD = prod(A, u); // u value or U vector?
	BoostDoubleVector DIFF = PROD - B;
	double square_norm_diff = norm_2(DIFF)*norm_2(DIFF);

	Q = Q - inner_prod(LAMBDA, DIFF) + (mu / 2)*square_norm_diff;
	return Q;
}

BoostDoubleVector alternating_minimisation(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{
	double n = U.size();
	double delta = 0.5;
	double rho = 0.5;
	double eta = 0.5;
	double Pl = 1;
	double C = Lagrangian(A, U, B, W, NU, LAMBDA, beta, mu); // NOT OK, Lagrangian function for all values and not the initialisation

	double tol = 0.01;

	BoostDoubleVector ALT_MIN (3); // ALT_MIN(0) and ALT_MIN(1) = w(i,0), ALT_MIN(2) = u(0)
	ALT_MIN(0) = 1; ALT_MIN(1) = 1, ALT_MIN(2) = 1;
	W(0, 0) = ALT_MIN(0); W(0, 1) = ALT_MIN(1);
	U(0) = ALT_MIN(2);
	U(1) = ; // Initialisation?
	do
	{
//*************************** "w sub-problem" ***************************
		for (unsigned long l = 1; l < n-1; ++l)
		{
			BoostDoubleVector Dl = Gradient2D(U, l);
			BoostDoubleVector DlUl = Dl * U(l);
			BoostDoubleVector NUl (2);
			BoostDoubleVector Wl (2);
			for (int j = 0; j < 2; ++j)
			{
				Wl(j) = W(l, j);
				NUl(j) = NU(l, j); 
			}
			Wl = Shrike(DlUl, NUl, beta); // NOT OK (how to fill W vector with Wl ?)
			ALT_MIN(0) = Wl(0);
			ALT_MIN(1) = Wl(1);

			double sl = U(l) - U(l-1);
			double dl = onestep_direction(A, U(l), B, Wl, NUl, LAMBDA, beta, mu);
			double dl_1 = onestep_direction(A, U(l-1), B, Wl, NUl, LAMBDA, beta, mu);
			double yl = dl - dl_1;

			//******** alpha = onestep_gradient ********
			double alpha = (sl*yl) / (yl*yl);
			do 
			{ 
				alpha = rho * alpha;
				double u_alphad = U(l) - alpha * dl;
				double Ql = u_subfunction(A, u_alphad, B, Dl, Wl, NUl, LAMBDA, beta, mu);
				double armijo_tol = C - delta*alpha*dl*dl; // trans(dl)? cpx?
			} while (Ql > armijo_tol);

//*************************** "u sub-problem" ***************************

			U(l+1) = U(l) - alpha * dl;
			ALT_MIN(2) = U(l+1);
			double innerstop = sqrt((U(l+1) - U(l))*(U(l+1) - U(l)));

			//********** Implement coefficents **********
			double Pl1 = eta*Pl + 1;
			double Ql1 = u_subfunction(A, U(l+1), B, Dl, Wl, NUl, LAMBDA, beta, mu);
			C = (eta*Pl*C + Ql1)/Pl1;
			Pl = Pl1;
		}
	} while (innerstop > tol);
return ALT_MIN;
}

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
