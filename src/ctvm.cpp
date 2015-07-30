#include "ctvm.h"

BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long rank) // rank = actual rank number -1 (first rank = U(0))
{
	/*
	* Function: Gradient2D (Di*u)
	* ----------------------------
	* Take a vector and give the right gradient
	* and the down gradient of the rank i 
	*
	*/
	unsigned long n = U.size();
	unsigned long pixel = rank;
	unsigned long l = sqrt(n); // !
	unsigned long rows = 0;
	unsigned long cols = 0;
	double q = 0;

	BoostDoubleVector Du (2);
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
		Du(0) = X(rows, cols) - X(rows, cols);
		Du(1) = X(rows, cols) - X(rows, cols);
	}
	else if (cols == l - 1)
	{
		Du(0) = X(rows, cols) - X(rows, cols);
		Du(1) = X(rows + 1, cols) - X(rows, cols);
	}
	else if (rows == l - 1)
	{
		Du(0) = X(rows, cols + 1) - X(rows, cols);
		Du(1) = X(rows, cols) - X(rows, cols);
	}
	else
	{
		Du(0) = X(rows, cols + 1) - X(rows, cols);
		Du(1) = X(rows + 1, cols) - X(rows, cols);
	}

	return Du;
}

BoostDoubleMatrix Gradient2DMatrix(BoostDoubleVector U)
{
	/*
	* Function: Gradient2DMatrix
	* ----------------------------
	* Give the matrix of all the right and down 
	* gradients for all the values of the vector
	*
	*/
	unsigned long n = U.size();
	BoostDoubleVector Du (2);
	BoostDoubleMatrix D (n, 2);

	for (unsigned long i = 0; i < n; ++i)
	{
		Du = Gradient2D(U, i);
		for (int j = 0; j < 2; ++j)
		{
			D(i,j) = Du(j);
		}
	}
	return D;
}

double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{
/*
* Function: Lagrangian
* ----------------------------
* Give the result (double) of the Augmented Lagrangian function 
* developped by Chengbo Li in his thesis "An efficient 
* algorithm for total variation regularization with 
* applications to the single pixel camera and compressive 
* sensing" 
*
*/
	double L = 0;
	unsigned long n = U.size();
	BoostDoubleVector Wi (2);
	BoostDoubleVector DiU (2);
	BoostDoubleVector NUi (2);

	for (unsigned long i = 0; i < n; ++i)
	{
		DiU = Gradient2D(U, i);
		for (int j = 0; j < 2; ++j)
		{
			Wi(j) = W(i, j);
			NUi(j) = NU(i, j);
		}
		BoostDoubleVector DIFFi = DiU - Wi;
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


BoostDoubleVector ShrikeIsotropic(BoostDoubleVector W, BoostDoubleVector Nu, double beta){
/* 
 * Function: Shrike Isotropic
 * -------------------------------
 * Impelements the Isotropic version of the shrinkage function to be applied
 * to the gradient vector at pixel `i`. Takes as input an "un-shrunk" gradient
 * vector and returns the "shriked" version.
 *
 * Input --
 *  W:    a (d x 1) gradient vector (i.e. for 2-d images, d=2)
 *  Nu:   a (d x 1) set of multipliers
 *  beta: a scalar scaling term
 *
 * Output -- a (d x 1) "shriked" version of the gradient vector
*/
 	unsigned int d = W.size();
 	BoostDoubleVector WShifted = W - Nu/beta;
	BoostDoubleVector WShriked = AbsoluteValueVector(WShifted) - BoostScalarDoubleVector(d,1/beta); 	

return HadamardProduct(MaxVector(WShriked,0.0),SignVector(WShifted));
}

BoostDoubleVector ShrikeAnisotropic(BoostDoubleVector W, BoostDoubleVector Nu, double beta){
/* 
 * Function: Shrike Anisotropic
 * -------------------------------
 * Impelements the Anisotropic version of the shrinkage function to be applied
 * to the gradient vector at pixel `i`. Takes as input an "un-shrunk" gradient
 * vector and returns the "shriked" version.
 *
 * Input --
 *  W:    a (d x 1) gradient vector (i.e. for 2-d images, d=2)
 *  Nu:   a (d x 1) set of multipliers
 *  beta: a scalar scaling term
 *
 * Output -- a (d x 1) "shriked" version of the gradient vector
*/
	BoostDoubleVector WShifted = W - Nu/beta;
	double WShiftedNorm = norm_2(WShifted);

return (fmax(WShiftedNorm - 1/beta,0.0) * (WShifted/WShiftedNorm));
}


BoostDoubleVector Onestep_direction(BoostDoubleMatrix A, BoostDoubleVector Uk, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{
	unsigned long n = Uk.size();
	BoostDoubleVector D (n);
	BoostDoubleMatrix Di = Gradient2DMatrix(Uk);
	BoostDoubleVector Wi (2);
	BoostDoubleVector NUi (2);

	for (unsigned long i = 0; i < n; ++i)
	{
		BoostDoubleVector DiU = Gradient2D(Uk, i); // Is there more pertinent method  than calculate this gradient at every loop?
		for (int j = 0; j < 2; ++j)
		{
			Wi(j) = W(i, j);
			NUi(j) = NU(i, j);
		}
		BoostDoubleVector DiU_W = -DiU - Wi;
		D = beta*prod(Di, DiU_W) - prod(Di, NUi); // Di = (2x1)Vector -> innerproduct NOT OK
	}

	BoostDoubleVector Au = prod(A, Uk);
	BoostDoubleVector DIFF = Au - B;

	D = D + mu*prod(A, DIFF) - prod(A, LAMBDA);

return D;
}

double U_subfunction(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{
	double Q;
	unsigned long n = U.size();
	BoostDoubleVector NUi (2);
	BoostDoubleVector Wi (2);

	for (unsigned long i = 0; i < n; ++i)
	{
		BoostDoubleVector DiU = Gradient2D(U, i);
		for (int j = 0; j < n; ++j)
		{
			NUi(j) = NU(i, j);
			Wi(j) = W(i, j);
		}
		BoostDoubleVector DIFFk = DiU - Wi;
		double square_norm_diffk = norm_2(DIFFk)*norm_2(DIFFk);

		Q = Q - inner_prod(NUi, DIFFk) + (beta / 2) * square_norm_diffk;
	}
	BoostDoubleVector PROD = prod(A, U);
	BoostDoubleVector DIFF = PROD - B;
	double square_norm_diff = norm_2(DIFF)*norm_2(DIFF);

	Q = Q - inner_prod(LAMBDA, DIFF) + (mu / 2)*square_norm_diff;
	return Q;
}

BoostDoubleMatrix alternating_minimisation(BoostDoubleMatrix A, BoostDoubleVector U, BoostDoubleVector B, BoostDoubleMatrix W, BoostDoubleMatrix NU, BoostDoubleVector LAMBDA, double beta, double mu)
{
	double n = U.size();
	double delta = 0.5;
	double rho = 0.5;
	double eta = 0.5;
	double Pk = 1;
	double C = Lagrangian(A, U, B, W, NU, LAMBDA, beta, mu);

	double armijo_tol, Qk, innerstop;
	double tol = 0.01;

	BoostDoubleMatrix AL_MIN (n, 3); // AL_MIN(:,0) and AL_MIN(:,1) = W(k+1), AL_MIN(:,2) = U(k+1)
	BoostDoubleVector Uk_1 (n);
	BoostDoubleVector Uk = U;
	do
	{
//*************************** "w sub-problem" ***************************
		for (unsigned long i = 0; i < n; ++i)
		{
			BoostDoubleVector DiUk = Gradient2D(Uk, i);
			BoostDoubleVector NUi (2);
			for (int j = 0; j < 2; ++j)
			{
				NUi(j) = NU(i, j);
				BoostDoubleVector Wi = ShrikeAnisotropic(DiUk, NUi, beta);
				W(j, i) = Wi(j);
			}
		}
//*************************** "u sub-problem" ***************************
		BoostDoubleVector Sk = Uk - Uk_1;
		BoostDoubleVector Dk_1 = Onestep_direction(A, Uk_1, B, W, NU, LAMBDA, beta, mu);
		BoostDoubleVector Dk = Onestep_direction(A, Uk, B, W, NU, LAMBDA, beta, mu);
		BoostDoubleVector Yk = Dk - Dk_1;

		//******** alpha = onestep_gradient ********
		double alpha = inner_prod(Sk, Yk) / inner_prod(Yk, Yk);
		do 
		{ 
			alpha = rho * alpha;

			BoostDoubleVector Uk_alphaD = Uk - alpha * Dk;
			Qk = U_subfunction(A, Uk_alphaD, B, W, NU, LAMBDA, beta, mu);
			armijo_tol = C - delta*alpha*inner_prod(Dk, Dk);
		} while (Qk > armijo_tol);

		BoostDoubleVector Uk1 = Uk - alpha * Dk;
		innerstop = norm_2(Uk1 - Uk);

//************************ Implement coefficents ************************
		double Pk1 = eta*Pk + 1;
		double Qk1 = U_subfunction(A, Uk1, B, W, NU, LAMBDA, beta, mu);
		C = (eta*Pk*C + Qk1)/Pk1;

		Pk = Pk1;
		Uk_1 = Uk;
		Uk = Uk1;
	} while (innerstop > tol);
	
	for (unsigned long pix = 0; pix < n; ++pix)
	{
		AL_MIN(pix, 0) = W(pix, 0);
		AL_MIN(pix, 1) = W(pix, 1);
		AL_MIN(pix, 2) = Uk(pix);
	}
return AL_MIN;
}

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix Sinogram, BoostDoubleVector TiltAngles) // Why TiltAngles?
{
	unsigned long l = Sinogram.size1(); // Size of the sample (in pixels)
	unsigned long o = Sinogram.size2(); // Numbers of tilt angles
	unsigned long m = l * o; // Numbers of measurements
	unsigned long n = l * l;
	
	double mu = 3;
	double beta = sqrt(2);
	double coef = 1.05;
	double outerstop;
	double tol = 0.01;

	BoostDoubleMatrix W (n, 2); // W(i,0) = 0 for all i
	BoostDoubleMatrix Wk (n, 2);
	BoostDoubleMatrix NU (n, 2);
	BoostDoubleMatrix A (m, n);
	BoostDoubleMatrix MIN (n, 3);
	BoostDoubleMatrix X (l, l);
	
	BoostDoubleVector LAMBDA (m);
	BoostDoubleVector B (m);
	BoostDoubleVector U (n); // U(0) = 0 for all i
	BoostDoubleVector Uk (n);

	A = CreateRandomMatrix(m, n);
	B = MatrixToVector(Sinogram); // u(0) ?
	
	do
	{
		Wk = W;
		Uk = U;
		MIN = alternating_minimisation(A, Uk, B, Wk, NU, LAMBDA, beta, mu);
		for (unsigned long i = 0; i < n; ++i)
		{
			W(i, 0) = MIN(i, 0);
			W(i, 1) = MIN(i, 1);
			U(i) = MIN(i, 2);
		}
		BoostDoubleMatrix DiU = Gradient2DMatrix(U);
		NU = NU - beta*(DiU - W);
		LAMBDA = LAMBDA - mu*(prod(A, U) - B);

		beta = coef*beta;
		mu = coef*beta;

		outerstop = norm_2(U - Uk);
	} while (outerstop > tol);

	X = VectorToMatrix(U, l, l);
return X;
}
