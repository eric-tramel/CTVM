#include "ctvm.h"


BoostDoubleVector PixelGradient(BoostDoubleVector X, unsigned long Index, 
								unsigned long SideLength){
/*
* Function: PixelGradient
* -----------------------
* Given a rasterized image vector, calculate the gradient vector at the
* specified index. As we are assuming images, this function is used to
* calculate the *two-dimensional* gradient vector.
*
* Input --
* X: a (N x 1) vector representing a rasterized image
* Index: the pixel at which to calculate the gradient
* SideLength: assuming square image dimensions, the length of the image side.
*             I.e. N = SideLength^2.
*
* Output -- A (2 x 1) pixel gradient vector.
*/
	BoostDoubleVector Gradient = BoostZeroVector(2);

	// Horizontal Gradient Component
	int RightIndex = RightNeighbor(Index,SideLength);
	Gradient(HORZ) = (RightIndex>0) ? (X[Index] - X[RightIndex]) : 0.0;

	// Vertical Gradient Component
	int DownIndex = DownNeighbor(Index,SideLength);
	Gradient(VERT) = (DownIndex>0) ? (X[Index] - X[DownIndex]) : 0.0;

return Gradient;
}

BoostDoubleMatrix AllPixelGradients(BoostDoubleVector X, unsigned long SideLength){
/*
* Function: AllPixelGradients
* ---------------------------
* Given a rasterized image vector, calculate the gradient vectors at every pixel
* and return the entire set of gradients as a matrix
*
* Input --
* X: a (N x 1) vector representing a rasterized image
* SideLength: assuming square image dimensions, the length of the image side.
*             I.e. N = SideLength^2.
*
* Output -- A (N x 2) pixel gradient vector.
*/	
	unsigned long N = X.size();
	BoostDoubleMatrix AllGradients(N,2);

	for(unsigned long i = 0; i < N; ++i){
		SetRow(AllGradients,PixelGradient(X,i,SideLength),i);
	}

return AllGradients;
}

BoostDoubleVector PixelGradientAdjointSum(BoostDoubleMatrix G, unsigned long SideLength){
/*
* Function: PixelGradientAdjointSum
* ---------------------------------
* Given a set of gradients at each pixel, calculate the adjoint sum, which is the 
* adjoint of the gradient operation summed over each of the pixels.
*
*    					X = sum_{i=1:N} D_i^T * G_i
*
* where D_i represents the gradient matrix at pixel i, T is the transpose operator, and
* G_i is the gradient vectore at pixel i. In other words, this function is a map from the
* space of pixel gradients back to the image space.
*
* Input --
* G: a (N x 2) matrix of pixel gradients
* SideLength: assuming square image dimensions, the length of the image side.
*             I.e. N = SideLength^2.
*
* Output -- A (N x 1) rasterized image vector.
*/	
	unsigned long N = G.size1();
	BoostDoubleVector ImageVector = BoostZeroVector(N);
	for(unsigned long i = 0; i < N; ++i){
		int thisRightNeighbor = RightNeighbor(i,SideLength);
		int thisDownNeighbor = DownNeighbor(i,SideLength);
		BoostDoubleVector thisGradient = GetRow(G,i);

		ImageVector(i) += thisGradient(HORZ) + thisGradient(VERT);
		if(thisRightNeighbor > 0){
			ImageVector(thisRightNeighbor) += -1* thisGradient(HORZ);
		}
		if(thisDownNeighbor > 0){
			ImageVector(thisDownNeighbor) += -1* thisGradient(VERT);
		}
	}

return ImageVector;
}

BoostDoubleVector ShrikeAnisotropic(BoostDoubleVector W, BoostDoubleVector Nu,
									double beta){
/* 
 * Function: Shrike Anisotropic
 * ----------------------------
 * Implements the Anisotropic version of the shrinkage function to be applied
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
 	unsigned long d = W.size();
 	BoostDoubleVector WShifted = W - Nu/beta;
	BoostDoubleVector WShriked = AbsoluteValueVector(WShifted) - BoostScalarDoubleVector(d,1/beta); 	
return HadamardProduct(MaxVector(WShriked,0.0),SignVector(WShifted));
}

BoostDoubleVector ShrikeIsotropic(BoostDoubleVector W, BoostDoubleVector Nu,
								  double beta){
/* 
 * Function: Shrike Isotropic
 * -------------------------------
 * Implements the Isotropic version of the shrinkage function to be applied
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
	BoostDoubleVector result;

	if(WShiftedNorm < 0.000000000001){
		result = BoostZeroVector(2);
	}
	else{
		result = (fmax(WShiftedNorm - 1/beta,0.0) * (WShifted/WShiftedNorm));
	}

return result;
}

BoostDoubleMatrix ApplyShrike(BoostDoubleMatrix AllW, BoostDoubleMatrix AllNu,
							  double beta, TVType ShrikeMode){
/* 
 * Function: ApplyShrike
 * ---------------------
 * Applies the Shrinkage-like operator to every gradient in W.
 *
 * Input --
 *  W:    a (N x d) matrix of gradients (i.e. for 2-d images, d=2)
 *  Nu:   a (N x d) set of multipliers
 *  beta: a scalar scaling term
 *  ShrikeMode: an Enum value of TVType whic specifies whether we use the
 *              isotropic or anisotropic Shrike operator.
 *
 * Output -- a (N x d) "shriked" version of the gradient matrix
*/	
 	/* Problem Dimensions */
 	unsigned long N = AllW.size1();

 	for(unsigned long i = 0; i < N; ++i){
 		BoostDoubleVector thisW = GetRow(AllW,i);
 		BoostDoubleVector thisNu = GetRow(AllNu,i);
 		BoostDoubleVector thisWShriked;

 		switch(ShrikeMode){
 			case ISOTROPIC:
 				thisWShriked = ShrikeIsotropic(thisW,thisNu,beta);
 				break;
 			case ANISOTROPIC:
 				thisWShriked = ShrikeAnisotropic(thisW,thisNu,beta);
 				break;
 		} 	
 		SetRow(AllW,thisWShriked,i);
 	}
return AllW;
}

double Lagrangian(BoostDoubleMatrix A, BoostDoubleVector U, 
					 BoostDoubleVector B, BoostDoubleMatrix W, 
					 BoostDoubleMatrix Nu, BoostDoubleVector Lambda, 
					 double beta, double mu,
					 unsigned long SideLength, TVType GradNorm){
/*
* Function: Lagrangian
* --------------------
* Calculate the Lagrangian cost function for the give problem state.
*
* Input --
* A: an (M x N) projection matrix
* U: an (N x 1) vector representing a rasterized image prediction
* B: an (M x 1) set of observations
* W: an (N x 2) set of dual variables corresponding to per-pixel (-voxel) gradients
* Nu: an (N x 2) set of Lagrangian multipliers
* Lambda: an (M x 1) set of Lagrangian multiplies
* beta: scaling term on the matching between W and the true gradients
* mu: scaling term on the matching between A*u and b
* SideLength: the side length for the target image, i.e. N = SideLength^2
* GradNorm: which TV norm to use on the gradients (Iso- or Anisotropic)
*
* Output -- a decimal value for the cost.
*/
	double L = 0.0;

	// Get all Gradients
	BoostDoubleMatrix Du = AllPixelGradients(U,SideLength);

	// Loop over pixels
	BoostDoubleVector Dui;
	BoostDoubleVector Wi;
	BoostDoubleVector Nui;
	for (unsigned long i = 0; i < U.size(); ++i){
		Dui = GetRow(Du,i);
		Wi  = GetRow(W,i);
		Nui = GetRow(Nu,i);

		BoostDoubleVector GradDiff = Dui-Wi;
		L+= - inner_prod(Nui, GradDiff) + (beta/2) * SquareNorm(GradDiff);

		switch(GradNorm){
			case ISOTROPIC:
				L += norm_2(Wi);
				break;
			case ANISOTROPIC:
				L += norm_1(Wi);
				break;
		}
	}

	// Residual Contribution
	BoostDoubleVector Residual = prod(A,U) - B;
	L += -inner_prod(Lambda, Residual) + (mu/2)*SquareNorm(Residual);

return L;
}

BoostDoubleVector Onestep_Direction(BoostDoubleMatrix A, BoostDoubleVector U, 
									BoostDoubleVector B, BoostDoubleMatrix W, 
									BoostDoubleMatrix Nu, BoostDoubleVector Lambda, 
									double beta, double mu,
									unsigned long SideLength){
	/*
	* Function: Onestep_Direction
	* ---------------------------
	* Calculate the direction for the gradient in the one-step steepest descent method.
	* To understand this following command, one can separate it in two steps
	* (that is, the sum over all the pixels) Du = AllPixelGradients(U,SideLength)
	*										 Dk = -PixelGradientAdjointSum(beta*Du + beta*W + Nu) + mu*A'*(A*u - b) - A'*lambda
	*
	* Input --
	* A: an (M x N) projection matrix
	* U: an (N x 1) vector representing a rasterized image prediction
	* B: an (M x 1) set of observations
	* W: an (N x 2) set of dual variables corresponding to per-pixel (-voxel) gradients
	* Nu: an (N x 2) set of Lagrangian multipliers
	* Lambda: an (M x 1) set of Lagrangian multiplies
	* beta: scaling term on the matching between W and the true gradients
	* mu: scaling term on the matching between A*u and b
	* SideLength: the side length for the target image, i.e. N = SideLength^2
	*
	* Output -- a decimal value for the cost.
	*/
	BoostDoubleVector Dk = BoostZeroVector(U.size());

	// BoostDoubleVector NegAdjSum = -PixelGradientAdjointSum(beta*AllPixelGradients(U, SideLength) + beta*W + NU, SideLength);
	// BoostDoubleVector TermTwo = mu*prod(trans(A), prod(A, U) - B);
	// BoostDoubleVector TermThree = - prod(trans(A), LAMBDA);
	// BoostDoubleMatrix grads  = AllPixelGradients(U,SideLength);

	// using namespace std;
	// cout<<" > Inside Onestep_Direction():"<<endl;
	// cout<<"    * PixelGradients: "<<grads<<endl;
	// cout<<"    * beta: "<<beta<<endl;
	// cout<<"    * W: "<<W<<endl;
	// cout<<"    * Nu: "<<NU<<endl;
	// cout<<"    * NegAdjSum :"<<NegAdjSum<<endl;
	// cout<<"    * Term2 :"<<TermTwo<<endl;
	// cout<<"    * Term3 :"<<TermThree<<endl;


	Dk = -PixelGradientAdjointSum(beta*AllPixelGradients(U, SideLength) + beta*W + Nu, SideLength) + mu*prod(trans(A), prod(A, U) - B) - prod(trans(A), Lambda);

return Dk;
}

double U_Subfunction(BoostDoubleMatrix A, BoostDoubleVector U,
					 BoostDoubleVector B, BoostDoubleMatrix W, 
					 BoostDoubleMatrix Nu, BoostDoubleVector Lambda,
					 double beta, double mu, 
					 unsigned long SideLength)
{
	/*
	* Function: U_Subfunction
	* -----------------------
	* Calculate the quadratic cost function for the give problem state.
	*
	* Input --
	* A: an (M x N) projection matrix
	* U: an (N x 1) vector representing a rasterized image prediction
	* B: an (M x 1) set of observations
	* W: an (N x 2) set of dual variables corresponding to per-pixel (-voxel) gradients
	* Nu: an (N x 2) set of Lagrangian multipliers
	* Lambda: an (M x 1) set of Lagrangian multiplies
	* beta: scaling term on the matching between W and the true gradients
	* mu: scaling term on the matching between A*u and b
	* SideLength: the side length for the target image, i.e. N = SideLength^2
	*
	* Output -- a decimal value for the cost.
	*/
	double Q = 0.0;

	// Get all Gradients
	BoostDoubleMatrix Du = AllPixelGradients(U, SideLength);

	// Loop over pixels
	BoostDoubleVector Dui;
	BoostDoubleVector Wi;
	BoostDoubleVector Nui;
	for (unsigned long i = 0; i < U.size(); ++i) {
		Dui = GetRow(Du, i);
		Wi = GetRow(W, i);
		Nui = GetRow(Nu, i);

		BoostDoubleVector GradDiff = Dui - Wi;
		Q += -inner_prod(Nui, GradDiff) + (beta / 2) * SquareNorm(GradDiff);
	}

	// Residual Contribution
	BoostDoubleVector Residual = prod(A, U) - B;
	Q += -inner_prod(Lambda, Residual) + (mu / 2)*SquareNorm(Residual);

	return Q;
}

void Alternating_Minimisation(BoostDoubleMatrix A, BoostDoubleVector &U,
							  BoostDoubleVector B, BoostDoubleMatrix &W,
							  BoostDoubleMatrix Nu, BoostDoubleVector Lambda,
							  double beta, double mu,
							  unsigned long SideLength)
{
	/*
	* Function: Alternating_Minimisation
	* ----------------------------------
	* Calculate the minima U* and Wi* of the augmented Lagrangian function.
	*
	* Input --
	* A: an (M x N) projection matrix
	* U: an (N x 1) vector representing a rasterized image prediction
	* B: an (M x 1) set of observations
	* W: an (N x 2) set of dual variables corresponding to per-pixel (-voxel) gradients
	* Nu: an (N x 2) set of Lagrangian multipliers
	* Lambda: an (M x 1) set of Lagrangian multiplies
	* beta: scaling term on the matching between W and the true gradients
	* mu: scaling term on the matching between A*u and b
	* SideLength: the side length for the target image, i.e. N = SideLength^2
	*
	* Output -- None.
	*/

	using namespace std;

	unsigned long N = U.size();
	double delta = 0.00001;
	double rho = 0.6;
	double eta = 0.9995;
	double Pk = 1; 
	double C = Lagrangian(A, U, B, W, Nu, Lambda, beta, mu, SideLength, ANISOTROPIC);

	double armijo_tol, Qk, innerstop;
	double tol = 0.001;

	unsigned int LoopCounter = 0;
	unsigned int MaxIterations = 5;

	unsigned int ArmijoLoopCounter = 0;
	unsigned int MaxArmijoIterations = 5;

	BoostDoubleVector Uk_1 = BoostZeroVector(N);

	do
	{
		std::cout<<"  * AM Loop Iter [" << LoopCounter +1 << "]" << flush << endl;
		// cout << "    U(k) = [" << U << "]" << endl;
//*************************** "w sub-problem" ***************************
		W = ApplyShrike(W,Nu,beta,ANISOTROPIC);

//*************************** "u sub-problem" ***************************
		BoostDoubleVector Sk = U - Uk_1;
		BoostDoubleVector Dk = Onestep_Direction(A, U, B, W, Nu, Lambda, beta, mu, SideLength);
		BoostDoubleVector Yk = Dk - Onestep_Direction(A, Uk_1, B, W, Nu, Lambda, beta, mu, SideLength);

		//******** alpha = onestep_gradient ********

		std::cout << "   * U(k) - U(k-1) ="<<Sk<<std::endl;
		// std::cout << "   * D(k) ="<<Dk<<std::endl;
		std::cout << "   * D(k) - D(k-1) ="<<Yk<<std::endl;
		double numerator = inner_prod(Sk, Yk);
		double denominator = inner_prod(Yk, Yk);	
		double alpha = numerator/denominator;
		cout << " alpha = [" << numerator << " / " << denominator << " = " << alpha << "]" << endl;
		ArmijoLoopCounter = 0;
		do 
		{ 
			alpha = rho * alpha;
			BoostDoubleVector U_alphad = U - alpha*Dk;
			Qk = U_Subfunction(A, U_alphad, B, W, Nu, Lambda, beta, mu, SideLength);
			armijo_tol = C - delta*alpha*inner_prod(Dk, Dk);
			cout << "       > Armijo Iter [" << ArmijoLoopCounter + 1 << "]" << endl;
			cout << "         U - alpha *d: " << U_alphad << endl;
			//cout << " alpha: [" << alpha << "]" << flush << endl;
			//cout << " Q(k): [" << Qk << "]" << flush << endl;
			//cout << " armijo tol: [" << armijo_tol << "]" << flush << endl;
			ArmijoLoopCounter++;
		} while ((Qk > armijo_tol) && (ArmijoLoopCounter < MaxArmijoIterations));

		Uk_1 = U;
		U -= alpha * Dk;
		innerstop = norm_2(U - Uk_1);

//************************ Implement coefficents ************************
		double Pk1 = eta*Pk + 1;
		C = (eta*Pk*C + U_Subfunction(A, U, B, W, Nu, Lambda, beta, mu, SideLength))/Pk1;
		Pk = Pk1;

		std::cout<<"  * End AM Loop Iter ["<<LoopCounter<<"]  U Step Size : "<<alpha<<std::flush<<std::endl;
		LoopCounter++;
	} while ((innerstop > tol) && (LoopCounter < MaxIterations));
}

BoostDoubleMatrix tval3_reconstruction(BoostDoubleMatrix A, BoostDoubleVector y,
									   unsigned long SideLength)
{
	/*
	* Function: tval3_reconstruction
	* ------------------------------
	* Calculate the reconstructed image of the sample by the TVAL3 method.
	*
	* Input --
	* A: an (M x N) projection matrix
	* y: an (M x 1) set of observations
	* SideLength: the side length for the target image, i.e. N = SideLength^2
	*
	* Output -- an (L x L) reconstructed matrix.
	*/

	// unsigned long L = Sinogram.size1(); // Size of the sample (in pixels)
	// unsigned long O = Sinogram.size2(); // Numbers of tilt angles
	// unsigned long M = L * O; // Numbers of measurements
	// unsigned long N = L * L; // Size of rasterized Image vector
	unsigned long M = A.size1();
	unsigned long N = A.size2();
	unsigned long L = SideLength; // allowing truncation
	
	double mu = 1024.0;
	double beta = 1024.0;
	double coef = 1.0;
	double outerstop;
	double tol = 0.001;
	unsigned int LoopCounter = 0;
	unsigned int MaxIterations = 5;
	
	// BoostDoubleVector U = BoostZeroVector(N); // U(0) = 0 for all i
	BoostDoubleVector U = prod(trans(A),y);
	BoostDoubleVector Uk_1 = BoostZeroVector(N); 
	BoostDoubleVector Lambda = BoostZeroVector(M);
	// BoostDoubleVector B = MatrixToVector(Sinogram);

	BoostDoubleMatrix Du = AllPixelGradients(U, L);
	BoostDoubleMatrix Nu = BoostZeroMatrix(N, 2);
	BoostDoubleMatrix W = ApplyShrike(Du, Nu, beta, ANISOTROPIC);	
	// BoostDoubleMatrix A = CreateRandomMatrix(M, N);

	using namespace std;
	
	do
	{
		cout << "Outer Iter [" << LoopCounter +1 << "]" << endl;
		// cout << " U* = [" << U << "]" << endl;
		Uk_1 = U;
		Alternating_Minimisation(A, U, y, W, Nu, Lambda, beta, mu, L);
		BoostDoubleMatrix Du = AllPixelGradients(U, L);
		Nu = Nu - beta*(Du - W);
		Lambda = Lambda - mu*(prod(A, U) - y);

		beta = coef*beta;
		mu = coef*mu;

		// cout<<"Outer Iter ["<<LoopCounter<<"] U: "<<U<<endl;
		cout<<"End Outer Iter ["<<LoopCounter<<"] mu : "<<mu<<" | beta : "<<beta<<endl;
		outerstop = norm_2(U - Uk_1);
		LoopCounter++;
	} while (outerstop > tol && LoopCounter < MaxIterations);

return VectorToMatrix(U, L, L);
}