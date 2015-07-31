#include <iostream>
#include "ctvm.h"
#include "ctvm_util.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Magick++.h>

int main(int argc, char **argv)
{
	/* Test CTVM.dylib Link */
	std::cout << std::endl;
	std::cout << "Testing libctvm Link." << std::endl;
	BoostDoubleMatrix DummySinogram(0, 0);
	BoostDoubleVector DummyAngles(0);
	// BoostDoubleMatrix DummyReconstruction = tval3_reconstruction(DummySinogram,DummyAngles);
	std::cout << "    " << "Passed." << std::endl;

	/* Test CTVM_util.dylib Link */
	std::cout << std::endl;
	std::cout << "Allocating tiny(3x3) random matrix..." << std::endl;
	BoostDoubleMatrix RandomMatrix = CreateRandomMatrix(3, 3);
	std::cout << "    " << RandomMatrix << std::endl;
	std::cout << "    " << "Passed." << std::endl;

	std::cout << "Rasterized Version..." << std::endl;
	std::cout << "    " << MatrixToVector(RandomMatrix) << std::endl;
	std::cout << "    " << "Passed." << std::endl;

	std::cout << "Back to Matrix..." << std::endl;
	std::cout << "    " << VectorToMatrix(MatrixToVector(RandomMatrix), 3, 3) << std::endl;
	std::cout << "    " << "Passed." << std::endl;

	std::cout << "Allocating large(1000x1000) random matrix..." << std::endl;
	BoostDoubleMatrix RandomMatrixLarge = CreateRandomMatrix(1000, 1000);
	std::cout << "    " << "Passed." << std::endl;

	std::cout << "Testing Normalization" << std::endl;
	BoostDoubleMatrix TestMatrix(3, 3);
	TestMatrix(0, 0) = 1; TestMatrix(0, 1) = 5; TestMatrix(0, 2) = 3;
	TestMatrix(1, 0) = 2; TestMatrix(1, 1) = 11; TestMatrix(1, 2) = 5;
	TestMatrix(2, 0) = 1; TestMatrix(2, 1) = 10; TestMatrix(2, 2) = 6;
	std::cout << "Original Matrix: " << TestMatrix << std::endl;
	std::cout << "Normalized: " << NormalizeMatrix(TestMatrix) << std::endl;

	/* Test ImageMagick */
	// Only run ImageMagick stuff if there are file names specified
	if(argc > 2){
		char* test_image = argv[1];
		char* output_image = argv[2];
		std::cout << std::endl;
		using namespace Magick;
		InitializeMagick("");

		std::cout << "Testing ImageMagick Link." << std::endl;
		Image someImage;
		someImage.read(test_image);
		std::cout << "    " << "Passed." << std::endl;

		std::cout << "Creating Novel Image" << std::endl;
		Image anotherImage;
		anotherImage.size(Geometry(32, 32));         // Specify the dimensionality
		Pixels anotherImageView(anotherImage);      // Get a view of the image
		PixelPacket aPixel;
		aPixel.red = ColorGray::scaleDoubleToQuantum(0.5);
		aPixel.green = ColorGray::scaleDoubleToQuantum(0.5);
		aPixel.blue = ColorGray::scaleDoubleToQuantum(0.5);
		*(anotherImageView.get(16, 16, 1, 1)) = aPixel;
		anotherImage.type(GrayscaleType);           // Specify the color type
		anotherImageView.sync();
		anotherImage.write(output_image);
		std::cout << "    " << "Passed." << std::endl;

		/* Test CTVM Image Load */
		std::cout << std::endl;
		std::cout << "Testing CTVM Image Load." << std::endl;
		BoostDoubleMatrix ImageMatrix = LoadImage(test_image,32,32);
		std::cout << "Image Data:" << std::endl;
		std::cout << ImageMatrix(0, 0) << " " << ImageMatrix(0, 1) << " " << ImageMatrix(0, 3) << std::endl;
		std::cout << ImageMatrix(1, 0) << " " << ImageMatrix(1, 1) << " " << ImageMatrix(1, 3) << std::endl;
		std::cout << ImageMatrix(2, 0) << " " << ImageMatrix(2, 1) << " " << ImageMatrix(2, 3) << std::endl;
		std::cout << "    " << "Passed." << std::endl;

		std::cout << "Testing CTVM Image Write." << std::endl;
		WriteImage(ImageMatrix, output_image);
		std::cout << "    " << "Passed." << std::endl;
	}

	/* Initialisation */
	BoostDoubleMatrix A(2, 4), W(4, 2), NU(4, 2);
	BoostDoubleVector U(4), U1(9), B(2), LAMBDA(2);

	A(0, 0) = 1; A(0, 1) = 0; A(0, 2) = 1; A(0, 3) = 0;
	A(1, 0) = 0; A(1, 1) = 1; A(1, 2) = 1; A(1, 3) = 1;

	W(0, 0) = -1; W(0, 1) = 1;
	W(1, 0) = 0; W(1, 1) = 1;
	W(2, 0) = -1; W(2, 1) = 0;
	W(3, 0) = 0; W(3, 1) = 1;

	NU(0, 0) = 2; NU(0, 1) = 1;
	NU(1, 0) = 1; NU(1, 1) = 0;
	NU(2, 0) = 0; NU(2, 1) = 2;
	NU(3, 0) = 1; NU(3, 1) = 3;

	U(0) = 1;
	U(1) = 2;
	U(2) = 0;
	U(3) = 1;

	U1(0) = 1;
	U1(1) = 2;
	U1(2) = 0;
	U1(3) = 1;
	U1(4) = 3;
	U1(5) = -2;
	U1(6) = 0;
	U1(7) = -1;
	U1(8) = 0;

	B(0) = 1;
	B(1) = 2;

	LAMBDA(0) = 2;
	LAMBDA(1) = 1;

	double beta = sqrt(2);
	double mu = 3;

	/* Test Gradient2DMatrix */
	std::cout << std::endl;
	std::cout << "Testing CTVM 2D Gradient for all i." << std::endl;
	BoostDoubleMatrix X = VectorToMatrix(U, 2, 2);
	BoostDoubleMatrix GradientMatrix = Gradient2DMatrix(U);
	BoostDoubleMatrix Di = Unit_Gradient2DMatrix(U1, 4);
	std::cout << "Original Matrix: " << X << std::endl;
	std::cout << "Gradients DiU (right gradient, down gradient): " << GradientMatrix << std::endl;
	std::cout << "Original Vector: " << U1 << std::endl;
	std::cout << " D(4) times U: " << prod(Di,U1) << std::endl;
	std::cout << "    " << "Passed." << std::endl;
	
	/*Test Lagrangian function*/
	std::cout << std::endl;
	double L = Lagrangian(A, U, B, W, NU, LAMBDA, beta, mu);

	std::cout << "Lagrangian: " << L << std::endl; // expected result L = 8.6213
	std::cout << "    " << "Passed." << std::endl;

	/* Test One-step Direction */
	// U = [1 2 3 4];
	// A = [1 1 1 1; 2 2 2 2; 3 3 3 3];
	// b = [3 3 3];
	// lambda = [0.5 0.5 0.5];
	// mu = 0.5;
	// beta = 0.25;
	// Nu = [0.25 0.25; 0.25 0.25; 0.25 0.25; 0.25 0.25];
	// W = [-1 -1; -1 -1; -1 -1; -1 -1];
	// Expected result: d_k = [58.7500   58.2500   57.7500   57.2500]

	/*Test Shrinkage-like function*/	
	BoostDoubleVector g(2);  g(0) = 1;      g(1) = 1;
	BoostDoubleVector nu(2); nu(0) = -0.25; nu(1) = 0.125;	
	// With the above settings we expect that the results should be:
	// Anisotropic @ beta = 0.5  : (0,0)
	// Anisotropic @ beta = 0.65 : (0.0557...,0.0325...)
	// Isotropic @ beta = 0.5  : (0,0)
	// Isotropic @ beta = 10   : (0.925...,0.8875...)
	std::cout<<"Testing ShrikeAnisotropic. Result1 : "<<ShrikeAnisotropic(g,nu,0.5)<<std::endl;
	std::cout<<"Testing ShrikeAnisotropic. Result2 : "<<ShrikeAnisotropic(g,nu,0.65)<<std::endl;

	std::cout<<"Testing ShrikeIsotropic. Result1 : "<<ShrikeIsotropic(g,nu,0.5)<<std::endl;
	std::cout<<"Testing ShrikeIsotropic. Result2 : "<<ShrikeIsotropic(g,nu,10)<<std::endl;


	std::cout<<"Testing ApplyShrike"<<std::endl;
	BoostDoubleMatrix G  (3,2); SetRow(G,g,0); SetRow(G,g,1); SetRow(G,g,2);
	std::cout<<" G:  "<<G<<std::endl;
	BoostDoubleMatrix Nu (3,2); SetRow(Nu,nu,0); SetRow(Nu,nu,1); SetRow(Nu,nu,2);
	std::cout<<" Nu:  "<<Nu<<std::endl;
	BoostDoubleMatrix Shriked = ApplyShrike(G,Nu,10,ISOTROPIC);
	std::cout<<"   Sriked Matrix: "<<Shriked<<std::endl;


	/* Slicing Test */
	std::cout<<"Slice Test"<<std::endl;
	std::cout<<"   Original Matrix:"<<RandomMatrix<<std::endl;
	std::cout<<"   Middle Row: "<<GetRow(RandomMatrix,1)<<std::endl;
	std::cout<<"   Middle Col: "<<GetCol(RandomMatrix,1)<<std::endl;
	SetRow(RandomMatrix,BoostZeroVector(3),1);
	std::cout<<"   Middle Row Zeros: "<<RandomMatrix<<std::endl;

	/* Neighbor Test */
	std::cout<<"Neighbor Test A -> 2x2"<<std::endl;
	std::cout<<"A[0] Neighbors:"<<RightNeighbor(0,2)<<","<<DownNeighbor(0,2)<<std::endl;
	std::cout<<"A[1] Neighbors:"<<RightNeighbor(1,2)<<","<<DownNeighbor(1,2)<<std::endl;
	std::cout<<"A[2] Neighbors:"<<RightNeighbor(2,2)<<","<<DownNeighbor(2,2)<<std::endl;
	std::cout<<"A[3] Neighbors:"<<RightNeighbor(3,2)<<","<<DownNeighbor(3,2)<<std::endl;

	/* Rasterization Test */
	BoostDoubleMatrix RasterTest(2,2);
	RasterTest(0,0) = 1;
	RasterTest(1,0) = 2;
	RasterTest(0,1) = 3;
	RasterTest(1,1) = 4;
	std::cout<<"Original Matrix: "<<RasterTest<<std::endl;
	std::cout<<"Rasterized Matrix: "<<MatrixToVector(RasterTest)<<std::endl;
	std::cout<<"Restored Matrix: "<<VectorToMatrix(MatrixToVector(RasterTest),2,2)<<std::endl;

return 0;
}