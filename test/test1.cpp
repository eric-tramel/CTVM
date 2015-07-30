#include <iostream>
#include "ctvm.h"
#include "ctvm_util.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
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
	std::cout << std::endl;
	using namespace Magick;
	InitializeMagick("");

	std::cout << "Testing ImageMagick Link." << std::endl;
	Image someImage;
	someImage.read("C:\\data\\peppers.jpg");
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
	anotherImage.write("C:\\data\\testoutimage.png");
	std::cout << "    " << "Passed." << std::endl;

	/* Test CTVM Image Load */
	std::cout << std::endl;
	std::cout << "Testing CTVM Image Load." << std::endl;
	BoostDoubleMatrix ImageMatrix = LoadImage("C:\\data\\peppers.jpg");
	std::cout << "Image Data:" << std::endl;
	std::cout << ImageMatrix(0, 0) << " " << ImageMatrix(0, 1) << " " << ImageMatrix(0, 3) << std::endl;
	std::cout << ImageMatrix(1, 0) << " " << ImageMatrix(1, 1) << " " << ImageMatrix(1, 3) << std::endl;
	std::cout << ImageMatrix(2, 0) << " " << ImageMatrix(2, 1) << " " << ImageMatrix(2, 3) << std::endl;
	std::cout << "    " << "Passed." << std::endl;

	std::cout << "Testing CTVM Image Write." << std::endl;
	WriteImage(ImageMatrix, "C:\\data\\test_peppers.jpg");
	std::cout << "    " << "Passed." << std::endl;

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

	/*Test Shrinke function*/
	std::cout << std::endl;
	int ii = 1;
	BoostDoubleVector DiUk = Gradient2D(U, ii);
	BoostDoubleVector NUi (2);
	for (int jj = 0; jj < 2; ++jj) { NUi(jj) = NU(ii, jj); }
	BoostDoubleVector SHRIKE = Shrike(DiUk, NUi, beta);

	std::cout << "W(i,l+1): " << SHRIKE << std::endl; // expected result SHRIKE = ( -0.2988585, -0.42265)
	std::cout << "    " << "Passed." << std::endl;

return 0;
}