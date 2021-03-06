#include <iostream>
#include "ctvm.h"
#include "ctvm_util.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Magick++.h>
#include <ctime>

#define prefix "   * "

double GetSeconds(clock_t elapsed) {
	return ((float)elapsed / CLOCKS_PER_SEC);
}

std::string ReportTime(clock_t elapsed) {
	return "<" + std::string(std::to_string(GetSeconds(elapsed))) + " sec.>";
}

void TestRandomMatrix() {
	using namespace std;
	cout << "Random Matrix Test" << endl;
	cout << "------------------" << endl;
	cout << prefix << "Generating tiny random matrix (3x3) matrix..." << flush;
	clock_t t = clock();
	BoostDoubleMatrix RandomMatrix = CreateRandomMatrix(3, 3);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Generating large random matrix (1000x1000)..." << flush;
	t = clock();
	RandomMatrix = CreateRandomMatrix(1000, 1000);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestNormalization() {
	using namespace std;

	cout << "Normalization Test" << endl;
	cout << "------------------" << endl;

	BoostDoubleMatrix TestMatrix(3, 3);
	TestMatrix(0, 0) = 1; TestMatrix(0, 1) = 5; TestMatrix(0, 2) = 3;
	TestMatrix(1, 0) = 2; TestMatrix(1, 1) = 11; TestMatrix(1, 2) = 5;
	TestMatrix(2, 0) = 1; TestMatrix(2, 1) = 10; TestMatrix(2, 2) = 6;
	cout << prefix << "Original Matrix: " << TestMatrix << endl;


	clock_t t = clock();
	BoostDoubleMatrix NormMatrix = NormalizeMatrix(TestMatrix);
	t = clock() - t;
	cout << prefix << "Normalized Matrix: " << NormMatrix << " " << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestImageMagick(char* test_image) {
	using namespace std;
	using namespace Magick;
	clock_t t;

	cout << "ImageMagick Test" << endl;
	cout << "----------------" << endl;

	cout << prefix << "Initializing ImageMagic..." << flush;
	t = clock();
	InitializeMagick("");
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;


	cout << prefix << "Reading image (" << test_image << ")..." << flush;
	t = clock();
	Image someImage;
	someImage.read(test_image);
	t = clock() - t;
	//cout<<"done. "<<ReportTime(t)<<endl;
	cout << "done. [" << someImage.rows() << "x" << someImage.columns() << "]. " << ReportTime(t) << endl;


	cout << prefix << "Passed." << endl << endl;
}

void TestMatrixIO(char* test_image, char* output_image) {
	using namespace std;
	clock_t t;

	cout << "Matrix I/O Test" << endl;
	cout << "---------------" << endl;

	/* Tiny Resolution */
	cout << prefix << "Loading image file into matrix at [32x32]..." << flush;
	t = clock();
	BoostDoubleMatrix TinyImageMatrix = LoadImage(test_image, 32, 32);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Saving to file from image matrix at [32x32]..." << flush;
	t = clock();
	WriteImage(TinyImageMatrix, output_image);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	/* Medium Resolution */
	cout << prefix << "Loading image file into matrix at [512,512]..." << flush;
	t = clock();
	BoostDoubleMatrix MediumImageMatrix = LoadImage(test_image, 512, 512);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Saving to file from image matrix at [512x512]..." << flush;
	t = clock();
	WriteImage(MediumImageMatrix, output_image);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	/* High Resolution */
	cout << prefix << "Loading image file into matrix at [2048,2048]..." << flush;
	t = clock();
	BoostDoubleMatrix LargeImageMatrix = LoadImage(test_image, 2048, 2048);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Saving to file from image matrix at [2048x2048]..." << flush;
	t = clock();
	WriteImage(LargeImageMatrix, output_image);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestNeighborCheck() {
	using namespace std;
	cout << "Neighbor Index Test" << endl;
	cout << "-------------------" << endl;
	cout << prefix << "Testing for a [2x2] Matrix" << endl;

	cout << prefix << "[0] Neighbors: H=" << RightNeighbor(0, 2) << ", V=" << DownNeighbor(0, 2) << endl;
	cout << prefix << "[1] Neighbors: H=" << RightNeighbor(1, 2) << ", V=" << DownNeighbor(1, 2) << endl;
	cout << prefix << "[2] Neighbors: H=" << RightNeighbor(2, 2) << ", V=" << DownNeighbor(2, 2) << endl;
	cout << prefix << "[3] Neighbors: H=" << RightNeighbor(3, 2) << ", V=" << DownNeighbor(3, 2) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestGradient() {
	using namespace std;
	cout << "Gradient Test" << endl;
	cout << "-------------" << endl;
	cout << prefix << "Testing for a [2x2] Matrix" << endl;

	/* Gradient Test */
	BoostDoubleMatrix GradTest(2, 2);
	GradTest(0, 0) = 1;
	GradTest(1, 0) = 2;
	GradTest(0, 1) = 3;
	GradTest(1, 1) = 4;
	std::cout << "Gradient Test on A = [1 3; 2 4]" << std::endl;
	std::cout << "  " << AllPixelGradients(MatrixToVector(GradTest), 2) << std::endl;

	/* Gradient Adjoint Test */
	// We should have a set of pixel gradients from A = [1 3; 2 4]
	// which look like
	// G = [-2 -1; -2 0; 0 -1; 0 0]
	// then the correct answer is a sum vector of
	// [-3, -1, 1, 3]
	BoostDoubleMatrix thisGradient = AllPixelGradients(MatrixToVector(GradTest), 2);
	std::cout << "Adjoint Sum Test." << std::endl;
	std::cout << "  Sum Vector = " << PixelGradientAdjointSum(thisGradient, 2) << std::endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestGradient(char* InputFile, char* OutputFile) {
	using namespace std;
	clock_t t;
	unsigned int L = 128;

	cout << "Gradient Test" << endl;
	cout << "-------------" << endl;

	/* Medium Resolution */
	cout << prefix << "Loading image file into matrix at [" << L << "," << L << "]..." << flush;
	t = clock();
	BoostDoubleMatrix MediumImageMatrix = LoadImage(InputFile, L, L);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	/* To Vector */
	BoostDoubleVector ImageVect = MatrixToVector(MediumImageMatrix);

	/* Compute Gradients */
	cout << prefix << "Computing gradients..." << flush;
	t = clock();
	BoostDoubleMatrix Gradients = AllPixelGradients(ImageVect, L);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	/* Get out just the horizontal */
	BoostDoubleVector HorzGradients = GetCol(Gradients, 0);

	/* Write as image */
	WriteImage(NormalizeMatrix(VectorToMatrix(HorzGradients, L, L)), OutputFile);
}

void TestShrike() {
	using namespace std;
	// With the above settings we expect that the results should be:
	// Isotropic @ beta = 0.5  : (0,0)
	// Isotropic @ beta = 0.65 : (0.0557...,0.0325...)
	// Anisotropic @ beta = 0.5  : (0,0)
	// Anisotropic @ beta = 10   : (0.925...,0.8875...)
	BoostDoubleVector g(2);  g(0) = 1;      g(1) = 1;
	BoostDoubleVector nu(2); nu(0) = -0.25; nu(1) = 0.125;
	BoostDoubleVector Result;
	clock_t t;

	cout << "Shrinkage-Like Function Test" << endl;
	cout << "----------------------------" << endl;
	cout << prefix << "Set g = " << g << endl;
	cout << prefix << "Set nu = " << nu << endl;

	cout << prefix << "Testing ShrikeAnisotropic @ beta = 0.5..." << flush;
	t = clock();
	Result = ShrikeAnisotropic(g, nu, 0.5);
	t = clock() - t;
	cout << "done. [" << Result << "]. " << ReportTime(t) << endl;

	cout << prefix << "Testing ShrikeAnisotropic @ beta = 10..." << flush;
	t = clock();
	Result = ShrikeAnisotropic(g, nu, 10);
	t = clock() - t;
	cout << "done. [" << Result << "]. " << ReportTime(t) << endl;

	cout << prefix << "Testing ShrikeIsotropic @ beta = 0.5..." << flush;
	t = clock();
	Result = ShrikeIsotropic(g, nu, 0.5);
	t = clock() - t;
	cout << "done. [" << Result << "]. " << ReportTime(t) << endl;

	cout << prefix << "Testing ShrikeIsotropic @ beta = 0.65..." << flush;
	t = clock();
	Result = ShrikeIsotropic(g, nu, 0.65);
	t = clock() - t;
	cout << "done. [" << Result << "]. " << ReportTime(t) << endl;

	BoostDoubleMatrix G(3, 2);
	SetRow(G, g, 0);
	SetRow(G, g, 1);
	SetRow(G, g, 2);
	BoostDoubleMatrix Nu(3, 2);
	SetRow(Nu, nu, 0);
	SetRow(Nu, nu, 1);
	SetRow(Nu, nu, 2);
	cout << prefix << "Set G = " << G << endl;
	cout << prefix << "Set Nu = " << Nu << endl;

	cout << prefix << "Testing ApplyShrike(Anisotropic) @ beta = 10..." << flush;
	t = clock();
	BoostDoubleMatrix Shriked = ApplyShrike(G, Nu, 10, ANISOTROPIC);
	t = clock() - t;
	cout << "done. [" << Shriked << "]. " << ReportTime(t) << endl;

	cout << prefix << "Testing ApplyShrike(Isotropic) @ beta = 0.65..." << flush;
	t = clock();
	Shriked = ApplyShrike(G, Nu, 0.65, ISOTROPIC);
	t = clock() - t;
	cout << "done. [" << Shriked << "]. " << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestRasterization() {
	using namespace std;
	/* Rasterization Test */
	BoostDoubleMatrix RasterTest(2, 2);
	RasterTest(0, 0) = 1;
	RasterTest(1, 0) = 2;
	RasterTest(0, 1) = 3;
	RasterTest(1, 1) = 4;
	clock_t t;

	cout << "Matrix Rasterization Test" << endl;
	cout << "-------------------------" << endl;

	cout << prefix << "Set Matrix as [" << RasterTest << "]" << endl;
	cout << prefix << "Rasterized Matrix: " << MatrixToVector(RasterTest) << endl;
	cout << prefix << "Restored Matrix: " << VectorToMatrix(MatrixToVector(RasterTest), 2, 2) << endl;

	BoostDoubleMatrix LargeMatrix(2048, 2048);
	cout << prefix << "Rasterizing large matrix [2048x2048]..." << flush;
	t = clock();
	BoostDoubleVector LargeVector = MatrixToVector(LargeMatrix);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Restoring large matrix [2048x2048]..." << flush;
	t = clock();
	RasterTest = VectorToMatrix(LargeVector, 2048, 2048);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestLagrangian() {
	using namespace std;

	BoostDoubleMatrix A(3, 4), W(4, 2), Nu(4, 2);
	BoostDoubleVector U(4), B(3), Lambda(3);
	double beta = 0.25;
	double mu = 0.5;
	clock_t t;

	/* Init Projections */
	A(0, 0) = 1; A(0, 1) = 1; A(0, 2) = 1; A(0, 3) = 1;
	A(1, 0) = 2; A(1, 1) = 2; A(1, 2) = 2; A(1, 3) = 2;
	A(2, 0) = 3; A(2, 1) = 3; A(2, 2) = 3; A(2, 3) = 3;

	/* Init Gradient Dual Variables  */
	W(0, 0) = -1; W(0, 1) = -1;
	W(1, 0) = -1; W(1, 1) = -1;
	W(2, 0) = -1; W(2, 1) = -1;
	W(3, 0) = -1; W(3, 1) = -1;

	/* Init Dual Lag. Multipliers */
	Nu(0, 0) = 0.25; Nu(0, 1) = 0.25;
	Nu(1, 0) = 0.25; Nu(1, 1) = 0.25;
	Nu(2, 0) = 0.25; Nu(2, 1) = 0.25;
	Nu(3, 0) = 0.25; Nu(3, 1) = 0.25;

	/* Init Image Vector */
	U(0) = 1;
	U(1) = 2;
	U(2) = 3;
	U(3) = 4;

	/* Init Observations */
	B(0) = 3;
	B(1) = 3;
	B(2) = 3;

	/* Init Observation Lag. Multipliers */
	Lambda(0) = 0.5;
	Lambda(1) = 0.5;
	Lambda(2) = 0.5;

	/* Test Lagrangian function */
	cout << "Lagrangian Test" << endl;
	cout << "---------------" << endl;

	cout << prefix << "Set A = " << A << endl;
	cout << prefix << "Set U = " << U << endl;
	cout << prefix << "Set B = " << B << endl;
	cout << prefix << "Set W = " << W << endl;
	cout << prefix << "Set Nu = " << Nu << endl;
	cout << prefix << "Set Lambda = " << Lambda << endl;
	cout << prefix << "Set beta = " << beta << endl;
	cout << prefix << "Set mu = " << mu << endl;

	cout << prefix << "Calculating Lagrangian..." << flush;
	t = clock();
	double L = Lagrangian(A, U, B, W, Nu, Lambda, beta, mu, 2, ISOTROPIC);
	t = clock() - t;
	cout << "done. [" << L << "]. [247.1569...] Expected. " << ReportTime(t) << endl;

	/* Calulating Langrangian on large dataset */
	A = BoostDoubleMatrix(64 * 64, 128 * 128);
	W = BoostDoubleMatrix(128 * 128, 2);
	Nu = BoostDoubleMatrix(128 * 128, 2);
	B = BoostDoubleVector(64 * 64);
	Lambda = BoostDoubleVector(64 * 64);
	U = BoostDoubleVector(128 * 128);

	cout << prefix << "Calculating Lagrangian(ISO) on large dataset..." << flush;
	t = clock();
	L = Lagrangian(A, U, B, W, Nu, Lambda, beta, mu, 128, ISOTROPIC);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Calculating Lagrangian(ANISO) on large dataset..." << flush;
	t = clock();
	L = Lagrangian(A, U, B, W, Nu, Lambda, beta, mu, 128, ANISOTROPIC);
	t = clock() - t;
	cout << "done. " << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestOnestep_Direction() {
	using namespace std;

	BoostDoubleMatrix A(3, 4), W(4, 2), Nu(4, 2);
	BoostDoubleVector U(4), B(3), Lambda(3);
	double beta = 0.25;
	double mu = 0.5;
	clock_t t;

	/* Init Projections */
	A(0, 0) = 1; A(0, 1) = 1; A(0, 2) = 1; A(0, 3) = 1;
	A(1, 0) = 2; A(1, 1) = 2; A(1, 2) = 2; A(1, 3) = 2;
	A(2, 0) = 3; A(2, 1) = 3; A(2, 2) = 3; A(2, 3) = 3;

	/* Init Gradient Dual Variables  */
	W(0, 0) = -2; W(0, 1) = -1;
	W(1, 0) = -2; W(1, 1) = 0;
	W(2, 0) = 0; W(2, 1) = -1;
	W(3, 0) = 0; W(3, 1) = 1;

	/* Init Dual Lag. Multipliers */
	Nu(0, 0) = 0.25; Nu(0, 1) = 0.25;
	Nu(1, 0) = 0.25; Nu(1, 1) = 0.25;
	Nu(2, 0) = 0.25; Nu(2, 1) = 0.25;
	Nu(3, 0) = 0.25; Nu(3, 1) = 0.25;

	/* Init Image Vector */
	U(0) = 1;
	U(1) = 2;
	U(2) = 3;
	U(3) = 4;

	/* Init Observations */
	B(0) = 3;
	B(1) = 3;
	B(2) = 3;

	/* Init Observation Lag. Multipliers */
	Lambda(0) = 0.5;
	Lambda(1) = 0.5;
	Lambda(2) = 0.5;

	/*Test Onestep direction*/
	cout << "One-step direction Test" << endl;
	cout << "-----------------------" << endl;

	cout << prefix << "Expected result: [58.75  58.25  57.75  57.25]" << endl;

	cout << prefix << "Calculating direction value..." << flush;
	t = clock();
	BoostDoubleVector D = Onestep_Direction(A, U, B, W, Nu, Lambda, beta, mu, 2);
	t = clock() - t;
	cout << "done. [" << D << "]." << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

void TestU_Subfunction() {
	using namespace std;

	BoostDoubleMatrix A(3, 9), W(9, 2), Nu(9, 2);
	BoostDoubleVector U(9), B(3), Lambda(3);
	double beta = 0.3333;
	double mu = sqrt(2);
	clock_t t;

	/* Init Projections */
	A(0, 0) = 1; A(0, 1) = 1; A(0, 2) = 1; A(0, 3) = 1; A(0, 4) = 1; A(0, 5) = 1; A(0, 6) = 1; A(0, 7) = 1; A(0, 8) = 1;
	A(1, 0) = 2; A(1, 1) = 2; A(1, 2) = 2; A(1, 3) = 2; A(1, 4) = 2; A(1, 5) = 2; A(1, 6) = 2; A(1, 7) = 2; A(1, 8) = 2;
	A(2, 0) = 3; A(2, 1) = 3; A(2, 2) = 3; A(2, 3) = 3;	A(2, 4) = 3; A(2, 5) = 3; A(2, 6) = 3; A(2, 7) = 3; A(2, 8) = 3;

	/* Init Gradient Dual Variables  */
	W(0, 0) = -1; W(0, 1) = 1;
	W(1, 0) = 1; W(1, 1) = -1;
	W(2, 0) = 1; W(2, 1) = 1;
	W(3, 0) = -1; W(3, 1) = -1;
	W(4, 0) = 0; W(4, 1) = 0;
	W(5, 0) = -1; W(5, 1) = 1;
	W(6, 0) = 1; W(6, 1) = -1;
	W(7, 0) = 1; W(7, 1) = 1;
	W(8, 0) = -1; W(8, 1) = -1;

	/* Init Dual Lag. Multipliers */
	Nu(0, 0) = 0.1; Nu(0, 1) = 0.2;
	Nu(1, 0) = 0.3; Nu(1, 1) = 0.4;
	Nu(2, 0) = 0.5; Nu(2, 1) = 0.6;
	Nu(3, 0) = 0.7; Nu(3, 1) = 0.8;
	Nu(4, 0) = 0.9; Nu(4, 1) = 1.0;
	Nu(5, 0) = 1.1; Nu(5, 1) = 1.2;
	Nu(6, 0) = 1.3; Nu(6, 1) = 1.4;
	Nu(7, 0) = 1.5; Nu(7, 1) = 1.6;
	Nu(8, 0) = 1.7; Nu(8, 1) = 1.8;

	/* Init Image Vector */
	U(0) = -2;
	U(1) = 0;
	U(2) = 1;
	U(3) = -2.5;
	U(4) = 3;
	U(5) = 4.3;
	U(6) = 1;
	U(7) = 0;
	U(8) = -1;

	/* Init Observations */
	B(0) = 1;
	B(1) = 0;
	B(2) = 3;

	/* Init Observation Lag. Multipliers */
	Lambda(0) = 0.25;
	Lambda(1) = 0.1;
	Lambda(2) = 0.5;

	/*Test quadratic function*/
	cout << "Quadratic function Test" << endl;
	cout << "-----------------------" << endl;
	cout << prefix << "Expected result: [111.841]" << endl;

	cout << prefix << "Calculating quadratic value for U=[9]..." << flush;
	t = clock();
	double Q = U_Subfunction(A, U, B, W, Nu, Lambda, beta, mu, 3);
	t = clock() - t;
	cout << "done. [" << Q << "]." << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;

}

void TestReconstruction(int argc, char **argv) {
	using namespace std;
	clock_t t;
	if (argc != 4) {
		cout << "Usage: test <imsize> <orig_image> <output_image>" << endl;
	}

	/* Read Input */
	char* ImageSizeStr = argv[1];
	char* OriginalImageFile = argv[2];
	char* OutputImageFile = argv[3];
	unsigned int L = atoi(ImageSizeStr);

	/* Specify Problem Settings */
	double MeasurementRate = 1.0;
	double NoiseVariance = 0.0;
	unsigned int N = L*L;
	unsigned int M = MeasurementRate*N;     // Allow truncation

	cout << "Running CS Experiment for L = " << L << " \\Alpha = " << MeasurementRate << ", \\Delta = " << NoiseVariance << endl;

	/* Load Image */
	cout << " * Loading image (" << OriginalImageFile << ")..." << flush;
	BoostDoubleMatrix XImage = LoadImage(OriginalImageFile, L, L);
	BoostDoubleVector XVect = MatrixToVector(XImage);
	cout << "done." << endl;

	/* Create Projection Matrix */
	cout << " * Creating Random Matrix (" << M << "x" << N << ")..." << flush;
	BoostDoubleMatrix A = CreateRandomMatrix(M, N);
	cout << "done." << endl;

	/* Create Measurements */
	cout << " * Generating Measurements..." << flush;
	BoostDoubleVector y = prod(A, XVect) + sqrt(NoiseVariance)*CreateRandomVector(M);
	cout << "done." << endl;

	/* Perform Reconstruction*/
	// Testing the tval3_reconstruction method --> need to call sinogram and tilt angles file...
	cout << " * Computing the recorded image by the TVAL3 method..." << endl;
	t = clock();
	BoostDoubleMatrix XRecImage = tval3_reconstruction(A, y, L);
	t = clock() - t;
	cout << "done." << ReportTime(t) << endl;

	/* Reshape Output */
	// XRecImage = VectorToMatrix(XRecImage,L,L);

	/* Write Result */
	cout << " * Writing result to image (" << OutputImageFile << ")..." << flush;
	WriteImage(NormalizeMatrix(XRecImage), OutputImageFile);
	cout << "done." << endl;
}

int main(int argc, char **argv) {
	using namespace std;
	cout << endl;

	if (argc < 3) {
		TestRasterization();
		TestRandomMatrix();
		TestNormalization();
		TestNeighborCheck();
		TestGradient();
		TestShrike();
		TestLagrangian();
		TestOnestep_Direction();
		TestU_Subfunction();
	}

	if (argc == 3) {
		// Only run tests requiring File I/O if the file names have been passed.
		TestImageMagick(argv[1]);
		TestMatrixIO(argv[1], argv[2]);
		TestGradient(argv[1], argv[2]);
	}

	if (argc == 4) {
		// Only run tests requiring Size and File I/O if the size length and the file names have been passed.
		TestReconstruction(argc, argv);
	}
	return 0;
}