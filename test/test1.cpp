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

double GetSeconds(clock_t elapsed){
return ((float)elapsed /CLOCKS_PER_SEC);
}

std::string ReportTime(clock_t elapsed){
return "<" + std::string(std::to_string(GetSeconds(elapsed))) + " sec.>";
}

void TestRandomMatrix(){
	using namespace std;
	cout<<"Random Matrix Test"<<endl;
	cout<<"------------------"<<endl;
	cout<<prefix<<"Generating tiny random matrix (3x3) matrix..."<<flush;
	clock_t t = clock();
	BoostDoubleMatrix RandomMatrix = CreateRandomMatrix(3,3);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Generating large random matrix (1000x1000)..."<<flush;
	t = clock();
	RandomMatrix = CreateRandomMatrix(1000,1000);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestNormalization(){
	using namespace std;

	cout<<"Normalization Test"<<endl;
	cout<<"------------------"<<endl;

	BoostDoubleMatrix TestMatrix(3, 3);
	TestMatrix(0, 0) = 1; TestMatrix(0, 1) = 5; TestMatrix(0, 2) = 3;
	TestMatrix(1, 0) = 2; TestMatrix(1, 1) = 11; TestMatrix(1, 2) = 5;
	TestMatrix(2, 0) = 1; TestMatrix(2, 1) = 10; TestMatrix(2, 2) = 6;
	cout<<prefix<<"Original Matrix: " << TestMatrix << endl;
	

	clock_t t = clock();
	BoostDoubleMatrix NormMatrix = NormalizeMatrix(TestMatrix);
	t = clock() - t;
	cout<<prefix<<"Normalized Matrix: "<<NormMatrix<<" "<<ReportTime(t)<< endl;

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestImageMagick(char* test_image){
	using namespace std;
	using namespace Magick;
	clock_t t;

	cout<<"ImageMagick Test"<<endl;
	cout<<"----------------"<<endl;	
	
	cout<<prefix<<"Initializing ImageMagic..."<<flush;
	t = clock();
	InitializeMagick("");
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;


	cout<<prefix<<"Reading image ("<<test_image<<")..."<<flush;
	t = clock();
	Image someImage;
	someImage.read(test_image);
	t = clock() - t;
	//cout<<"done. "<<ReportTime(t)<<endl;
	cout<<"done. ["<<someImage.rows()<<"x"<<someImage.columns()<<"]. "<<ReportTime(t)<<endl;


	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestMatrixIO(char* test_image, char* output_image){
	using namespace std;
	clock_t t;

	cout<<"Matrix I/O Test"<<endl;
	cout<<"---------------"<<endl;

	/* Tiny Resolution */
	cout<<prefix<<"Loading image file into matrix at [32x32]..."<<flush;
	t = clock();
	BoostDoubleMatrix TinyImageMatrix = LoadImage(test_image,32,32);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Saving to file from image matrix at [32x32]..."<<flush;
	t = clock();
	WriteImage(TinyImageMatrix,output_image);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	/* Medium Resolution */
	cout<<prefix<<"Loading image file into matrix at [512,512]..."<<flush;
	t = clock();
	BoostDoubleMatrix MediumImageMatrix = LoadImage(test_image,512,512);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Saving to file from image matrix at [512x512]..."<<flush;
	t = clock();
	WriteImage(MediumImageMatrix,output_image);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	/* High Resolution */
	cout<<prefix<<"Loading image file into matrix at [2048,2048]..."<<flush;
	t = clock();
	BoostDoubleMatrix LargeImageMatrix = LoadImage(test_image,2048,2048);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Saving to file from image matrix at [2048x2048]..."<<flush;
	t = clock();
	WriteImage(LargeImageMatrix,output_image);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestNeighborCheck(){
	using namespace std;
	cout<<"Neighbor Index Test"<<endl;
	cout<<"-------------------"<<endl;
	cout<<prefix<<"Testing for a [2x2] Matrix"<<endl;

	cout<<prefix<<"[0] Neighbors: H="<<RightNeighbor(0,2)<<", V="<<DownNeighbor(0,2)<<endl;
	cout<<prefix<<"[1] Neighbors: H="<<RightNeighbor(1,2)<<", V="<<DownNeighbor(1,2)<<endl;
	cout<<prefix<<"[2] Neighbors: H="<<RightNeighbor(2,2)<<", V="<<DownNeighbor(2,2)<<endl;
	cout<<prefix<<"[3] Neighbors: H="<<RightNeighbor(3,2)<<", V="<<DownNeighbor(3,2)<<endl;

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestShrike(){
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

	cout<<"Shrinkage-Like Function Test"<<endl;
	cout<<"----------------------------"<<endl;
	cout<<prefix<<"Set g = "<<g<<endl;
	cout<<prefix<<"Set nu = "<<nu<<endl;

	cout<<prefix<<"Testing ShrikeAnisotropic @ beta = 0.5..."<<flush;
	t = clock();
	Result = ShrikeAnisotropic(g,nu,0.5);
	t = clock() - t;
	cout<<"done. ["<<Result<<"]. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Testing ShrikeAnisotropic @ beta = 10..."<<flush;
	t = clock();
	Result = ShrikeAnisotropic(g,nu,10);
	t = clock() - t;
	cout<<"done. ["<<Result<<"]. "<<ReportTime(t)<<endl;
	
	cout<<prefix<<"Testing ShrikeIsotropic @ beta = 0.5..."<<flush;
	t = clock();
	Result = ShrikeIsotropic(g,nu,0.5);
	t = clock() - t;
	cout<<"done. ["<<Result<<"]. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Testing ShrikeIsotropic @ beta = 0.65..."<<flush;
	t = clock();
	Result = ShrikeIsotropic(g,nu,0.65);
	t = clock() - t;
	cout<<"done. ["<<Result<<"]. "<<ReportTime(t)<<endl;

	BoostDoubleMatrix G  (3,2); 
		SetRow(G,g,0); 
		SetRow(G,g,1); 
		SetRow(G,g,2);
	BoostDoubleMatrix Nu (3,2); 
		SetRow(Nu,nu,0); 
		SetRow(Nu,nu,1); 
		SetRow(Nu,nu,2);
	cout<<prefix<<"Set G = "<<G<<endl;
	cout<<prefix<<"Set Nu = "<<Nu<<endl;

	cout<<prefix<<"Testing ApplyShrike(Anisotropic) @ beta = 10..."<<flush;
	t = clock();
	BoostDoubleMatrix Shriked = ApplyShrike(G,Nu,10,ANISOTROPIC);
	t = clock() - t;
	cout<<"done. ["<<Shriked<<"]. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Testing ApplyShrike(Isotropic) @ beta = 0.65..."<<flush;
	t = clock();
	Shriked = ApplyShrike(G,Nu,0.65,ISOTROPIC);
	t = clock() - t;
	cout<<"done. ["<<Shriked<<"]. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestRasterization(){
	using namespace std;
	/* Rasterization Test */
	BoostDoubleMatrix RasterTest(2,2);
	RasterTest(0,0) = 1;
	RasterTest(1,0) = 2;
	RasterTest(0,1) = 3;
	RasterTest(1,1) = 4;
	clock_t t;

	cout<<"Matrix Rasterization Test"<<endl;
	cout<<"-------------------------"<<endl;

	cout<<prefix<<"Set Matrix as ["<<RasterTest<<"]"<<endl;
	cout<<prefix<<"Rasterized Matrix: "<<MatrixToVector(RasterTest)<<endl;
	cout<<prefix<<"Restored Matrix: "<<VectorToMatrix(MatrixToVector(RasterTest),2,2)<<endl;
	
	BoostDoubleMatrix LargeMatrix(2048,2048);
	cout<<prefix<<"Rasterizing large matrix [2048x2048]..."<<flush;
	t = clock();
	BoostDoubleVector LargeVector = MatrixToVector(LargeMatrix);
	t = clock()-t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Restoring large matrix [2048x2048]..."<<flush;
	t = clock();
	RasterTest = VectorToMatrix(LargeVector,2048,2048);
	t = clock()-t;
	cout<<"done. "<<ReportTime(t)<<endl;	

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestLagrangian(){
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
	W(0, 0) =  -1; W(0, 1) = -1;
	W(1, 0) =  -1; W(1, 1) = -1;
	W(2, 0) =  -1; W(2, 1) = -1;
	W(3, 0) =  -1; W(3, 1) = -1;

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
	
	/*Test Lagrangian function*/	
	cout<<"Lagrangian Test"<<endl;
	cout<<"---------------"<<endl;

	cout<<prefix<<"Set A = "<<A<<endl;
	cout<<prefix<<"Set U = "<<U<<endl;
	cout<<prefix<<"Set B = "<<B<<endl;
	cout<<prefix<<"Set W = "<<W<<endl;
	cout<<prefix<<"Set Nu = "<<Nu<<endl;
	cout<<prefix<<"Set Lambda = "<<Lambda<<endl;
	cout<<prefix<<"Set beta = "<<beta<<endl;
	cout<<prefix<<"Set mu = "<<mu<<endl;

	cout<<prefix<<"Calculating Lagrangian..."<<flush;
	t = clock();
	double L = LagrangianNew(A, U, B, W, Nu, Lambda, beta, mu,2,ISOTROPIC);
	t = clock() - t;
	cout<<"done. ["<<L<<"]. [247.1569...] Expected. "<<ReportTime(t)<<endl;



	A = BoostDoubleMatrix(64*64,128*128);
	W = BoostDoubleMatrix(128*128,2);
	Nu = BoostDoubleMatrix(128*128,2);
	B = BoostDoubleVector(64*64);
	Lambda = BoostDoubleVector(64*64);
	U = BoostDoubleVector(128*128);

	cout<<prefix<<"Calculating Lagrangian(ISO) on large dataset..."<<flush;
	t = clock();
	L = LagrangianNew(A, U, B, W, Nu, Lambda, beta, mu,128,ISOTROPIC);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Calculating Lagrangian(ANISO) on large dataset..."<<flush;
	t = clock();
	L = LagrangianNew(A, U, B, W, Nu, Lambda, beta, mu,128,ANISOTROPIC);
	t = clock() - t;
	cout<<"done. "<<ReportTime(t)<<endl;

	cout<<prefix<<"Passed."<<endl<<endl;
}

void TestOnestep_Direction() {
	using namespace std;

	BoostDoubleMatrix Ad(3, 4), Wd(4, 2), NUd(4, 2);
	BoostDoubleVector Ud(4), Bd(3), LAMBDAd(3);
	double betad = 0.25;
	double mud = 0.5;
	clock_t t;
	
	/* Init Projections */
	Ad(0, 0) = 1; Ad(0, 1) = 1; Ad(0, 2) = 1; Ad(0, 3) = 1;
	Ad(1, 0) = 2; Ad(1, 1) = 2; Ad(1, 2) = 2; Ad(1, 3) = 2;
	Ad(2, 0) = 3; Ad(2, 1) = 3; Ad(2, 2) = 3; Ad(2, 3) = 3;
	
	/* Init Gradient Dual Variables  */
	Wd(0, 0) = -1; Wd(0, 1) = -1;
	Wd(1, 0) = -1; Wd(1, 1) = -1;
	Wd(2, 0) = -1; Wd(2, 1) = -1;
	Wd(3, 0) = -1; Wd(3, 1) = -1;

	/* Init Dual Lag. Multipliers */
	NUd(0, 0) = 0.25; NUd(0, 1) = 0.25;
	NUd(1, 0) = 0.25; NUd(1, 1) = 0.25;
	NUd(2, 0) = 0.25; NUd(2, 1) = 0.25;
	NUd(3, 0) = 0.25; NUd(3, 1) = 0.25;

	/* Init Image Vector */
	Ud(0) = 1;
	Ud(1) = 2;
	Ud(2) = 3;
	Ud(3) = 4;

	/* Init Observations */
	Bd(0) = 3;
	Bd(1) = 3;
	Bd(2) = 3;

	/* Init Observation Lag. Multipliers */
	LAMBDAd(0) = 0.5;
	LAMBDAd(1) = 0.5;
	LAMBDAd(2) = 0.5;

	/*Test One-step direction*/
	cout << "One-step direction Test" << endl;
	cout << "-----------------------" << endl;

	cout << prefix << "Set A = " << Ad << endl;
	cout << prefix << "Set U = " << Ud << endl;
	cout << prefix << "Set B = " << Bd << endl;
	cout << prefix << "Set W = " << Wd << endl;
	cout << prefix << "Set Nu = " << NUd << endl;
	cout << prefix << "Set Lambda = " << LAMBDAd << endl;
	cout << prefix << "Set beta = " << betad << endl;
	cout << prefix << "Set mu = " << mud << endl;

	cout << prefix << "Calculating direction value... [58.75  58.25  57.75  57.25] Expected. " << flush;
	t = clock();
	BoostDoubleVector D = Onestep_Direction(Ad, Ud, Bd, Wd, NUd, LAMBDAd, betad, mud, 2);
	t = clock() - t;
	cout << "done. [" << D << "]." << ReportTime(t) << endl;

	cout << prefix << "Passed." << endl << endl;
}

int main(int argc, char **argv){
	using namespace std;
	cout<<endl;

	TestRasterization();
	TestRandomMatrix();
	TestNormalization();
	TestNeighborCheck();
	TestShrike();
	TestLagrangian();
	TestOnestep_Direction();


	if(argc > 2){
	// Only run tests requiring File I/O if the file names
	// have been passed.
		TestImageMagick(argv[1]);
		TestMatrixIO(argv[1],argv[2]);
	}

	/* Gradient Test */
	BoostDoubleMatrix GradTest(2,2);
	GradTest(0,0) = 1;
	GradTest(1,0) = 2;
	GradTest(0,1) = 3;
	GradTest(1,1) = 4;
	std::cout<<"Gradient Test on A = [1 3; 2 4]"<<std::endl;
	std::cout<<"  "<<AllPixelGradients(MatrixToVector(GradTest),2)<<std::endl;

	/* Gradient Adjoint Test */
	// We should have a set of pixel gradients from A = [1 3; 2 4]
	// which look like
	// G = [-2 -1; -2 0; 0 -1; 0 0]
	// then the correct answer is a sum vector of
	// [-3, -1, 1, 3]
	BoostDoubleMatrix thisGradient = AllPixelGradients(MatrixToVector(GradTest),2);
	std::cout<<"Adjoint Sum Test"<<std::endl;
	std::cout<<"  Sum Vector = "<<PixelGradientAdjointSum(thisGradient,2)<<std::endl;
return 0;
}