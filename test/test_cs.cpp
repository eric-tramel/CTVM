#include <iostream>
#include <cmath>
#include "ctvm.h"
#include "ctvm_util.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Magick++.h>

int main(int argc, char **argv)
{
    using namespace std;
    if(argc != 4){
        cout<<"Usage: test-cs <imsize> <orig_image> <output_image>"<<endl;
        return 0;
    }

    /* Read Input */
    char* ImageSizeStr = argv[1];
    char* OriginalImageFile = argv[2];
    char* OutputImageFile = argv[3];
    unsigned int L = atoi(ImageSizeStr);

    /* Specify Problem Settings */
    double MeasurementRate = 0.9;
    double NoiseVariance = 0.0;
    unsigned int N = L*L;
    unsigned int M = MeasurementRate*N;     // Allow truncation

    /* Load Image */
    BoostDoubleMatrix XImage = LoadImage(OriginalImageFile,L,L);
    BoostDoubleVector XVect = MatrixToVector(XImage);

    /* Create Projection Matrix */
    BoostDoubleMatrix A = CreateRandomMatrix(M,N);

    /* Create Measurements */
    BoostDoubleVector y = prod(A,XVect) + sqrt(NoiseVariance)*CreateRandomVector(M);

    /* Perform Reconstruction */
    // Just testing a silly transpose operation for the moment
    BoostDoubleVector XRecVect = prod(trans(A),y);
    BoostDoubleMatrix XRecImage = VectorToMatrix(XRecVect,L,L);

    /* Write Result */
    WriteImage(NormalizeMatrix(XRecImage),OutputImageFile);

return 0;
}