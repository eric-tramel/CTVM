#include "ctvm_util.h"

BoostDoubleMatrix CreateRandomMatrix(int rows, int cols){
/*
* Function: CreateRandomMatrix
* ------------------------------
* Allocates a matrix of the specified dimensionality and fills it with random 
* entries drawn from a Normal distribution of variance 1.
* 
* rows: number of rows in generated matrix
* cols: number of columns in generated matrix
*
* return -- A random matrix of type `matrix<double>` 
*/
    /* Matrix Allocation */
    BoostDoubleMatrix RandomMatrix (rows,cols);

    /* Random Number Generation */
    boost::posix_time::ptime CurrentTick = boost::posix_time::microsec_clock::local_time();
    boost::mt19937 RNG (static_cast<unsigned int>(CurrentTick.time_of_day().total_milliseconds()));     // Set the integer RNG engine. Time Based Seed.
    boost::normal_distribution<> NormDist (0.0,1.0);   // Specify the distribution with zero mean and unit variance
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > RandNormValue(RNG,NormDist);    // Finally, build the number generator itself.

    /* Loop and Assign */
    for(unsigned i=0; i < RandomMatrix.size1(); ++i){
        for(unsigned j=0; j < RandomMatrix.size2(); ++j){
            RandomMatrix(i,j) = RandNormValue();
        }
    }

    /* Return Filled Random Matrix */
    return RandomMatrix;
}

// BoostDoubleMatrix LoadImage(const char* ImageFileName,ImageFileType FileType){
// /*
// * Function: LoadImage
// * ----------------------------
// * Attempts to load an image file from the specified location using GIL. Subsequently,
// * the GIL image file is convereted to a matrix<double> and returned.
// *
// * return -- A random matrix of type `matrix<double>`.
// */
//     // using namespace boost::gil;

//     // rgb8_image_t Image;

//     // /* Attempt to Load Image */
//     // switch(FileType){
//     //     case JPG: 
//     //         jpeg_read_image(ImageFileName,Image);
//     //         break;
//     //     case TIFF:
//     //         tiff_read_image(ImageFileName,Image);
//     //         break;
//     //     case PNG:
//     //         png_read_image(ImageFileName,Image);
//     //         break;     
//     // }
    

//     Image RunTimeImage;
//     // try{
//         RunTimeImage.read(ImageFileName);        
//     // }
//     // catch(Exception &error_){
//     //     std::cout << "Caught exception: "<<error_.what() <<std::endl;
//     // }
    
//     BoostDoubleMatrix DBMImage (10,10);

// return DBMImage;
// }