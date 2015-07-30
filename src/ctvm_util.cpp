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
    boost::normal_distribution<> NormDist (0.0,1.0);  // Specify the distribution with zero mean and unit variance
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > RandNormValue(RNG,NormDist);    // Finally, build the number generator itself.

    /* Loop and Assign */
    for(unsigned long i=0; i < RandomMatrix.size1(); ++i){
        for(unsigned long j=0; j < RandomMatrix.size2(); ++j){
            RandomMatrix(i,j) = RandNormValue();
        }
    }

    /* Return Filled Random Matrix */
    return RandomMatrix;
}

BoostDoubleVector CreateRandomVector(int length){
/*
* Function: CreateRandomVector
* ------------------------------
* Allocates a vector of the specified length and fills it with random 
* entries drawn from a Normal distribution of variance 1.
* 
* length: number of elements in the generated vector
*
* return -- A random vector of type BoostDoubleVector
*/
    /* Matrix Allocation */
    BoostDoubleVector RandomVector (length);

    /* Random Number Generation */
    boost::posix_time::ptime CurrentTick = boost::posix_time::microsec_clock::local_time();
    boost::mt19937 RNG (static_cast<unsigned int>(CurrentTick.time_of_day().total_milliseconds()));     // Set the integer RNG engine. Time Based Seed.
    boost::normal_distribution<> NormDist (0.0,1.0);  // Specify the distribution with zero mean and unit variance
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > RandNormValue(RNG,NormDist);    // Finally, build the number generator itself.

    /* Loop and Assign */
    for(unsigned long i=0; i < RandomVector.size(); ++i){
        RandomVector(i) = RandNormValue();
    }

    /* Return Filled Random Matrix */
    return RandomVector;
}

BoostDoubleMatrix ImageToMatrix(Magick::Image AnImage){
/*
* Function: ImageToMatrix
* ---------------------------
* Convert from the Magick++ defined Image type to the BoostDoubleMatrix
* type.
*/
    using namespace Magick;
    
    /* Image Properties */
    int rows = AnImage.rows();
    int cols = AnImage.columns();

    /* Convert to Grayscale */
    AnImage.quantizeColorSpace( GRAYColorspace );
    AnImage.quantizeColors(256);
    AnImage.quantize();

    /* Assignment to BoostDoubleMatrix */
    BoostDoubleMatrix DBMImage (rows,cols);     // Allocate the matrix
    Pixels View(AnImage);  // Specify the "view" into the image
    PixelPacket *AllPixels = View.get(0,0,cols,rows);   // Get the pointer to the pixel array

    /* Loop over Pixels */
    for(int i=0; i<rows; ++i){
        for(int j=0; j< cols; ++j){
            // The entire image should be cast into grayscale, so all of the
            // RGB components should be equal to the same values. Hence, we
            // take the red pixel value.
            // Additionally, it is necessary to increment the pixel pointer.
            Quantum thisPixel = (*AllPixels++).red;
            // Finally, we have to cast the Quantum type into a double.
            DBMImage(i,j) = ColorGray::scaleQuantumToDouble(thisPixel);
        }
    }

return DBMImage;
}

BoostDoubleMatrix LoadImage(const char* ImageFileName){
/*
* Function: LoadImage
* ----------------------------
* Attempts to load an image file from the specified location using ImageMagick. 
* Subsequently,
* the ImageMagick image file is convereted to a matrix<double> and returned.
*
* return -- A random matrix of type `matrix<double>`.
*/
    using namespace Magick;
    

    Image RunTimeImage;
    try{
        RunTimeImage.read(ImageFileName);        
    }
    catch(Exception &error_){
        std::cout << "Caught exception: "<<error_.what() <<std::endl;
    }    

return ImageToMatrix(RunTimeImage);
}

BoostDoubleMatrix LoadImage(const char* ImageFileName, int newRows, int newCols){
/*
* Function: LoadImage
* ----------------------------
* Attempts to load an image file from the specified location using ImageMagick. 
* Subsequently,
* the ImageMagick image file is convereted to a BoostMatrixDouble and returned.
*
* return -- A random matrix of type BoostMatrixDouble.
*/
    using namespace Magick;
    

    Image RunTimeImage;
    try{
        RunTimeImage.read(ImageFileName);        
    }
    catch(Exception &error_){
        std::cout << "Caught exception: "<<error_.what() <<std::endl;
    }    

    /* Attept Image Resize */
    // RunTimeImage = AdaptiveResizeImage(RunTimeImage,newRows,newCols);
    RunTimeImage.adaptiveResize(Geometry(newCols,newRows));

return ImageToMatrix(RunTimeImage);   
}

void WriteImage(BoostDoubleMatrix AMatrix, const char* OutputFile){
/*
* Function: WriteImage
* ----------------------------
* Given a UBlas matrix as well as an output file name, attempt to write
* that matrix as an image. This assumes that the matrix is already scaled
* in the range [0,1].
*/
    using namespace Magick;
    
    int rows = AMatrix.size1();
    int cols = AMatrix.size2();

    /* Create Image */
    Image OutputImage;                                // Specify the output file
    OutputImage.size(Geometry(cols,rows));            // Specify image dimensions
    Pixels OutputImageView(OutputImage);

    /* Pixel Assignment Loop */
    // This procedure is, I'm sure, one of the least efficient ways
    // of accomplishing this task.
    PixelPacket thisPixel;
    for(int i = 0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            thisPixel.red = ColorGray::scaleDoubleToQuantum(AMatrix(i,j));
            thisPixel.green = ColorGray::scaleDoubleToQuantum(AMatrix(i,j));
            thisPixel.blue = ColorGray::scaleDoubleToQuantum(AMatrix(i,j));

            *(OutputImageView.get(j,i,1,1)) = thisPixel;
        }
    }

    /* Final Update */
    OutputImage.type(GrayscaleType);
    OutputImageView.sync();

    /* Write to Disk */
    OutputImage.write(OutputFile);
}

BoostDoubleMatrix NormalizeMatrix(BoostDoubleMatrix AMatrix){
/*
* Function: NormalizeMatrix
* ----------------------------
* Given a UBlas matrix, normalize it within the range [0,1].
*/  
    double AMatrixMin = MinimumEntry(AMatrix);
    double AMatrixMax = MaximumEntry(AMatrix);
    double AMatrixRange = AMatrixMax - AMatrixMin;

    BoostDoubleMatrix AMatrixNormalized(AMatrix);
    if(AMatrixRange != 0){
        // Boost/UBLAS does not support operations of a scalar against a matrix
        // natively.
        AMatrixNormalized -= BoostScalarDoubleMatrix(AMatrix.size1(),AMatrix.size2(),AMatrixMin);
        AMatrixNormalized /= AMatrixRange;
    }else{
        AMatrixNormalized = BoostDoubleMatrix(AMatrix.size1(),AMatrix.size2(),1);
    }

return AMatrixNormalized;
}



BoostDoubleVector MatrixToVector(BoostDoubleMatrix AMatrix){
/*
* Function: MatrixToVector
* ----------------------------
* Take a matrix and rasterize into a vector in a column-by-column 
* manner.
*
*/
    unsigned long N = AMatrix.size1()*AMatrix.size2();

    BoostDoubleVector AVector (N);
    
    unsigned long VectorIndex = 0;
    for(unsigned long i = 0; i < AMatrix.size1(); ++i){
        for(unsigned long j = 0; j < AMatrix.size2(); ++j){
            AVector(VectorIndex++) = AMatrix(i,j);
        }
    }

return AVector;
}

BoostDoubleMatrix VectorToMatrix(BoostDoubleVector AVector,unsigned long rows, unsigned long cols){
/*
* Function: VectorToMatrix
* ----------------------------
* Take a vector and convert it to a matrix given the specified 
* dimensions.
*
*/
    BoostDoubleMatrix AMatrix (rows,cols);    

    unsigned long VectorIndex = 0;
    for(unsigned long i = 0; i < AMatrix.size1(); ++i){
        for(unsigned long j = 0; j < AMatrix.size2(); ++j){
            AMatrix(i,j) = AVector(VectorIndex++);
        }
    }
return AMatrix;
}

BoostDoubleVector ReadTiltAngles(char* TiltAngleFile){
    using namespace std;
    BoostDoubleVector TiltAngles (1000);

    /* Define File Buffer */
    filebuf TAFileBuffer;
    
    /* Attempt to Open file */
    if(TAFileBuffer.open(TiltAngleFile,ios::in)){
        istream TAInputStream(&TAFileBuffer);

        /* Loop through file */
        double thisTiltAngle;
        unsigned int VectorIndex = 0;
        while(TAInputStream){
            TAInputStream>>thisTiltAngle;
            TiltAngles(VectorIndex++) = thisTiltAngle;
        }

        /* Check for Double Count */
        // It can happen that, if we have an additional new lines at the
        // end of the datafile we can end up double-counting the last entry
        // of the TiltAngle series. Since we should have no duplicates in
        // the list of angle files, anyway, we will simply remove this
        // repeated value.
        if(TiltAngles(VectorIndex-1) == TiltAngles(VectorIndex-2)){
            TiltAngles.resize(VectorIndex-1,1);
        }else{
            TiltAngles.resize(VectorIndex,1);
        }

        TAFileBuffer.close();
    }else{
        /* File Not Found */
        cout<<"Specified tilt angle file ("<<TiltAngleFile<<") could not be opened."<<endl;
    }

return TiltAngles;
}


double MaximumEntry(BoostDoubleMatrix AMatrix){
/*  Function: MaximumEntry
    ----------------------
    Return the maximal value in a given BoostDoubleMatrix
*/    
    double MaxValue = 0;
    for(unsigned int i=0;i<AMatrix.size1();++i){
        for(unsigned int j=0;j<AMatrix.size2();++j){
            MaxValue = std::max(MaxValue,AMatrix(i,j));
        }
    }
return MaxValue;
}

double MaximumEntry(BoostDoubleVector AVector){
/*  Function: MaximumEntry
    ----------------------
    Return the maximal value in a given BoostDoubleVector
*/    
    double MaxValue = 0;
    for(unsigned int i=0;i<AVector.size();++i){
        MaxValue = std::max(MaxValue,AVector(i));
    }
return MaxValue;
}


double MinimumEntry(BoostDoubleMatrix AMatrix){
/*  Function: MinimumEntry
    ----------------------
    Return the minimum value in a given BoostDoubleMatrix
*/    
    double MinValue = std::numeric_limits<double>::infinity();
    for(unsigned int i=0;i<AMatrix.size1();++i){
        for(unsigned int j=0;j<AMatrix.size2();++j){
            MinValue = std::min(MinValue,AMatrix(i,j));
        }
    }
return MinValue;
}

double MinimumEntry(BoostDoubleVector AVector){
/*  Function: MinimumEntry
    ----------------------
    Return the minimum value in a given BoostDoubleVector
*/    
    double MinValue = std::numeric_limits<double>::infinity();
    for(unsigned int i=0;i<AVector.size();++i){
            MinValue = std::min(MinValue,AVector(i));
    }
return MinValue;
}