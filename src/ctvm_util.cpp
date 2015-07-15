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

    /* Set Namespace */
    // using namespace boost::numeric::ublas;

    /* Matrix Allocation */
    BoostDoubleMatrix RandomMatrix (rows,cols);

    /* Random Number Generation */
    boost::posix_time::ptime CurrentTick = boost::posix_time::microsec_clock::local_time();
    boost::mt19937 RNG (static_cast<unsigned int>(CurrentTick.time_of_day().total_milliseconds()));     // Set the integer RNG engine. Time Based Seed.
    boost::normal_distribution<> NormDist (0.0,1.0);   // Specify the distribution with zero mean and unit variance
    boost::variate_generator<boost::mt19937&,
                             boost::normal_distribution<> > RandNormValue(RNG,NormDist);

    /* Loop and Assign */
    for(unsigned i=0; i < RandomMatrix.size1(); ++i){
        for(unsigned j=0; j < RandomMatrix.size2(); ++j){
            RandomMatrix(i,j) = RandNormValue();
        }
    }

    /* Return Filled Random Matrix */
    return RandomMatrix;
}