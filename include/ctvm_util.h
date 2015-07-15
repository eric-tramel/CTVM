#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;

BoostDoubleMatrix CreateRandomMatrix(int rows, int cols);