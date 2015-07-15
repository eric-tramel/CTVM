#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
// #include <Magick++.h>

typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
enum ImageFileType {PNG,JPG,TIFF};

BoostDoubleMatrix CreateRandomMatrix(int rows, int cols);
// BoostDoubleMatrix LoadImage(const char* ImageFileName,ImageFileType FileType);
