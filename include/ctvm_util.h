#ifndef CTVM_UTIL_H
#define CTVM_UTIL_H
#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Magick++.h>

typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
typedef boost::numeric::ublas::vector<double> BoostDoubleVector;
enum ImageFileType {PNG,JPG,TIFF};

BoostDoubleMatrix CreateRandomMatrix(int rows, int cols);
BoostDoubleMatrix LoadImage(const char* ImageFileName);
BoostDoubleVector MatrixToVector(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix VectorToMatrix(BoostDoubleVector AVector,unsigned int rows, unsigned int cols);
BoostDoubleVector ReadTiltAngles(char* TiltAngleFile);

#endif