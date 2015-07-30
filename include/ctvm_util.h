#ifndef CTVM_UTIL_H
#define CTVM_UTIL_H
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Magick++.h>

typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
typedef boost::numeric::ublas::scalar_matrix<double> BoostScalarDoubleMatrix;
typedef boost::numeric::ublas::vector<double> BoostDoubleVector;
typedef boost::numeric::ublas::scalar_vector<double> BoostScalarDoubleVector;

enum ImageFileType {PNG,JPG,TIFF};

double MaximumEntry(BoostDoubleMatrix AMatrix);
double MinimumEntry(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix NormalizeMatrix(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix CreateRandomMatrix(int rows, int cols);
BoostDoubleMatrix LoadImage(const char* ImageFileName);
void WriteImage(BoostDoubleMatrix AMatrix, const char* OutputFile);
BoostDoubleVector MatrixToVector(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix VectorToMatrix(BoostDoubleVector AVector,unsigned long rows, unsigned long cols);
BoostDoubleVector ReadTiltAngles(char* TiltAngleFile);

#endif