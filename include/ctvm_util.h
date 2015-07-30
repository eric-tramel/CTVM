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
typedef boost::numeric::ublas::zero_matrix<double> BoostZeroMatrix;
typedef boost::numeric::ublas::vector<double> BoostDoubleVector;
typedef boost::numeric::ublas::scalar_vector<double> BoostScalarDoubleVector;
typedef boost::numeric::ublas::zero_vector<double> BoostZeroVector;

double MaximumEntry(BoostDoubleMatrix AMatrix);
double MaximumEntry(BoostDoubleVector AVector);

double MinimumEntry(BoostDoubleMatrix AMatrix);
double MinimumEntry(BoostDoubleVector AVector);

BoostDoubleVector HadamardProduct(BoostDoubleVector A, BoostDoubleVector B);
BoostDoubleVector SignVector(BoostDoubleVector AVector);
BoostDoubleVector AbsoluteValueVector(BoostDoubleVector AVector);
BoostDoubleVector MaxVector(BoostDoubleVector A, BoostDoubleVector B);
BoostDoubleVector MaxVector(BoostDoubleVector A, double B);

BoostDoubleVector MakeUnitVector(BoostDoubleVector AVector);
BoostDoubleMatrix NormalizeMatrix(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix CreateRandomMatrix(int rows, int cols);
BoostDoubleVector CreateRandomVector(int length);

BoostDoubleMatrix ImageToMatrix(Magick::Image AnImage);
BoostDoubleMatrix LoadImage(const char* ImageFileName);
BoostDoubleMatrix LoadImage(const char* ImageFileName, int newRows, int newCols);

void WriteImage(BoostDoubleMatrix AMatrix, const char* OutputFile);
BoostDoubleVector MatrixToVector(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix VectorToMatrix(BoostDoubleVector AVector,unsigned long rows, unsigned long cols);
BoostDoubleVector ReadTiltAngles(char* TiltAngleFile);

#endif