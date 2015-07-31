#ifndef CTVM_UTIL_H
#define CTVM_UTIL_H
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Magick++.h>

/* UBLAS Typecasts */
typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
typedef boost::numeric::ublas::scalar_matrix<double> BoostScalarDoubleMatrix;
typedef boost::numeric::ublas::zero_matrix<double> BoostZeroMatrix;
typedef boost::numeric::ublas::vector<double> BoostDoubleVector;
typedef boost::numeric::ublas::scalar_vector<double> BoostScalarDoubleVector;
typedef boost::numeric::ublas::zero_vector<double> BoostZeroVector;

/* Matrix Manipulation */
BoostDoubleVector GetRow(BoostDoubleMatrix AMatrix,unsigned int row);
BoostDoubleVector GetCol(BoostDoubleMatrix AMatrix,unsigned int col);
void SetRow(BoostDoubleMatrix &AMatrix, BoostDoubleVector RowVect, unsigned int row);
void SetCol(BoostDoubleMatrix &AMatrix, BoostDoubleVector ColVect, unsigned int col);
BoostDoubleVector MatrixToVector(BoostDoubleMatrix AMatrix);
BoostDoubleMatrix VectorToMatrix(BoostDoubleVector AVector,unsigned long rows, unsigned long cols);

/* Matrix Search */
double MaximumEntry(BoostDoubleMatrix AMatrix);
double MaximumEntry(BoostDoubleVector AVector);
double MinimumEntry(BoostDoubleMatrix AMatrix);
double MinimumEntry(BoostDoubleVector AVector);

/* Linear Operations */
BoostDoubleVector HadamardProduct(BoostDoubleVector A, BoostDoubleVector B);
BoostDoubleVector SignVector(BoostDoubleVector AVector);
BoostDoubleVector AbsoluteValueVector(BoostDoubleVector AVector);
BoostDoubleVector MaxVector(BoostDoubleVector A, BoostDoubleVector B);
BoostDoubleVector MaxVector(BoostDoubleVector A, double B);
BoostDoubleMatrix NormalizeMatrix(BoostDoubleMatrix AMatrix);

/* Linear Algebra */
BoostDoubleVector MakeUnitVector(BoostDoubleVector AVector);

/* Matrix Generation */
BoostDoubleMatrix CreateRandomMatrix(int rows, int cols);
BoostDoubleVector CreateRandomVector(int length);

/* Image Operations */
unsigned int RightNeighbor(BoostDoubleVector ImageVector,unsigned int index, unsigned int SideLength);
unsigned int DownNeighbor(BoostDoubleVector ImageVector,unsigned int index, unsigned int SideLength);

/* File I/O */
// Image
BoostDoubleMatrix ImageToMatrix(Magick::Image AnImage);
BoostDoubleMatrix LoadImage(const char* ImageFileName);
BoostDoubleMatrix LoadImage(const char* ImageFileName, int newRows, int newCols);
void WriteImage(BoostDoubleMatrix AMatrix, const char* OutputFile);
// Raw Data
BoostDoubleVector ReadTiltAngles(char* TiltAngleFile);

#endif