#include <iostream>
#include <chrono>
#include <thread>
#include <ctime>
#include <cmath>
#include <random>

#define CHECK_TIME_START start = std::chrono::high_resolution_clock::now()
#define CHECK_TIME_END stop = std::chrono::high_resolution_clock::now(); \
                       duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
#define OUTPUT_DURATION(msg) std::cout << msg << " = " << (double) duration.count() / 1000 << "(s)" << std::endl
#define LINE_BREAK std::cout << "\n--------------------------------------\n" << std::endl

// Loop unrolling
void loopUnrolling();
void initMat();
void transposeMatrix(double *pMat, int MatSize);
void multiplySquareMatrices1(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);
void multiplySquareMatrices2(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);
void multiplySquareMatrices3(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);
void multiplySquareMatrices4(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);

// Polynomial evaluation
void polynomialEvaluation();
void initPolynomial();
void polynomialEvaluationOriginal(double *y, double *x, int n_x, double *a, int deg);
void polynomialEvaluationHorner(double *y, double *x, int n_x, double *a, int deg);
void checkDifference(double *y_n, double *y_h, int n_x);

// Floating-point error mitigation
void floatingPointErrorMitigation();
double taylorSeries(double x, int n);
double taylorSeriesImprovement(double x, int n);
