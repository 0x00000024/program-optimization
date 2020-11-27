#include "my_solver.h"

#define DEGREE 10
#define N_X 1048576

double *a, *x, *y_n, *y_h;

void polynomialEvaluation() {
  initPolynomial();

  std::chrono::high_resolution_clock::time_point start, stop;
  std::chrono::milliseconds duration;

  CHECK_TIME_START;
  polynomialEvaluationOriginal(y_n, x, N_X, a, DEGREE + 1);
  CHECK_TIME_END;
  OUTPUT_DURATION("Original method for polynomial evaluation");

  CHECK_TIME_START;
  polynomialEvaluationHorner(y_h, x, N_X, a, DEGREE + 1);
  CHECK_TIME_END;
  OUTPUT_DURATION("Horner's method for polynomial evaluation");

  checkDifference(y_n, y_h, N_X);

  delete[] a;
  delete[] x;
  delete[] y_n;
  delete[] y_h;
}

void initPolynomial() {
  a = new double[DEGREE + 1];
  x = new double[N_X];
  y_n = new double[N_X];
  y_h = new double[N_X];

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for (int i = 0; i < DEGREE + 1; ++i)
    a[i] = (double) dis(gen);
  for (int i = 0; i < N_X; ++i)
    x[i] = (double) dis(gen);
}

// std::pow
void polynomialEvaluationOriginal(double *y, double *x, int n_x, double *a, int deg) {
  int i, j;
  for (i = 0; i < n_x; ++i)
    for (j = 0, y[i] = 0.0; j <= deg; ++j)
      y[i] += a[j] * std::pow(x[i], j);
}

// Horner’s method for polynomial evaluation - reduce multiplication
void polynomialEvaluationHorner(double *y, double *x, int n_x, double *a, int deg) {
  int i, j;
  for (i = 0; i < n_x; ++i)
    for (j = deg - 1, y[i] = a[deg]; j >= 0; --j)
      y[i] = y[i] * x[i] + a[j];
}

void checkDifference(double *y_n, double *y_h, int n_x) {
  for (int i = 0; i < n_x; i++) {
    if ((float) y_n[i] != (float) y_h[i]) {
      printf("Value y_n[%d] (%.15lf) is not equal with y_h[%d] (%.15lf)\n\n", i, y_n[i], i, y_h[i]);
      return;
    }
  }
  std::cout << "All values are equal." << std::endl;
}
