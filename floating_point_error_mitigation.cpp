#include "my_solver.h"

#define N 25

void floatingPointErrorMitigation() {
  // e^(-8.3)
  double x = -8.3;
  std::cout << "*** f<-8.3> = " << std::scientific << taylorSeries(x, N) << std::endl;
  std::cout << "*** f<-8.3> = " << std::scientific << taylorSeriesImprovement(x, N) << std::endl;
  std::cout << "*** f<-8.3> = " << std::scientific << std::pow(2.71828182846, x) << std::endl;
}

double taylorSeries(double x, int n) {
  // Initialize sum of series
  double sum = 1.0;
  for (int i = n - 1; i > 0; --i)
    sum = 1.0 + x * sum / (double) i;

  return sum;
}

// Floating-point error mitigation - avoid subtraction
double taylorSeriesImprovement(double x, int n) {
  if (x < 0)
    return 1 / taylorSeries(-x, n);
  else
    return taylorSeries(x, n);
}
