#ifndef FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_BLACK_SCHOLE_H_
#define FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_BLACK_SCHOLE_H_

#include <cmath>

namespace FinancialMathematics {
  double d1(double t, double x, double K, double mu, double sigma, double T);
  double d2(double t, double x, double K, double mu, double sigma, double T);

  double norm_cdf(double value);

  double BS_call_price(double t, double x, double K, double mu, double sigma, double T);
  double BS_put_price(double t, double x, double K, double mu, double sigma, double T);
}

#endif //FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_BLACK_SCHOLE_H_
