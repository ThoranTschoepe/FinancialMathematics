#include <FinancialMathematics/black_scholes.h>

namespace FinancialMathematics {

  double d1(double t, double x, double K, double mu, double sigma, double T) {
    return (std::log(x / K + 1e-10) + (mu + 0.5 * sigma * sigma) * (T - t)) / (sigma * std::sqrt(T - t) + 1e-10);
  }

  double d2(double t, double x, double K, double mu, double sigma, double T) {
    return (std::log(x / K + 1e-10) + (mu - 0.5 * sigma * sigma) * (T - t)) / (sigma * std::sqrt(T - t) + 1e-10);
  }

  double norm_cdf(double value) {
    return 0.5 * erfc(-value * M_SQRT1_2);
  }

  double BS_call_price(double t, double x, double K, double mu, double sigma, double T) {
    return x * norm_cdf(d1(t, x, K, mu, sigma, T)) - K * std::exp(-mu * (T - t)) * norm_cdf(d2(t, x, K, mu, sigma, T));
  }

  double BS_put_price(double t, double x, double K, double mu, double sigma, double T) {
    return K * std::exp(-mu * (T - t)) * norm_cdf(-d2(t, x, K, mu, sigma, T)) - x * norm_cdf(-d1(t, x, K, mu, sigma, T));
  }
}