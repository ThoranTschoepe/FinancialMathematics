#include <FinancialMathematics/mc.h>
#include <FinancialMathematics/brownian_motion.h>
#include <FinancialMathematics/black_scholes.h>

#include <cmath>

#include <omp.h>
#include <numeric>
#include <utility>

namespace FinancialMathematics {
  double basicMonteCarlo(double t, const std::function<double(const std::vector<double> &)> &payoff,
                         const AssetPaths &asset_paths, const Config &config) {
    double sum = 0.0;
#pragma omp parallel for default(none) shared(asset_paths, payoff, config) reduction(+:sum)
    for (int i = 0; i < asset_paths.paths(); ++i) {
      sum += payoff(asset_paths[i]);
    }

    return std::exp(-config.mu * (config.T - t)) * sum / asset_paths.paths();
  }

  double MonteCarloAntithetic(double t, const std::function<double(const std::vector<double> &)> &payoff,
                              const AssetPaths &ap, const AssetPaths &ap_antithetic, const Config &config) {
    double sum = 0.0;

#pragma omp parallel for default(none) shared(ap, ap_antithetic, payoff, config) reduction(+:sum)
    for (int i = 0; i < ap.paths(); ++i) {
      sum += (payoff(ap[i]) + payoff(ap_antithetic[i])) / 2.0;
    }

    return std::exp(-config.mu * (config.T - t)) * sum / ap.paths();
  }

  double calc_control_variate_coefficient(const std::vector<double> &x, const std::vector<double> &y) {
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

    double cov = 0.0;
    double var = 0.0;

    for (int i = 0; i < x.size(); ++i) {
      cov += (x[i] - mean_x) * (y[i] - mean_y);
      var += (y[i] - mean_y) * (y[i] - mean_y);
    }

    return cov / var;
  }

  double MonteCarloControlVariateUnderlying(double t, const std::function<double(const std::vector<double> &)> &payoff,
                                            const AssetPaths &ap, const AssetPaths &ap_optimize, const Config &config) {
    std::vector<double> payoffs_optimize(ap_optimize.paths());
    std::vector<double> underlying_optimize(ap_optimize.paths());

    double precalc1 = std::exp(-config.mu * (config.T - t));
    double precalc2 = config.x0 * std::exp(config.mu * (config.T - t));

#pragma omp parallel for default(none) shared(ap_optimize, payoff, payoffs_optimize, underlying_optimize, precalc1)
    for (int i = 0; i < ap_optimize.paths(); ++i) {
      payoffs_optimize[i] = precalc1 * payoff(ap_optimize[i]);
      underlying_optimize[i] = ap_optimize[i].back();
    }

    double control_variate_coefficient = calc_control_variate_coefficient(payoffs_optimize, underlying_optimize);

    double sum = 0.0;

#pragma omp parallel for default(none) shared(t, ap, payoff, config, control_variate_coefficient, precalc2) reduction(+:sum)
    for (int i = 0; i < ap.paths(); ++i) {
      sum += payoff(ap[i]) - control_variate_coefficient * (ap[i].back() - precalc2);
    }

    return std::exp(-config.mu * (config.T - t)) * sum / ap.paths();
  }

  double MonteCarloImportanceSampling(double t, const std::function<double(const std::vector<double> &)> &payoff,
                                      double const_drift, double m_drift,
                                      const EtaMatrix &eta, const Config &config) {
    std::vector<double> drift(eta.steps(), const_drift);
    for (int i = 0; i < eta.paths(); ++i) {
      m_drift += m_drift * i;
    }

    auto g_h = [&](const std::vector<double> &x, const std::vector<double> &drift) {
      return std::exp(-std::inner_product(drift.begin(), drift.end(), x.begin(), 0.0) -
                      0.5 * std::inner_product(drift.begin(), drift.end(), drift.begin(), 0.0));
    };

    AssetPaths ap_drifted = BrownianMotion(eta, config, [&](double x, uint32_t i) { return x + drift[i]; });

    double sum = 0;

//#pragma omp parallel for default(none) shared(ap_drifted, payoff, sum, g_h, drift, config, t)
    for (int i = 0; i < ap_drifted.paths(); ++i) {
      sum += std::exp(-config.mu * (config.T - t)) * payoff(ap_drifted[i]) * g_h(eta[i], drift);
    }

    return sum / ap_drifted.paths();
  }

  namespace American {
    double golden_search(double a, double b, const double tol, const std::function<double(double)> &f) {
      double phi = (1 + std::sqrt(5)) / 2;
      double c = b - (b - a) / phi;
      double d = a + (b - a) / phi;
      double fc = f(c);
      double fd = f(d);

      while (std::abs(c - d) > tol) {
        if (fc < fd) {
          b = d;
          d = c;
          fd = fc;
          c = b - (b - a) / phi;
          fc = f(c);
        } else {
          a = c;
          c = d;
          fc = fd;
          d = a + (b - a) / phi;
          fd = f(d);
        }
      }

      return (c + d) / 2;
    }

    double american_put_martingale_MC(const AssetPaths &asset_paths, const Config &config) {
      uint32_t paths = asset_paths.paths();
      uint32_t steps = asset_paths.steps();

      // Initialize the 2D arrays
      std::vector<std::vector<double>> R(paths, std::vector<double>(steps, 0.0));
      std::vector<std::vector<double>> discounted_european_put_mt(paths, std::vector<double>(steps, 0.0));

      double K = 100;
      auto american_put = [&](double S) {
        return std::max(K - S, 0.0);
      };

// Fill the arrays
#pragma omp parallel for
      for (int i = 0; i < paths; ++i) {
        for (int j = 0; j < steps; ++j) {
          double t = static_cast<double>(j) * config.T / steps;
          double BS_put = BS_put_price(t, asset_paths[i][j], K, config.mu, config.sigma, config.T);

          R[i][j] = std::exp(-config.mu * t) * american_put(asset_paths[i][j]);
          discounted_european_put_mt[i][j] = std::exp(-config.mu * t) * BS_put;
        }
      }

      std::vector<double> sup(paths, 0);

      auto to_optimize = [&](double x) {
#pragma omp parallel for
        for (int i = 0; i < paths; ++i) {
          double max = -std::numeric_limits<double>::infinity();

          for (int j = 0; j < steps; ++j) {
            double val = R[i][j] + x * (discounted_european_put_mt[i].back() - discounted_european_put_mt[i][j]);
            if (val > max) {
              max = val;
            }
          }

          sup[i] = max;
        }

        return std::accumulate(sup.begin(), sup.end(), 0.0) / sup.size();
      };

      return to_optimize(golden_search(0.0, 100, 0.0001, to_optimize));
    }

    //Barone-Adesi and Whaley approximation

    double american_put_binomial_tree(const Config &config, double K, uint32_t n) {
      double deltaT = config.T / n;
      double up = std::exp(config.sigma * std::sqrt(deltaT));
      double down = 1 / up;
      double p = (std::exp(config.mu * deltaT) - down) / (up - down);
      double discount = std::exp(-config.mu * deltaT);
      double p0 = discount * p;
      double p1 = discount * (1 - p0);

      std::vector<double> values(n + 1);

      for (int i = 0; i <= n; ++i) {
        values[i] = std::max(0.0, K - config.x0 * std::pow(up, n - i) * std::pow(down, i));
      }

      for (int j = n - 1; j >= 0; --j) {
        for (int i = 0; i <= j; ++i) {
          // holding values vs exercise values
          values[i] = std::max(p0 * values[i] + p1 * values[i + 1], K - config.x0 * std::pow(up, j - i) * std::pow(down, i));
        }
      }

      return values[0];
    }

    //Korn Kreer lenssen (trinomial tree)
  }
}