#include <FinancialMathematics/brownian_motion.h>

#include <cmath>
#include <cstdint>

#include <omp.h>

namespace FinancialMathematics {
  AssetPaths BrownianMotion(const EtaMatrix& eta, const Config& config, const std::function<double(double, uint32_t)>& eta_op) {
    const uint32_t n = eta.steps();
    const uint32_t m = eta.paths();

    const double dt = config.T / static_cast<double>(n);
    const double dt_sqrt = std::sqrt(dt);

    const double precalc1 = (config.mu - 0.5 * config.sigma * config.sigma) * dt;
    const double precalc2 = config.sigma * dt_sqrt;

    std::vector<std::vector<double>> x(m, std::vector<double>(n + 1));

#pragma omp parallel for default(none) shared(x, eta, eta_op, config, dt, dt_sqrt, n, m, precalc1, precalc2) schedule(static)
    for (int j = 0; j < m; ++j) {
      x[j][0] = config.x0;
      for (uint32_t i = 1; i <= n; ++i) {
        x[j][i] = x[j][i - 1] * std::exp(precalc1 + precalc2 * eta_op(eta[j][i - 1], i - 1));
      }
    }

    return AssetPaths(x);
  }

  AssetPaths BrownianMotionAntitheticEta(const EtaMatrix& eta, const Config& config, const std::function<double(double, uint32_t)>& eta_op) {
    const uint32_t n = eta.steps();
    const uint32_t m = eta.paths();

    const double dt = config.T / static_cast<double>(n);
    const double dt_sqrt = std::sqrt(dt);

    const double precalc1 = (config.mu - 0.5 * config.sigma * config.sigma) * dt;
    const double precalc2 = config.sigma * dt_sqrt;

    std::vector<std::vector<double>> x(m, std::vector<double>(n + 1));

#pragma omp parallel for default(none) shared(x, eta, eta_op, config, dt, dt_sqrt, n, m, precalc1, precalc2) schedule(static)
    for (int j = 0; j < m; ++j) {
      x[j][0] = config.x0;
      for (uint32_t i = 1; i <= n; ++i) {
        x[j][i] = x[j][i - 1] * std::exp(precalc1 - precalc2 * eta_op(eta[j][i - 1], i - 1));
      }
    }

    return AssetPaths(x);
  }
}