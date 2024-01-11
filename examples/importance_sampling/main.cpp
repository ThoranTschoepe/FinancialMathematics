#include <FinancialMathematics/FinancialMathematics.h>

#include <iostream>

using namespace FinancialMathematics;

int main() {
  SetSeed(13);

  double t = 0.0;
  Config config;
  config.x0 = 100;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1;

  uint32_t paths = 100000;
  uint32_t steps = 1;

  EtaMatrix eta = EtaMatrix(paths, steps);

  double K = 220;
  auto payoff = [K](const std::vector<double> &x) {
    return std::max(x.back() - K, 0.0);
  };

  AssetPaths asset_paths = BrownianMotion(eta, config);

  std::cout << basicMonteCarlo(t, payoff, asset_paths, config) << std::endl;
  std::cout << MonteCarloImportanceSampling(t, payoff, 5, 0, eta, config) << std::endl;

  std::cout << "Analytical solution: " << BS_call_price(t, config.x0, K, config.mu, config.sigma, config.T) << std::endl;

  return 0;
}