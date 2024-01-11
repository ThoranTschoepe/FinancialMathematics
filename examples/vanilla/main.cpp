#include <FinancialMathematics/FinancialMathematics.h>

#include <iostream>

using namespace FinancialMathematics;

int main() {
  SetSeed(13);

  double t = 0.0;
  Config config;
  config.x0 = 1.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1;

  uint32_t paths = 100000;
  uint32_t steps = 1;

  EtaMatrix eta = EtaMatrix(paths, steps);
  EtaMatrix eta_opt = EtaMatrix(paths / 10, 1);
  AssetPaths asset_paths = BrownianMotion(eta, config);
  AssetPaths asset_paths_anti = BrownianMotionAntitheticEta(eta, config);
  AssetPaths asset_paths_opt = BrownianMotion(eta_opt, config);

  std::cout << "Asset paths:" << std::endl;
  std::cout << asset_paths << std::endl;

  double K = 1.1;
  auto payoff = [K](const std::vector<double> &x) {
    return std::max(x.back() - K, 0.0);
  };

    std::cout << "Values of European call option with x0 = " << config.x0 << ", mu = " << config.mu << ", sigma = "
                << config.sigma << ", T = " << config.T << ", K = " << K << std::endl;
  std::cout << basicMonteCarlo(t, payoff, asset_paths, config) << std::endl;
  std::cout << MonteCarloAntithetic(t, payoff, asset_paths, asset_paths_anti, config) << std::endl;
  std::cout << MonteCarloControlVariateUnderlying(t, payoff, asset_paths, asset_paths_opt, config) << std::endl;
  std::cout << "Analytical solution: " << BS_call_price(t, config.x0, K, config.mu, config.sigma, config.T) << std::endl;

  return 0;
}