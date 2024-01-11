#include <FinancialMathematics/FinancialMathematics.h>

#include <iostream>
#include <iomanip>

using namespace FinancialMathematics;

int main() {
  SetSeed(13);

  double t = 0.0;
  Config config;
  config.x0 = 100;
  config.mu = 0.06;
  config.sigma = 0.4;
  config.T = 0.5;

  uint32_t paths = 1000;
  uint32_t steps = 500;

  EtaMatrix eta = EtaMatrix(paths, steps);

  AssetPaths asset_paths = BrownianMotion(eta, config);

  std::vector<double> x0 = {80, 100, 120, 140, 160};

  std::cout << "Values of American put option with mu = " << config.mu << ", sigma = " << config.sigma << ", T = "
            << config.T << ", K = 100" << std::endl;
  std::cout << "---------------------------------------------------------------------------------------------"
            << std::endl;
  std::cout << std::setw(20) << std::setprecision(2) << std::fixed << "x0";
  for (auto x : x0) {
    std::cout << std::setw(10) << std::setprecision(2) << std::fixed << x << " ";
  }
  std::cout << std::endl;

  std::cout << std::setw(20) << std::setprecision(2) << std::fixed << "Martingale sim";

  for (auto x : x0) {
    config.x0 = x;
    asset_paths = BrownianMotion(eta, config);
    std::cout << std::setw(10) << std::setprecision(2) << std::fixed << American::american_put_martingale_MC(
        asset_paths,
        config)
              << " ";
  }
  std::cout << std::endl;

  std::cout << std::setw(20) << std::setprecision(2) << std::fixed << "Tree(100)";

  for (auto x : x0) {
    config.x0 = x;
    asset_paths = BrownianMotion(eta, config);
    std::cout << std::setw(10) << std::setprecision(2) << std::fixed
              << American::american_put_binomial_tree(config, 100, 100)
              << " ";
  }
  std::cout << std::endl;

  std::cout << std::setw(20) << std::setprecision(2) << std::fixed << "Tree(1000)";

  for (auto x : x0) {
    config.x0 = x;
    asset_paths = BrownianMotion(eta, config);
    std::cout << std::setw(10) << std::setprecision(2) << std::fixed
              << American::american_put_binomial_tree(config, 100, 1000)
              << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Values of European put option with mu = " << config.mu << ", sigma = " << config.sigma << ", T = "
            << config.T << ", K = 100" << std::endl;

  std::cout << "---------------------------------------------------------------------------------------------"
            << std::endl;

  std::cout << std::setw(20) << std::setprecision(2) << std::fixed << "x0";
  for (auto x : x0) {
    std::cout << std::setw(10) << std::setprecision(2) << std::fixed << x << " ";
  }
  std::cout << std::endl;
  std::cout << std::setw(20) << std::setprecision(2) << std::fixed << "BlackScholes";
  for (auto x : x0) {
    std::cout << std::setw(10) << std::setprecision(2) << std::fixed << BS_put_price(t, x, 100, config.mu, config.sigma,
                                                                                   config.T) << " ";
  }
  std::cout << std::endl;

  return 0;
}