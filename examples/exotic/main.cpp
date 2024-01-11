#include <FinancialMathematics/FinancialMathematics.h>

#include <iostream>
#include <numeric>

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
  uint32_t steps = 1000;

  EtaMatrix eta = EtaMatrix(paths, steps);
  EtaMatrix eta_opt = EtaMatrix(paths / 10, steps / 10);
  AssetPaths asset_paths = BrownianMotion(eta, config);
  AssetPaths asset_paths_anti = BrownianMotionAntitheticEta(eta, config);
  AssetPaths asset_paths_opt = BrownianMotion(eta_opt, config);

  std::cout << "Asset paths:" << std::endl;
  std::cout << asset_paths << std::endl;

  double K = 1.1;
  auto average = [K](const std::vector<double> &x) {
    double avg = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    return std::max(avg - K, 0.0);
  };

  std::cout << "Values of Asian call option with x0 = " << config.x0 << ", mu = " << config.mu << ", sigma = "
            << config.sigma << ", T = " << config.T << ", K = " << K << std::endl;
  std::cout << basicMonteCarlo(t, average, asset_paths, config) << std::endl;
  std::cout << MonteCarloAntithetic(t, average, asset_paths, asset_paths_anti, config) << std::endl;
  std::cout << MonteCarloControlVariateUnderlying(t, average, asset_paths, asset_paths_opt, config) << std::endl;

  //knock in option
  double B = 1.2;
  K = 1.1;
  auto knock_in = [K, B](const std::vector<double> &x) {
    bool knock_in = false;
    for (int i = 0; i < x.size(); ++i) {
      if (x[i] > B) {
        knock_in = true;
        break;
      }
    }

    if (knock_in) {
      return std::max(x.back() - K, 0.0);
    } else {
      return 0.0;
    }

  };

  std::cout << "Values of knock-in option with x0 = " << config.x0 << ", mu = " << config.mu << ", sigma = "
            << config.sigma << ", T = " << config.T << ", K = " << K << ", B = " << B << std::endl;
  std::cout << basicMonteCarlo(t, knock_in, asset_paths, config) << std::endl;
  std::cout << MonteCarloAntithetic(t, knock_in, asset_paths, asset_paths_anti, config) << std::endl;
  std::cout << MonteCarloControlVariateUnderlying(t, knock_in, asset_paths, asset_paths_opt, config) << std::endl;

  return 0;
}