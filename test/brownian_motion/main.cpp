#include <gtest/gtest.h>
#include <FinancialMathematics/FinancialMathematics.h>

using namespace FinancialMathematics;

double max_relative_error = 0.01;

TEST(MonteCarloMethods, TestBasicMonteCarlo) {
  SetSeed(13);
  // Define the parameters for the Monte Carlo method
  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  double K = 110.0;

  // Define the payoff function
  auto payoff = [K](const std::vector<double> &x) {
    return std::max(x.back() - K, 0.0);
  };

  // Generate the asset paths
  EtaMatrix eta = EtaMatrix(100000, 1);
  AssetPaths ap = BrownianMotion(eta, config);

  // Call the basicMonteCarlo method
  double result = basicMonteCarlo(t, payoff, ap, config);

  // compare to the analytical result BS_call_price
  double BS_call_price = FinancialMathematics::BS_call_price(t, config.x0, K, config.mu, config.sigma, config.T);

  //check relative error max 0.05
  ASSERT_LT(std::abs(result - BS_call_price) / BS_call_price, max_relative_error);
}

TEST(MonteCarloMethods, TestMonteCarloAntithetic) {
  SetSeed(13);
  // Define the parameters for the Monte Carlo method
  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  double K = 110.0;

  // Define the payoff function
  auto payoff = [K](const std::vector<double> &x) {
    return std::max(x.back() - K, 0.0);
  };

  // Generate the asset paths
  EtaMatrix eta = EtaMatrix(100000, 1);
  AssetPaths ap = BrownianMotion(eta, config);
  AssetPaths ap_antithetic = BrownianMotionAntitheticEta(eta, config);

  // Call the MonteCarloAntithetic method
  double result = MonteCarloAntithetic(t, payoff, ap, ap_antithetic, config);

  // compare to the analytical result BS_call_price
  double BS_call_price = FinancialMathematics::BS_call_price(t, config.x0, K, config.mu, config.sigma, config.T);

  //check relative error max 0.05
  ASSERT_LT(std::abs(result - BS_call_price) / BS_call_price, max_relative_error);
}

TEST(MonteCarloMethods, TestMonteCarloControlVariateUnderlying) {
  SetSeed(13);
  // Define the parameters for the Monte Carlo method
  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  double K = 110.0;

  // Define the payoff function
  auto payoff = [K](const std::vector<double> &x) {
    return std::max(x.back() - K, 0.0);
  };

  // Generate the asset paths
  uint32_t paths = 100000;
  uint32_t paths_optimize = paths/ 10;
  EtaMatrix eta = EtaMatrix(paths, 1);
  EtaMatrix eta_optimize = EtaMatrix(paths_optimize, 1);
  AssetPaths ap = BrownianMotion(eta, config);
  AssetPaths ap_optimize = BrownianMotion(eta_optimize, config);

  // Call the MonteCarloControlVariateUnderlying method
  double result = MonteCarloControlVariateUnderlying(t, payoff, ap, ap_optimize, config);

  // compare to the analytical result BS_call_price
  double BS_call_price = FinancialMathematics::BS_call_price(t, config.x0, K, config.mu, config.sigma, config.T);

  //check relative error max 0.05
  ASSERT_LT(std::abs(result - BS_call_price) / BS_call_price, max_relative_error);
}

#include <gtest/gtest.h>
#include <FinancialMathematics/FinancialMathematics.h>

using namespace FinancialMathematics;

TEST(FiniteDifferencesTest, TestFiniteDifferences) {
  SetSeed(13);

  double t = 0.0;
  Config config;
  config.x0 = 1.0;
  config.mu = 0;
  config.sigma = 1;
  config.T = 1;

  double K = 100;

  auto payoff = [K](double x) {
    return std::max(x - K, 0.0);
  };

  //finite differences
  int time_steps = 10000;
  int price_steps = 1000;
  constexpr bool quadratic_payoff = false;

  Eigen::MatrixXd pi = finiteDifferences<quadratic_payoff>(time_steps,
                                         price_steps,
                                         payoff,
                                         config.mu,
                                         config.sigma,
                                         config.T);

  //calculate analytical solution on grid
  Eigen::VectorXd ygrid = Eigen::VectorXd::LinSpaced(price_steps, std::log(0.0001), std::log(10000));
  Eigen::VectorXd analytical = Eigen::VectorXd::Zero(price_steps);

  for (int i = 0; i < price_steps; ++i) {
    analytical(i) = BS_call_price(t, std::exp(ygrid(i)), K, config.mu, config.sigma, config.T);
  }

  Eigen::VectorXd error = Eigen::VectorXd::Zero(price_steps);
  for (int i = 0; i < price_steps; ++i) {
    error(i) = std::abs(pi(i, time_steps) - analytical(i));
  }

  // Check that the L2 error is within an acceptable range
  EXPECT_LT(error.norm(), 0.01);
}