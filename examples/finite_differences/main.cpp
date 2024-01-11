#include <FinancialMathematics/FinancialMathematics.h>

#include <iostream>

using namespace FinancialMathematics;

int main() {
  SetSeed(13);

  double t = 0.0;
  Config config;
  config.x0 = 1.0;
  config.mu = 0;
  config.sigma = 0.5;
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

  std::cout << "L2 error analytical: " << error.norm() << std::endl;

  return 0;
}