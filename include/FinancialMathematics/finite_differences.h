#ifndef FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_FINITE_DIFFERENCES_H_
#define FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_FINITE_DIFFERENCES_H_

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace FinancialMathematics {

  template<bool quadratic_payoff>
  Eigen::MatrixXd finiteDifferences(int time_steps,
                                    int price_steps,
                                    const std::function<double(double)> &payoff,
                                    double r,
                                    double sigma,
                                    double T) {
    Eigen::setNbThreads(6);

    double xmax = 10000;
    double xmin = 0.0001;

    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(price_steps, time_steps + 1);

    double a = 0.5 * (2 * r / (sigma * sigma) - 1);
    double b = 0.25 * std::pow((2 * r / (sigma * sigma) + 1), 2);
    double Tau = 0.5 * sigma * sigma * T;
    double ymax = std::log(xmax);
    double ymin = std::log(xmin);
    double delta_tau = Tau / time_steps;
    double delta_y = (ymax - ymin) / (price_steps - 1);

    if (delta_tau > 0.5 * delta_y * delta_y) {
      std::cout << "Error: delta_tau is too large" << std::endl;
      std::cout << "delta_tau = " << delta_tau << std::endl;
      std::cout << "0.5*delta_y^2 = " << 0.5 * delta_y * delta_y << std::endl;
    }

    Eigen::VectorXd ygrid = Eigen::VectorXd::LinSpaced(price_steps, ymin, ymax);

    for (int i = 1; i < price_steps; ++i) {
      u(i, 0) = std::exp(a * ygrid(i)) * payoff(std::exp(ygrid(i)));
    }

    double precalc1 = b - r * 2 / (sigma * sigma);
    double precalc2 = b + r * 2 / (sigma * sigma) + 2;
    double precalc3 = a * ymin;
    double precalc4 = a * ymax;
    double precalc5 = payoff(0);
    double precalc6 = payoff(std::exp(ymax));

    Eigen::VectorXd taugrid = Eigen::VectorXd::LinSpaced(time_steps + 1, 0, Tau);

    for (int n = 0; n <= time_steps; ++n) {
      if constexpr (quadratic_payoff) {
        u(0, n) = std::exp(precalc1 * taugrid[n] + precalc3) * precalc5;
        u(price_steps - 1, n) = std::exp(precalc4 + precalc2 * taugrid[n]) * precalc6;
      } else {
        u(0, n) = std::exp(precalc1 * taugrid[n] + precalc3) * precalc5;
        u(price_steps - 1, n) = std::exp(precalc4 + b * taugrid[n]) * precalc6;
      }
    }

    double s = delta_tau / (delta_y * delta_y);
    Eigen::VectorXd offdiag = Eigen::VectorXd::Constant(price_steps - 1, s);
    Eigen::VectorXd diag = Eigen::VectorXd::Constant(price_steps, 1 - 2 * s);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(price_steps, price_steps);
    A.diagonal() = diag;
    A.diagonal(-1) = offdiag;
    A.diagonal(1) = offdiag;

    for (int n = 1; n <= time_steps; ++n) {
      u.block(1, n, u.rows() - 2, 1).noalias() = (A * u.col(n - 1)).block(1, n, u.rows() - 2, 1);
    }

    Eigen::MatrixXd meshtau = Eigen::MatrixXd::Zero(price_steps, time_steps + 1);
    Eigen::MatrixXd meshy = Eigen::MatrixXd::Zero(price_steps, time_steps + 1);

    for (int i = 0; i < price_steps; ++i) {
      meshtau.row(i).noalias() = taugrid.transpose();
    }
    for (int i = 0; i < time_steps + 1; ++i) {
      meshy.col(i) = ygrid;
    }

    Eigen::MatrixXd pi = -a * meshy - b * meshtau;
    pi = pi.array().exp() * u.array();

    return pi;
  }
}

#endif //FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_FINITE_DIFFERENCES_H_
