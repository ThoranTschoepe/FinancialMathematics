#ifndef FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_MC_H_
#define FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_MC_H_

#include <vector>
#include <functional>
#include <FinancialMathematics/AssetPaths.h>

namespace FinancialMathematics {
  double basicMonteCarlo(double t,
                         const std::function<double(const std::vector<double> &)> &payoff,
                         const AssetPaths &asset_paths, const Config &config);

  double MonteCarloAntithetic(double t, const std::function<double(const std::vector<double> &)> &payoff,
                              const AssetPaths &ap, const AssetPaths &ap_antithetic, const Config &config);

  double MonteCarloControlVariateUnderlying(double t, const std::function<double(const std::vector<double> &)> &payoff,
                                            const AssetPaths &ap, const AssetPaths &ap_optimize, const Config &config);

  double MonteCarloImportanceSampling(double t, const std::function<double(const std::vector<double> &)> &payoff,
                                      double const_drift, double m_drift,
                                      const EtaMatrix &eta, const Config &config);

  namespace American {
    double golden_search(double a, double b, double tol, const std::function<double(double)> &f);
    double american_put_martingale_MC(const AssetPaths &asset_paths, const Config &config);

    double american_put_binomial_tree(const Config &config, double K, uint32_t n);
  }
}

#endif //FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_MC_H_
