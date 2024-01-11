#ifndef FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_ASSETPATHS_H_
#define FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_ASSETPATHS_H_

#include <vector>
#include <cstdint>
#include <iostream>
#include <functional>

#include <FinancialMathematics/random_utils.h>

namespace FinancialMathematics {
  class EtaMatrix {
    public:
    EtaMatrix(uint32_t paths, uint32_t steps) : eta_(GenerateSTDNormalRandomNumbers(paths, steps)) {}

    EtaMatrix(std::vector<std::vector<double>> eta) : eta_(std::move(eta)) {}

    [[nodiscard]] uint32_t paths() const { return eta_.size(); }
    [[nodiscard]] uint32_t steps() const { return eta_[0].size(); }

    std::vector<double> operator[](uint32_t i) const {
      return eta_[i];
    }

   private:
    std::vector<std::vector<double>> eta_;
  };

  struct Config{
    double x0;
    double mu;
    double sigma;
    double T;
  };

  class AssetPaths {
  public:
    std::vector<double> operator[](uint32_t i) const {
      return x_[i];
    }

    [[nodiscard]] uint32_t paths() const { return x_.size(); }
    [[nodiscard]] uint32_t steps() const { return x_[0].size(); }

    friend std::ostream& operator<<(std::ostream& os, const AssetPaths& asset_paths);

   private:
    explicit AssetPaths(std::vector<std::vector<double>> x);

  private:
    std::vector<std::vector<double>> x_;

    friend AssetPaths BrownianMotion(const EtaMatrix& eta, const Config& config, const std::function<double(double, uint32_t)>& eta_op);
    friend AssetPaths BrownianMotionAntitheticEta(const EtaMatrix& eta, const Config& config, const std::function<double(double, uint32_t)>& eta_op);
  };
}

#endif //FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_ASSETPATHS_H_
