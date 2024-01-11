#ifndef FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_BROWNIAN_MOTION_H_
#define FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_BROWNIAN_MOTION_H_

#include <FinancialMathematics/AssetPaths.h>

#include <vector>

namespace FinancialMathematics {
  AssetPaths BrownianMotion(const EtaMatrix& eta, const Config& config, const std::function<double(double, uint32_t)>& eta_op = [](double x, uint32_t i) { return x; });
  AssetPaths BrownianMotionAntitheticEta(const EtaMatrix& eta, const Config& config, const std::function<double(double, uint32_t)>& eta_op = [](double x, uint32_t i) { return x; });
}

#endif //FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_BROWNIAN_MOTION_H_
