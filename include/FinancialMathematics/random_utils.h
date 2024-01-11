#ifndef FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_RANDOM_UTILS_H_
#define FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_RANDOM_UTILS_H_

#include <vector>
#include <cstdint>

namespace FinancialMathematics {
  void SetSeed(int seed);

  std::vector<std::vector<double>> GenerateSTDNormalRandomNumbers(uint32_t n, uint32_t m);
}

#endif //FINANCIALMATHEMATICS_INCLUDE_FINANCIALMATHEMATICS_RANDOM_UTILS_H_
