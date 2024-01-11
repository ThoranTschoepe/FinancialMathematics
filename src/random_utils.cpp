#include <FinancialMathematics/random_utils.h>

#include <random>

namespace FinancialMathematics {
  std::mt19937_64 generator;
  std::normal_distribution<double> normal_distribution(0.0, 1.0);

  void SetSeed(int seed) {
    generator.seed(seed);
  }

  std::vector<std::vector<double>> GenerateSTDNormalRandomNumbers(uint32_t n, uint32_t m) {
    std::vector<std::vector<double>> eta(n, std::vector<double>(m));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; ++j) {
        eta[i][j] = normal_distribution(generator);
      }
    }
    return eta;
  }
}