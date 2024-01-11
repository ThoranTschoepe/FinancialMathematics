#include <FinancialMathematics/FinancialMathematics.h>

#include <benchmark/benchmark.h>

using namespace FinancialMathematics;

namespace bm = benchmark;

static void BM_BasicMC(bm::State &state) {
  SetSeed(1);

  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  auto payoff = [](const std::vector<double> &x) {
    return std::max(x.back() - 110.0, 0.0);
  };

  EtaMatrix eta = EtaMatrix(1000000, 1);
  AssetPaths asset_paths = BrownianMotion(eta, config);

  for (auto _ : state) {
    bm::DoNotOptimize(basicMonteCarlo(t, payoff, asset_paths, config));
  }
}

static void BM_MC_Anti(bm::State &state) {
  SetSeed(1);

  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  auto payoff = [](const std::vector<double> &x) {
    return std::max(x.back() - 110.0, 0.0);
  };

  EtaMatrix eta = EtaMatrix(1000000, 1);
  AssetPaths asset_paths = BrownianMotion(eta, config);
  AssetPaths asset_paths_anti = BrownianMotion(eta, config);

  for (auto _ : state) {
    bm::DoNotOptimize(MonteCarloAntithetic(t, payoff, asset_paths, asset_paths_anti, config));
  }
}

static void BM_MC_ControlVariate(bm::State &state) {
  SetSeed(1);

  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  auto payoff = [](const std::vector<double> &x) {
    return std::max(x.back() - 110.0, 0.0);
  };

  uint32_t paths = 1000000;
  uint32_t paths_opt = paths / 10;

  EtaMatrix eta = EtaMatrix(paths, 1);
  EtaMatrix eta_opt = EtaMatrix(paths_opt, 1);
  AssetPaths asset_paths = BrownianMotion(eta, config);
  AssetPaths asset_paths_opt = BrownianMotion(eta_opt, config);

  for (auto _ : state) {
    bm::DoNotOptimize(MonteCarloControlVariateUnderlying(t, payoff, asset_paths, asset_paths_opt, config));
  }
}

static void BM_MC_ImportanceSamplingConstDrift(bm::State &state) {
  SetSeed(1);

  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  auto payoff = [](const std::vector<double> &x) {
    return std::max(x.back() - 220.0, 0.0);
  };

  EtaMatrix eta = EtaMatrix(1000000, 1);

  for (auto _ : state) {
    bm::DoNotOptimize(MonteCarloImportanceSampling(t, payoff, 5, 0, eta, config));
  }
}

static void BM_FiniteDifference(bm::State &state) {
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
  int time_steps = 1000;
  int price_steps = 200;
  constexpr bool quadratic_payoff = false;

  for (auto _ : state) {
    bm::DoNotOptimize(finiteDifferences<quadratic_payoff>(time_steps,
                                        price_steps,
                                        payoff,
                                        config.mu,
                                        config.sigma,
                                        config.T));
  }
}

static void BM_AmericanPut(bm::State &state) {
  SetSeed(1);

  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  uint32_t paths = 1000;
  uint32_t steps = 500;

  EtaMatrix eta = EtaMatrix(paths, steps);

  AssetPaths asset_paths = BrownianMotion(eta, config);

  for (auto _ : state) {
    bm::DoNotOptimize(American::american_put_martingale_MC(asset_paths, config));
  }
}

static void BM_AmericanPutTree(bm::State &state) {
  SetSeed(1);

  double t = 0.0;
  Config config;
  config.x0 = 100.0;
  config.mu = 0.05;
  config.sigma = 0.2;
  config.T = 1.0;

  uint32_t n = 1000;

  for (auto _ : state) {
    bm::DoNotOptimize(American::american_put_binomial_tree(config, 100, n));
  }
}


BENCHMARK(BM_BasicMC)->Unit(bm::kMillisecond);
/*BENCHMARK(BM_MC_Anti)->Unit(bm::kMillisecond);
BENCHMARK(BM_MC_ControlVariate)->Unit(bm::kMillisecond);
BENCHMARK(BM_MC_ImportanceSamplingConstDrift)->Unit(bm::kMillisecond);
BENCHMARK(BM_FiniteDifference)->Unit(bm::kMillisecond);*/

BENCHMARK(BM_AmericanPut)->Unit(bm::kMillisecond);
BENCHMARK(BM_AmericanPutTree)->Unit(bm::kMillisecond);

BENCHMARK_MAIN();