#include <iostream>

#include <benchmark/benchmark.h>

#include "rti/ray/rectangle_origin_z.hpp"
#include "rti/rng/cstdlib_rng.hpp"
#include "rti/rng/mt64_rng.hpp"

using namespace rti;
using nt = float;

  auto zval = (nt) 2.0;
  // Corners of domain for benchmarking
  auto c1 = util::pair<nt> {2.0, 2.0}; // x, y
  auto c2 = util::pair<nt> {12.0, 24.0};

void rectangle_origin_z_mt64_rng(benchmark::State& pState)
{
  // auto zval = (nt) 2.0;
  // // Corners of domain for benchmarking
  // auto c1 = util::pair<nt> {2.0, 2.0}; // x, y
  // auto c2 = util::pair<nt> {12.0, 24.0};
  auto org = ray::rectangle_origin_z<nt> {zval, c1, c2};
  auto rng = rng::mt64_rng {};
  auto rngstate1 = rng::mt64_rng::state {1234567890};
  auto rngstate2 = rng::mt64_rng::state { 987654321};
  for (auto _ : pState) {
    auto sample = org.get(rng, rngstate1, rngstate2);
    benchmark::DoNotOptimize(sample);
  }
}

void std_uniform_real_distribution_z_mt64(benchmark::State& pState)
{
  auto zval = (nt) 2.0;
  auto rng1 = std::mt19937_64 {std::mt19937_64::default_seed};
  auto dist1 = std::uniform_real_distribution<nt> {c1[0], c2[0]}; // xmin, xmax
  auto rng2 = std::mt19937_64 {std::mt19937_64::default_seed / 2};
  auto dist2 = std::uniform_real_distribution<nt> {c1[1], c2[1]}; // ymin, ymax
  for (auto _ : pState) {
      auto xx = dist1(rng1);
      auto yy = dist2(rng2);
      auto sample = util::triple<nt> {xx, yy, zval};
      benchmark::DoNotOptimize(sample);
  }
}

BENCHMARK(rectangle_origin_z_mt64_rng)->UseRealTime();
BENCHMARK(std_uniform_real_distribution_z_mt64)->UseRealTime();
