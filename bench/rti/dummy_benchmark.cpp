#include <benchmark/benchmark.h>

static void dummy_b_fun(benchmark::State& pState) {
  for (auto _ : pState)
    auto str = std::string {"foo"};
}
BENCHMARK(dummy_b_fun);
