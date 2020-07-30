#include <benchmark/benchmark.h>

#include <compare>
#include <vector>

inline long lsgnNumeric(long x) {return (x > 0) ? 1 : ((x < 0) ? -1 : 0);}

inline long lsgnBinary(long x) {return ((x >> 31) | (x >> 31)^1) & (~(!(!x)) + 1);}

inline std::strong_ordering lsgnSpace(long x) {return x <=> 0;}

inline long lsgnCmp(long x) {return (0 < x) - (0 > x);}

static std::vector<long> nn = {-14657815, 8137456};

static void CheckLsgnNumeric(benchmark::State& state)
{
  int idx = 0;
  for (auto _ : state) {
    auto result = lsgnNumeric(nn[idx]);
    idx = (idx + 1) % 1;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(CheckLsgnNumeric);

static void CheckLsgnBinary(benchmark::State& state)
{
  int idx = 0;
  for (auto _ : state) {
    auto result = lsgnBinary(nn[idx]);
    idx = (idx + 1) % 1;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(CheckLsgnBinary);

static void CheckLsgnSpace(benchmark::State& state)
{
  int idx = 0;
  for (auto _ : state) {
    auto result = lsgnSpace(nn[idx]);
    idx = (idx + 1) % 1;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(CheckLsgnSpace);

static void CheckLsgnCmp(benchmark::State& state)
{
  int idx = 0;
  for (auto _ : state) {
    auto result = lsgnCmp(nn[idx]);
    idx = (idx + 1) % 1;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(CheckLsgnCmp);
