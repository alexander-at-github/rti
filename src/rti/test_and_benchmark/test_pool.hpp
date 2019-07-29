#pragma once

#include "rti/test_and_benchmark/test_run.hpp"

namespace rti { namespace test_and_benchmark {
  class test_pool {
  public:
    test_pool() :
      mTestRuns() {}
    test_pool(std::vector<rti::test_and_benchmark::test_run> pTestRuns) :
      mTestRuns(pTestRuns) {}
    void add_test_run(rti::test_and_benchmark::test_run pTestRun) {
      mTestRuns.push_back(pTestRun);
    }
    std::vector<rti::test_and_benchmark::test_result> run() {
      for (auto & test : mTestRuns) {
        rti::test_and_benchmark::test_result rr = test.run();
        mTestResults.push_back(rr);
      }
      return mTestResults;
    }
  private:
    std::vector<rti::test_and_benchmark::test_run> mTestRuns;
    std::vector<rti::test_and_benchmark::test_result> mTestResults;
  };
}} // namespace
