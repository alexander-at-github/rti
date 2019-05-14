#pragma once

#include "rti/test_run.hpp"

namespace rti {
  class test_pool {
  public:
    test_pool() :
      mTestRuns() {}
    test_pool(std::vector<rti::test_run> pTestRuns) :
      mTestRuns(pTestRuns) {}
    void add_test_run(rti::test_run pTestRun) {
      mTestRuns.push_back(pTestRun);
    }
    std::vector<test_result> run() {
      for (auto & test : mTestRuns) {
        test_result rr = test.run();
        mTestResults.push_back(rr);
      }
      return mTestResults;
    }
  private:
    std::vector<test_run> mTestRuns;
    std::vector<test_result> mTestResults;
  };
} // namespace rti
