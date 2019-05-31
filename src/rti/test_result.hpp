#pragma once

#include <chrono>

namespace rti {
  class test_result {
  public:
    // std::chrono::high_resolution_clock::time_point startTime;
    // std::chrono::high_resolution_clock::time_point endTime;
    boost::int_least64_t timeNanoseconds = 0;
    std::string geometryClassName;
    std::string inputFilePath;
    size_t numRays;
    size_t hitc;
    size_t nonhitc;

    std::string to_string() {
      std::stringstream strstream;
      strstream << "(:class test_result " << inputFilePath << " "
                << geometryClassName << " "
                << hitc << "hits "
                << nonhitc << "non-hits "
                << numRays << "rays "
                << timeNanoseconds << "ns"
                << ")";
      return strstream.str();
    }
  };
} // namespace rti
