#pragma once

#include <chrono>

namespace rti {
  class test_result {
  public:
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point endTime;
    boost::int_least64_t timeNanoseconds = 0;
    std::string geometryClassName;
    std::string inputFilePath;
    size_t numRays;

    std::string to_string() {
      std::stringstream strstream;
      strstream << "(:class test_result " << inputFilePath << " "
                << geometryClassName << " "
                << numRays << "rays "
                << timeNanoseconds << "ns"
                << ")";
      //          << " startTime " << std::chrono::duration_cast<std::chrono::milliseconds>(startTime);
      //                << " endTime " << endTime ")";
      return strstream.str();
    }
  };
} // namespace rti
