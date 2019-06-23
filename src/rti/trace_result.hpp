#pragma once

#include <chrono>
#include <ostream>

// include ostream overload template to provide out stream functionality
// by means of the print function.
#include "rti/ostream_overload_template.hpp"

namespace rti {
  class trace_result {
  public:
    boost::int_least64_t timeNanoseconds = 0;
    std::string geometryClassName;
    std::string inputFilePath;
    size_t numRays;
    size_t hitc;
    size_t nonhitc;

    void print(std::ostream& pOs) const {
      pOs
        << geometryClassName << " "
        << inputFilePath << " "
        << numRays << "rays "
        << hitc << "hits "
        << nonhitc << "nonhits" << std::endl;
    }
  };
} // namespace rti
