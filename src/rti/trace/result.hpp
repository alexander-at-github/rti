#pragma once

//#include <chrono>

//#include "rti/trace/i_hit_counter.hpp"
#include "rti/trace/i_hit_accumulator.hpp"
// include ostream overload template to provide out stream functionality
// by means of the print function.
#include "rti/util/ostream_overload_template.hpp"

namespace rti { namespace trace {
  // import name from ostream_overload_template.hpp into local namespace
  using rti::util::operator<<;
  template<typename Ty>
  class result {
  public:
    std::unique_ptr<rti::trace::i_hit_accumulator<Ty> > hitAccumulator;
    uint64_t timeNanoseconds = 0;
    std::string geometryClassName;
    std::string tracerFunctionName;
    std::string inputFilePath;
    size_t numRays;
    size_t hitc = 0;
    size_t nonhitc = 0;

    void print(std::ostream& pOs) const {
      pOs
        << "===[" << boost::core::demangle(typeid(this).name()) << "]=== "
        << geometryClassName << " "
        << inputFilePath << " "
        << numRays << "rays "
        // << hitc << "hits "
        // << nonhitc << "nonhits "
        << timeNanoseconds*1e-9 << "seconds"
        << std::endl;
    }
  };
}} // namespace
