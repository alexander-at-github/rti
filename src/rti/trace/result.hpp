#pragma once

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
    std::string inputFilePath;
    size_t numRays;
    size_t hitc;
    size_t nonhitc;

    void print(std::ostream& pOs) const {
      pOs
        << "===[" << typeid(this).name() << "]=== "
        << geometryClassName << " "
        << inputFilePath << " "
        << numRays << "rays "
        // << hitc << "hits "
        // << nonhitc << "nonhits "
        << timeNanoseconds*1e-9 << "seconds"
        << std::endl;
    }
  };
}}
