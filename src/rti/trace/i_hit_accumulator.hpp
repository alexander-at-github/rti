#pragma once

#include "rti/util/ostream_overload_template.hpp"

namespace rti { namespace trace {
  using rti::util::operator<<;
  template<typename Ty>
  class i_hit_accumulator {
  public:
    virtual ~i_hit_accumulator() {}
    virtual void use(unsigned int pPrimID, Ty value) = 0;
    virtual std::vector<Ty> get_counts() = 0;
    virtual void print(std::ostream& pOs) const = 0;
  };
}}
