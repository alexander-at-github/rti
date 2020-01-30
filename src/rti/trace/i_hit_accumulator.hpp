#pragma once

#include "rti/util/ostream_overload_template.hpp"

namespace rti { namespace trace {
  using rti::util::operator<<;
  template<typename Ty>
  class i_hit_accumulator {
  public:
    virtual ~i_hit_accumulator() {}
    virtual void use(unsigned int pPrimID, Ty value) = 0;
    virtual std::vector<Ty> get_values() = 0;
    virtual std::vector<size_t> get_cnts() = 0;
    virtual std::vector<Ty> get_relative_errors() = 0;
    virtual std::vector<Ty> get_vov() = 0;
    virtual void print(std::ostream& pOs) const = 0;
    virtual Ty get_relative_error_for_id(unsigned int pPrimID) = 0;
  };
}}
