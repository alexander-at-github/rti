#pragma once

#include "rti/util/ostream_overload_template.hpp"

namespace rti { namespace trace {
  using rti::util::operator<<;
  template<typename numeric_type>
  class i_hit_accumulator {
  protected:
    using internal_numeric_type = double;
  public:
    virtual ~i_hit_accumulator() {}
    virtual void use(unsigned int pPrimID, numeric_type value) = 0;
    virtual std::vector<internal_numeric_type> get_energy_values() = 0;
    virtual std::vector<size_t> get_cnts() = 0;
    virtual size_t get_cnts_sum() = 0;
    virtual std::vector<internal_numeric_type> get_relative_error() = 0;
    virtual std::vector<internal_numeric_type> get_vov() = 0;
    virtual void set_exposed_areas(std::vector<internal_numeric_type>&) = 0;
    virtual std::vector<internal_numeric_type> get_exposed_areas() = 0;
    virtual void print(std::ostream& pOs) const = 0;
    virtual numeric_type get_total_energy_used() = 0;
    virtual size_t get_size() = 0;
  };
}}
