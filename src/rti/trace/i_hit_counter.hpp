#pragma once

#include "../util/ostream_overload_template.hpp"

namespace rti { namespace trace {
  using rti::util::operator<<;
  class i_hit_counter {
    // Pure Virtual Class
  public:
    // A hit counter counts hits. That is, it also is not able to account for
    // fractions. We suggest that if one wants to accumulate rational numbers
    // then one should write a new interface (e.g., i_hit_accumulator).
    virtual ~i_hit_counter() {}
    virtual void use(const RTCRayHit& pRayhit) = 0;
    virtual std::vector<size_t> get_counts() = 0;
    virtual void print(std::ostream& pOs) const = 0;
  };
}} // namespace
