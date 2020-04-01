#pragma once

#include"rti/ray/i_direction.hpp"

namespace rti { namespace ray {
  template<typename numeric_type>
  class constant_direction : public rti::ray::i_direction<numeric_type> {

  public:
    constant_direction(rti::util::triple<numeric_type> direction) :
      direction(direction) {
      rti::util::normalize(direction);
    }

    rti::util::triple<numeric_type>
    get(rti::rng::i_rng&, rti::rng::i_rng::i_state&) const
    {
      return direction;
    }

  private:
    rti::util::triple<numeric_type> direction;
  };
}}
