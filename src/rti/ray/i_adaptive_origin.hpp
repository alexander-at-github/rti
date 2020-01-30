#pragma once

namespace rti { namespace ray {
  template<typename Ty>
  class i_adaptive_origin : public rti::ray::i_origin<Ty> {
  public:
    virtual void consider(rti::util::triple<Ty> xyz, double relativeerror) = 0;
    virtual void update_adaptive_sampling_state() = 0;
  };
}}
