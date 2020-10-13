#pragma once

#include "../rng/i_rng.hpp"

namespace rti { namespace reflection {
  template<typename Ty>
  class i_reflection_model {
  public:
    // Pure Virtual Class
    virtual ~i_reflection_model() {}
    // Decides whether or not to reflect. If a reflection should happen, it sets
    // the origin and direction in the RTCRayHit object and returns true. If no
    // reflection should happen, then it does not change pRayhit and returns
    // false.
    virtual rti::util::pair<rti::util::triple<Ty> >
    use(RTCRay& rayin, RTCHit& hitin, rti::geo::i_abs_geometry<Ty>&,
        rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) = 0;
  };
}} // namespace

