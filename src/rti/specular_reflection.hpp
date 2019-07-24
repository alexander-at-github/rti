#pragma once

#include "rti/i_reflection_model.hpp"

namespace rti {
  template<typename Ty>
  class specular_reflection : public i_reflection_model<Ty> {
  public:
    bool use(RTCRayHit& pRayhit, rti::i_rng& pRng, rti::i_rng::i_state& pRngState, const i_geometry<Ty>& pGeometry, i_hit_counter& pHitcounter) const override final {
      assert(false && "Not implemented");
      return false;
    }
  };
} // namespace rti
