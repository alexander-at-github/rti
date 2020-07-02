#pragma once

#include <cmath>

#include "rti/ray/cos_hemi.hpp"
#include "rti/reflection/i_reflection_model.hpp"
#include "rti/rng/i_rng.hpp"

namespace rti { namespace reflection {
  template<typename Ty>
  class diffuse : public rti::reflection::i_reflection_model<Ty> {
  public:
    rti::util::pair<rti::util::triple<Ty> >
    use(RTCRay& pRayIn, RTCHit& pHitIn, rti::geo::i_abs_geometry<Ty>& pGeometry,
        rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {

      auto primID = pHitIn.primID;
      // Get an origin for the refelcted ray from the i_geometry implementation
      auto newOrigin = pGeometry.get_new_origin(pRayIn, primID);

      auto normal = pGeometry.get_normal(primID);
      /* Compute lambertian reflection with respect to surface normal */
      auto orthonormalBasis = rti::util::get_orthonormal_basis(normal);
      auto direction = rti::ray::cos_hemi::get<Ty>(orthonormalBasis, pRng, pRngState);

      return {newOrigin, direction};
    }
  };
}}
