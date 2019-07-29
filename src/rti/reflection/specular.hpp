#pragma once

#include "rti/reflection/i_reflection_model.hpp"

namespace rti { namespace reflection {
  template<typename Ty>
  class specular : public rti::reflection::i_reflection_model<Ty> {
  public:
    bool use(RTCRayHit& pRayhit,
             rti::rng::i_rng& pRng,
             rti::rng::i_rng::i_state& pRngState,
             rti::geo::i_abs_geometry<Ty> const& pGeometry,
             rti::trace::i_hit_counter& pHitcounter) const override final {

      auto primID = pRayhit.hit.primID;
      auto normal = pGeometry.get_normal(primID);
      // Instead of querying the geometry object for the surface normal one could used
      // the (unnormalized) surface normal provided by the rayhit.hit struct.
      auto dirOldInv = rti::util::inv ( rti::util::triple<Ty> {pRayhit.ray.dir_x, pRayhit.ray.dir_y, pRayhit.ray.dir_z} );
      // For computing the the specular refelction direction we need the vectors to be normalized.
      assert(rti::util::is_normalized(normal) && "surface normal vector is supposed to be normalized");
      assert(rti::util::is_normalized(dirOldInv) && "direction vector is supposed to be normalized");
      // Compute new direction
      auto direction = rti::util::diff(rti::util::scale(2 * rti::util::dot_product(normal, dirOldInv), normal), dirOldInv);

      auto epsilon = 1e-6;
      auto ox = pRayhit.ray.org_x + pRayhit.ray.dir_x * pRayhit.ray.tfar + normal[0] * epsilon;
      auto oy = pRayhit.ray.org_y + pRayhit.ray.dir_y * pRayhit.ray.tfar + normal[1] * epsilon;
      auto oz = pRayhit.ray.org_z + pRayhit.ray.dir_z * pRayhit.ray.tfar + normal[2] * epsilon;
      auto newOrigin = rti::util::triple<Ty> {ox, oy, oz};

      pRayhit.ray.org_x = newOrigin[0];
      pRayhit.ray.org_y = newOrigin[1];
      pRayhit.ray.org_z = newOrigin[2];

      pRayhit.ray.dir_x = direction[0];
      pRayhit.ray.dir_y = direction[1];
      pRayhit.ray.dir_z = direction[2];

      return true;
    }
  };
}} // namespace
