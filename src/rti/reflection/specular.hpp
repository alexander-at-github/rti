#pragma once

#include "i_reflection_model.hpp"

namespace rti { namespace reflection {
  template<typename Ty>
  class specular : public rti::reflection::i_reflection_model<Ty> {
  public:
    rti::util::pair<rti::util::triple<Ty> >
    use(RTCRay& pRayIn, RTCHit& pHitIn, rti::geo::meta_geometry<Ty>& pGeometry,
        rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {

      auto primID = pHitIn.primID;
      auto normal = pGeometry.get_normal(primID);
      // Instead of querying the geometry object for the surface normal one could used
      // the (unnormalized) surface normal provided by the rayhit.hit struct.
      auto dirOldInv = rti::util::inv ( rti::util::triple<Ty> {pRayIn.dir_x, pRayIn.dir_y, pRayIn.dir_z} );
      // For computing the specular refelction direction we need the vectors to be normalized.
      assert(rti::util::is_normalized(normal) && "surface normal vector is supposed to be normalized");
      assert(rti::util::is_normalized(dirOldInv) && "direction vector is supposed to be normalized");
      // Compute new direction
      auto direction =
        rti::util::diff(rti::util::scale(2 * rti::util::dot_product(normal, dirOldInv), normal), dirOldInv);

      // // instead of using this epsilon one could set tnear to a value other than zero
      // auto epsilon = 1e-6;
      // auto ox = pRayIn.org_x + pRayIn.dir_x * pRayIn.tfar + normal[0] * epsilon;
      // auto oy = pRayIn.org_y + pRayIn.dir_y * pRayIn.tfar + normal[1] * epsilon;
      // auto oz = pRayIn.org_z + pRayIn.dir_z * pRayIn.tfar + normal[2] * epsilon;
      // auto newOrigin = rti::util::triple<Ty> {(Ty) ox, (Ty) oy, (Ty) oz};
      auto newOrigin = pGeometry.get_new_origin(pRayIn, primID);

      return {newOrigin, direction};
    }

  private:
  };
}} // namespace
