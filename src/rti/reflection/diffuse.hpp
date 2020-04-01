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

      //Ty epsilon = 1e-6;

      /* Reflect */
      //
      /* Get surface normal at intersection */
      //
      // // With the new surface representations we cannot use the surface normal
      // // provided by Embree anymore.
      // //
      // //
      // //rti::util::triple<Ty> pGeometry.get_normal(pRayhit.hit.geomID);
      // rti::util::triple<Ty> normalO {pRayhit.hit.Ng_x, pRayhit.hit.Ng_y, pRayhit.hit.Ng_z};

      auto primID = pHitIn.primID;
      // Get an origin for the refelcted ray from the i_geometry implementation
      auto newOrigin = pGeometry.get_new_origin(pRayIn, primID);


      // Ty xxahit = pRayhit.ray.org_x + pRayhit.ray.dir_x * pRayhit.ray.tfar;
      // Ty yyahit = pRayhit.ray.org_y + pRayhit.ray.dir_y * pRayhit.ray.tfar;
      // Ty zzahit = pRayhit.ray.org_z + pRayhit.ray.dir_z * pRayhit.ray.tfar;

      // { // DEBUG
      //   RLOG_DEBUG << std::endl;
      //   RLOG_DEBUG << "hit location: " << xxahit << " " << yyahit << " " << zzahit << ")" << std::endl;
      //   RLOG_DEBUG << "normal from rayhit.hit: ";
      //   for (size_t idx = 0; idx < normalO.size(); ++idx)
      //     RLOG_DEBUG << normalO.at(idx) << " ";
      //   RLOG_DEBUG << std::endl;
      //   auto geoNormal = pGeometry.get_normal(pRayhit.hit.primID);
      //   RLOG_DEBUG << "normal from geometry object: ";
      //   for (size_t idx = 0; idx < geoNormal.size(); ++idx)
      //     RLOG_DEBUG << geoNormal.at(idx) << " ";
      //   RLOG_DEBUG << std::endl;
      // } // DEBUG
      // // Sanity check; DOES NOT WORK ANYMORE WITH NEW SURFACE REPRESENTATION
      // // auto diff = rti::util::diff(normalO, pGeometry.get_normal(pRayhit.hit.primID));
      // // assert (diff[0] < 1e-5 && diff[1] < 1e-5 && diff[2] < 1e-5 && "Caluclation of surface normal");
      // // // Is that the normal at the hit location? Yes.

      // // CAUTION: For now we use inverted geometries from GMSH, that is, we have
      // // to invert the normal.
      // auto normal = rti::util::inv(normalO);

      auto normal = pGeometry.get_normal(primID);
      /* Compute lambertian reflection with respect to surface normal */
      auto orthonormalBasis = rti::util::get_orthonormal_basis(normal);
      auto direction = rti::ray::cos_hemi::get<Ty>(orthonormalBasis, pRng, pRngState);

      return {newOrigin, direction};
    }
  };
}} // namespace
