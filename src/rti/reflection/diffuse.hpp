#pragma once

#include <cmath>

#include "rti/ray/cos_hemi.hpp"
#include "rti/reflection/i_reflection_model.hpp"
#include "rti/rng/cstdlib_rng.hpp"
#include "rti/rng/i_rng.hpp"

namespace rti { namespace reflection {
  template<typename Ty>
  class diffuse : public rti::reflection::i_reflection_model<Ty> {
  public:
    diffuse(Ty pStickingC) :
      mStickingC(pStickingC) {
    }

    bool use(RTCRayHit& pRayhit,
             rti::rng::i_rng& pRng,
             rti::rng::i_rng::i_state& pRngState,
             const rti::geo::i_abs_geometry<Ty>& pGeometry, // covariant parameter; does that work?
             rti::trace::i_hit_counter& pHitcounter) const override final {

      Ty epsilon = 1e-6;
      // THIS CODE IS HAND-CRAFTED FOR THE CYLINDER WITH SOURCE PLANE AT X == 0.

      /* Get random number and decide whether or not to reflect. */
      auto rndm = pRng.get(pRngState);
      if (rndm <= (pRng.max() * mStickingC)) {
        // Do not reflect
        pHitcounter.use(pRayhit);
        return false;
      }

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

      auto primID = pRayhit.hit.primID;
      // Get an origin for the refelcted ray from the i_geometry implementation
      auto newOrigin = pGeometry.get_new_origin(pRayhit, primID);


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
      auto orthonormalBasis = get_orthonormal_basis(normal);
      auto direction = rti::ray::cos_hemi::get<Ty>(orthonormalBasis, pRng, pRngState);

      pRayhit.ray.org_x = newOrigin[0];
      pRayhit.ray.org_y = newOrigin[1];
      pRayhit.ray.org_z = newOrigin[2];

      pRayhit.ray.dir_x = direction[0];
      pRayhit.ray.dir_y = direction[1];
      pRayhit.ray.dir_z = direction[2];

      return true;
    }

  private:
    // The sticking coefficient
    Ty mStickingC;

    // Returns some orthonormal basis containing a the input vector pVector
    // (possibly scaled) as the first element of the return value.
    // This function is deterministic, i.e., for one input it will return always
    // the same result.
    rti::util::triple<rti::util::triple<Ty> > get_orthonormal_basis(const rti::util::triple<Ty> pVector) const {
      rti::util::triple<rti::util::triple<Ty> > rr;
      rr[0] = pVector;
      // rr[0][0] = pVector[0];
      // rr[0][1] = pVector[1];
      // rr[0][2] = pVector[2];

      // Calculate a vector (rr[1]) which is perpendicular to rr[0]
      // https://math.stackexchange.com/questions/137362/how-to-find-perpendicular-vector-to-another-vector#answer-211195
      rti::util::triple<Ty> candidate0 {rr[0][2], rr[0][2], -(rr[0][0] + rr[0][1])};
      rti::util::triple<Ty> candidate1 {rr[0][1], -(rr[0][0] + rr[0][2]), rr[0][1]};
      rti::util::triple<Ty> candidate2 {-(rr[0][1] + rr[0][2]), rr[0][0], rr[0][0]};
      // We choose the candidate which maximizes the sum of its components, because we want to avoid
      // numeric errors and that the result is (0, 0, 0).
      std::array<rti::util::triple<Ty>, 3> cc = {candidate0, candidate1, candidate2};
      auto sumFun = [](rti::util::triple<Ty> oo){return oo[0] + oo[1] + oo[2];};
      int maxIdx = 0;
      for (size_t idx = 1; idx < cc.size(); ++idx) {
        if (sumFun(cc[idx]) > sumFun(cc[maxIdx])) {
          maxIdx = idx;
        }
      }
      assert (maxIdx < 3 && "Error in computation of perpenticular vector");
      rr[1] = cc[maxIdx];

      // Calculat cross product of rr[0] and rr[1] to form orthogonal basis
      // rr[2][0] = rr[0][1] * rr[1][2] - rr[0][2] * rr[1][1];
      // rr[2][1] = rr[0][2] * rr[1][0] - rr[0][0] * rr[1][2];
      // rr[2][2] = rr[0][0] * rr[1][1] - rr[0][1] * rr[1][0];
      rr[2] = rti::util::cross_product(rr[0], rr[1]);

      // Normalize the length of these vectors.
      rti::util::normalize(rr[0]);
      rti::util::normalize(rr[1]);
      rti::util::normalize(rr[2]);

      // Sanity check
      Ty epsilon = 1e-6;
      assert(std::abs(rti::util::dot_product(rr[0], rr[1])) < epsilon && "Error in orthonormal basis computation");
      assert(std::abs(rti::util::dot_product(rr[1], rr[2])) < epsilon && "Error in orthonormal basis computation");
      assert(std::abs(rti::util::dot_product(rr[2], rr[0])) < epsilon && "Error in orthonormal basis computation");

      return rr;
    }

  };
}} // namespace
