#pragma once

#include <cmath>

#include "rti/cos_hemi.hpp"
#include "rti/cstdlib_rng.hpp"
#include "rti/i_reflection_model.hpp"
#include "rti/i_rng.hpp"

namespace rti {
  template<typename Ty>
  class diffuse_reflection : public i_reflection_model<Ty> {
  public:
    diffuse_reflection(Ty pStickingC) :
      mStickingC(pStickingC) {
    }

    bool use(
             RTCRayHit& pRayhit,
             rti::i_rng& pRng,
             rti::i_rng::i_state& pRngState,
             const i_geometry<Ty>& pGeometry,
             i_hit_counter& pHitcounter) const override final {

      Ty epsilon = 1e-6;
      // THIS CODE IS HAND-CRAFTED FOR THE CYLINDER WITH SOURCE PLANE AT X == 0.
      if (pRayhit.ray.org_x + pRayhit.ray.dir_x * pRayhit.ray.tfar <= epsilon) {
        // Don't reflect and don't count
        return false;
      }

      /* Get random number and decide whether or not to reflect. */
      uint64_t rndm = pRng.get(pRngState);
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
      // //rti::triple<Ty> pGeometry.get_normal(pRayhit.hit.geomID);
      // rti::triple<Ty> normalO {pRayhit.hit.Ng_x, pRayhit.hit.Ng_y, pRayhit.hit.Ng_z};

      // TODO: Where to start the new ray?

      Ty xxahit = pRayhit.ray.org_x + pRayhit.ray.dir_x * pRayhit.ray.tfar;
      Ty yyahit = pRayhit.ray.org_y + pRayhit.ray.dir_y * pRayhit.ray.tfar;
      Ty zzahit = pRayhit.ray.org_z + pRayhit.ray.dir_z * pRayhit.ray.tfar;
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
      // // auto diff = rti::diff(normalO, pGeometry.get_normal(pRayhit.hit.primID));
      // // assert (diff[0] < 1e-5 && diff[1] < 1e-5 && diff[2] < 1e-5 && "Caluclation of surface normal");
      // // // Is that the normal at the hit location? Yes.

      // // CAUTION: For now we use inverted geometries from GMSH, that is, we have
      // // to invert the normal.
      // auto normal = rti::inv(normalO);

      auto normal = pGeometry.get_normal(pRayhit.hit.primID);


      /* Compute lambertian reflection with respect to surface normal */
      auto orthonormalBasis = get_orthonormal_basis(normal);

      // Set new origin and direction in pRayhit
      // We add a small epsilon to the origin to make sure that we do not intersect the same surface again.
      // Without that the simulation does not work.
      Ty epsilonOrg = 1e-6;
      auto obFirst = orthonormalBasis[0]; // normalized surface normal
      pRayhit.ray.org_x = xxahit + obFirst[0] * epsilonOrg;
      pRayhit.ray.org_y = yyahit + obFirst[1] * epsilonOrg;
      pRayhit.ray.org_z = zzahit + obFirst[2] * epsilonOrg;
      auto direction = rti::cos_hemi::get<Ty>(orthonormalBasis, pRng, pRngState);
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
    rti::triple<rti::triple<Ty> > get_orthonormal_basis(const rti::triple<Ty> pVector) const {
      rti::triple<rti::triple<Ty> > rr;
      rr[0] = pVector;
      // rr[0][0] = pVector[0];
      // rr[0][1] = pVector[1];
      // rr[0][2] = pVector[2];

      // Calculate a vector (rr[1]) which is perpendicular to rr[0]
      // https://math.stackexchange.com/questions/137362/how-to-find-perpendicular-vector-to-another-vector#answer-211195
      rti::triple<Ty> candidate0 {rr[0][2], rr[0][2], -(rr[0][0] + rr[0][1])};
      rti::triple<Ty> candidate1 {rr[0][1], -(rr[0][0] + rr[0][2]), rr[0][1]};
      rti::triple<Ty> candidate2 {-(rr[0][1] + rr[0][2]), rr[0][0], rr[0][0]};
      // We choose the candidate which maximizes the sum of its components, because we want to avoid
      // numeric errors and that the result is (0, 0, 0).
      std::array<rti::triple<Ty>, 3> cc = {candidate0, candidate1, candidate2};
      auto sumFun = [](rti::triple<Ty> oo){return oo[0] + oo[1] + oo[2];};
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
      rr[2] = rti::cross_product(rr[0], rr[1]);

      // Normalize the length of these vectors.
      rti::normalize(rr[0]);
      rti::normalize(rr[1]);
      rti::normalize(rr[2]);

      // Sanity check
      Ty epsilon = 1e-6;
      assert(std::abs(rti::dot_product(rr[0], rr[1])) < epsilon && "Error in orthonormal basis computation");
      assert(std::abs(rti::dot_product(rr[1], rr[2])) < epsilon && "Error in orthonormal basis computation");
      assert(std::abs(rti::dot_product(rr[2], rr[0])) < epsilon && "Error in orthonormal basis computation");

      return rr;
    }

  };
} // namespace rti
