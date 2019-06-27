#pragma once

#include <cmath>

#include "rti/i_reflection_model.hpp"

namespace rti {
  class lambertian_reflection : public i_reflection_model {
  public:
    bool use(RTCRayHit& pRayhit, const i_geometry& pGeometry, i_hit_counter& pHitcounter) const override final {
      // TODO: Get random number and decide whether or not to reflect.
      // TODO: get surface normal at intersection
      // TODO: Compute lambertian reflection with respect to surface normal

      // TODO: FIX
      pHitcounter.use(pRayhit);
      return false;
    }

  private:
    // Returns some orthonormal basis containing a the input vector pVector.
    // Is deterministic, i.e., for one input it will return always the same
    // result.
    template<typename T>
    rti::triple<rti::triple<T> > get_orthonormal_basis(rti::triple<T> pVector) {
      rti::triple<rti::triple<T> > rr;
      rr.frst.frst = pVector.frst;
      rr.frst.scnd = pVector.scnd;
      rr.frst.thrd = pVector.thrd;

      // Calculate a vector (rr.scnd) which is perpendicular to rr.frst
      rr.scnd.frst = 1; // magic number; will be normalized later
      rr.scnd.scnd = 0; // magic number
      rr.scnd.thrd =
        ((rr.frst.frst * rr.scnd.frst) + (rr.frst.scnd * rr.scnd.scnd)) / rr.frst.thrd;

      // Calculat cross product of rr.frst and rr.scnd to form orthogonal basis
      rr.thrd.frst = rr.frst.scnd * rr.scnd.thrd - rr.frst.thrd * rr.scnd.scnd;
      rr.thrd.scnd = rr.frst.thrd * rr.scnd.frst - rr.frst.frst * rr.scnd.thrd;
      rr.thrd.thrd = rr.frst.frst * rr.scnd.scnd - rr.frst.scnd * rr.scnd.frst;

      // Normalize the length of these vectors.
      T frstNorm = std::sqrt(
          rr.frst.frst * rr.frst.frst +
          rr.frst.scnd * rr.frst.scnd +
          rr.frst.thrd * rr.frst.thrd);
      rr.frst.frst *= frstNorm;
      rr.frst.scnd *= frstNorm;
      rr.frst.thrd *= frstNorm;
      T scndNorm = std::sqrt(
          rr.scnd.frst * rr.scnd.frst +
          rr.scnd.scnd * rr.scnd.scnd +
          rr.scnd.thrd * rr.scnd.thrd);
      rr.scnd.frst *= scndNorm;
      rr.scnd.scnd *= scndNorm;
      rr.scnd.thrd *= scndNorm;
      T thrdNorm = std::sqrt(
          rr.thrd.frst * rr.thrd.frst +
          rr.thrd.scnd * rr.thrd.scnd +
          rr.thrd.thrd * rr.thrd.thrd);
      rr.thrd.frst *= thrdNorm;
      rr.thrd.scnd *= thrdNorm;
      rr.thrd.thrd *= thrdNorm;

      // Sanity check
      assert(std::abs(dot_product(rr.frst, rr.scnd)) < 1e-6 && "Error in orthonormal basis computation");
      assert(std::abs(dot_product(rr.frst, rr.thrd)) < 1e-6 && "Error in orthonormal basis computation");
      assert(std::abs(dot_product(rr.scnd, rr.thrd)) < 1e-6 && "Error in orthonormal basis computation");

      return rr;
    }

    template<typename T>
    T dot_product(rti::triple<T> pFrst, rti::triple pScnd) {
      return pFrst.frst * pScnd.frst + pFrst.scnd * pScnd.scnd + pFrst.thrd * pScnd.thrd;
    }
  };
} // namespace rti
