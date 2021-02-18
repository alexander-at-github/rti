#pragma once

#include "i_source.hpp"
#include "../rng/i_rng.hpp"
#include "../rng/non_rng.hpp"

namespace rti {
namespace ray {
  template<typename numeric_type>
  class nonrng_loc_rng_dir_source : public i_source {
    
public:
    // interesting settings:
    // (A) 46080000 rays with mNumXOrgs == 144 and mNumYOrgs == 32 noise level: ~0.98929 (0.01071)
    //     with Chi-Squared == 251.797 and degrees of freedom == 303 (p-value of 0.985126)
    // (B) 299635200 rays with mNumXOrgs == 144 and mNumYOrgs == 32 noise level: ~0.99583 (0.00417)
    //     with Chi-Squared == 278.035 and degrees of freedom == 303 (p-value of 0.845575)

    using pair = util::pair<numeric_type>;
    
    nonrng_loc_rng_dir_source(numeric_type pZval, pair pC1, pair pC2, size_t pNumRays,
                  ray::i_direction<numeric_type>& pDirection) :
      mZval(pZval), mC1(pC1), mC2(pC2), mNumRays(pNumRays), mDirection(pDirection)
    {
      std::cout << "USING nonrng_loc_rng_dir_source" << std::endl;
      mNumXOrgs = 144;
      mNumYOrgs = 32;
      mNumXOrgs = 288;
      mNumYOrgs = 64;
      dx = (mC2[0] - mC1[0]) / mNumXOrgs;
      dy = (mC2[1] - mC1[1]) / mNumYOrgs;
      std::cout << "mNumRays == " << mNumRays << " mNumRays / mNumXOrgs / mNumYOrgs == " << (double) mNumRays / mNumXOrgs / mNumYOrgs << std::endl;
      assert(mNumRays % (mNumXOrgs * mNumYOrgs) == 0 && "Correctness Assumption");
      mNumRaysPerOrg = mNumRays / mNumXOrgs / mNumYOrgs;
    }
    
  ~nonrng_loc_rng_dir_source() {}

    void fill_ray(RTCRay& pRay, rng::i_rng& pRng,
                  rng::i_rng::i_state& pRngState1, rng::i_rng::i_state& pRngState2,
                  rng::i_rng::i_state& pRngState3, rng::i_rng::i_state& pRngState4
                  ) override final {
      // std::cerr << "HERE" << std::endl;
      if (mNumYOidx >= mNumYOrgs) {
        // indicates an error
        assert(false && "ERROR");
      }

      auto xx = mC1[0] + dx/2 + dx * mNumXOidx;
      auto yy = mC1[1] + dy/2 + dy * mNumYOidx;
      auto dir = mDirection.get(pRng, pRngState3, pRngState4);
      pRay.org_x = (numeric_type) xx;
      pRay.org_y = (numeric_type) yy;
      pRay.org_z = (numeric_type) mZval;
      pRay.dir_x = dir[0];
      pRay.dir_y = dir[1];
      pRay.dir_z = dir[2];
      if (pRay.dir_z == 0) {
        // pRay.dir_z can lead to infinitely bouncing rays.
        // We do a quick fix.
        pRay.dir_z = -1e-3;
      }
      // advance index
      mNumRaysPerOrgIdx += 1;
      if (mNumRaysPerOrgIdx >= mNumRaysPerOrg) {
        mNumRaysPerOrgIdx = 0;
        mNumXOidx += 1;
        if (mNumXOidx >= mNumXOrgs) {
          mNumXOidx = 0;
          mNumYOidx += 1;
          if (mNumYOidx >= mNumYOrgs) {
            // don't do anything but instead recognize error on next iteration.
            //mNumYOidx = 0;
          }
        }
      }
      // std::cerr << "HERE END" << std::endl;
    }
    
  private:

    numeric_type mZval;
    pair mC1;
    pair mC2;
    size_t mNumRays;
    ray::i_direction<numeric_type>& mDirection;
      
    size_t mNumXOrgs;
    size_t mNumYOrgs;

    size_t mNumXOidx = 0;
    size_t mNumYOidx = 0;
    double dx;
    double dy;

    size_t mNumRaysPerOrgIdx = 0;
    size_t mNumRaysPerOrg = 0;
};
} // namespace ray
} // namespace rti
