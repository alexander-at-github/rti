#pragma once

#include "i_source.hpp"
#include "../rng/i_rng.hpp"
#include "../rng/non_rng.hpp"

namespace rti {
namespace ray {
  template<typename numeric_type>
  class non_mc_source : public i_source {
    
public:
    // interesting settings:
    // mNonRng1Max, mNumXOrgs, mNumYOrgs
    // (A) 99u, 80, 4
    // (B) 99u, 80, 8
    //
    // (C) 99u, 144, 32 seems to roduce a result with low noise levels: ~0.99873 (0.00127)
    //     with Chi-Squared == 7.73513 and degrees of freedom == 303
    //     that is 46080000 rays
    // (D) 255u, 144, 32 seems to roduce a result with low noise levels: ~0.0.99926 (0.0007)
    //     with Chi-Squared == 17.9828 and degrees of freedom == 303
    //     that is 299635200 rays
    // (E) 999u, 144, 32 noise levels: ~0.9992 (0.0008)
    //     with Chi-Squared == 263.218 and degrees of freedom == 303
    //     that is 4608000000 rays

    using pair = util::pair<numeric_type>;
    
    non_mc_source(numeric_type pZval, pair pC1, pair pC2, size_t pNumRays,
                  ray::i_direction<numeric_type>& pDirection) :
      mZval(pZval), mC1(pC1), mC2(pC2), mNumRays(pNumRays), mDirection(pDirection),
      mNonRng1Max(999u),
      //mNonRng2Max(10),
      mNonRng1(mNonRng1Max)
      //mNonRng2(mNonRng2Max)
    {
      mNumXOrgs = 144;
      mNumYOrgs = 32;
      dx = (mC2[0] - mC1[0]) / mNumXOrgs;
      dy = (mC2[1] - mC1[1]) / mNumYOrgs;
      std::cout << "mNumRays == " << mNumRays << " non-mc-count == " << mNumXOrgs * mNumYOrgs * (mNonRng1Max+1) * (mNonRng1Max+1) << std::endl;
      // The non-rngs iterate from 0 to its max value (including).
      assert(mNumRays == mNumXOrgs * mNumYOrgs * (mNonRng1Max+1) * (mNonRng1Max+1) && "Correctness Assumption");
    }
    
  ~non_mc_source() {}


    // class non_rng : public rng::i_rng {
    //   struct state : public i_state {
    //     std::unique_ptr<rti::rng::i_rng::i_state> clone() const override final {
    //       assert(false);
    //       return std::make_unique<state> ();
    //     }
    //   };
    //   non_rng(size_t pMax) : mMax(pMax) {}
    //   uint64_t get(rti::rng::i_rng::i_state& pState) const override final {
    //     if (mCurr > mMax) {
    //       std::cerr << "ERROR non_rng" << std::endl;
    //       assert(false);
    //     }
    //     return mCurr++;
    //   }
    //   uint64_t min() const override final {
    //     return 0;
    //   }

    //   uint64_t max() const override final {
    //     return mMax;
    //   }
    // private:
    //   size_t mMax;
    //   size_t mCurr;
    // };

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
      auto dir = mDirection.get(mNonRng1, mNonRngState1, mNonRngState2);
      pRay.org_x = (numeric_type) xx;
      pRay.org_y = (numeric_type) yy;
      pRay.org_z = (numeric_type) mZval;
      pRay.dir_x = dir[0];
      pRay.dir_y = dir[1];
      pRay.dir_z = dir[2];
      if (pRay.dir_z == 0) {
        // This condition can lead to infinitely bouncing rays.
        // We do a quick fix.
        pRay.dir_z = -1e-3;
      }
      // advance indices
      // fun fact: the code is ugly
      mNonRngState1.mCurr += 1;
      if (mNonRngState1.mCurr > mNonRng1.mMax) {
        mNonRngState1.mCurr = 0;
        mNonRngState2.mCurr += 1;
        if (mNonRngState2.mCurr > mNonRng1.mMax) {
          mNonRngState2.mCurr = 0;
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
    size_t mNonRng1Max;
    //size_t mNonRng2Max;
    rng::non_rng mNonRng1;
    //rng::non_rng mNonRng2;
    rng::non_rng::state mNonRngState1;
    rng::non_rng::state mNonRngState2;

    size_t mNumXOidx = 0;
    size_t mNumYOidx = 0;
    double dx;
    double dy;
};
} // namespace ray
} // namespace rti
