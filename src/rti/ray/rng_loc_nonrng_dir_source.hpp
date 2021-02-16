#pragma once

#include "../rng/i_rng.hpp"
#include "../rng/non_rng.hpp"
#include "i_source.hpp"

namespace rti {
namespace ray {
template <typename numeric_type>
class rng_loc_nonrng_dir_source : public i_source {

public:
  // interesting settings:
  // mNonRng1Max == 99u and a total of 20000000 rays
  // flux values: ~0.98791 (0.01209)
  // Chi-Squared: Chi-Squared == Chi-Squared == 241.196 degrees of freedom == 303
  // mNonRng1Max == 99u and a total of 46080000 rays
  // flux values: ~0.99002 (0.00998)
  // Chi-Squared: Chi-Squared == 280.523 degrees of freedom == 303




  using pair = util::pair<numeric_type>;

  rng_loc_nonrng_dir_source(size_t pNumRays,
                            ray::i_origin<numeric_type> &pOrigin,
                            ray::i_direction<numeric_type> &pDirection)
      : mNumRays(pNumRays), mOrigin(pOrigin), mDirection(pDirection),
        mNonRng1Max(99u), mNonRng1(mNonRng1Max) {
    std::cout << "USING rng_loc_nonrng_dir_source" << std::endl;
    std::cout << "mNumRays == " << mNumRays
              << " mNumRays / (mNonRng1Max + 1) / (mNonRng1Max + 1) == "
              << (double)mNumRays / (mNonRng1Max + 1) / (mNonRng1Max + 1) << std::endl;
    assert(mNumRays % ((mNonRng1Max + 1) * (mNonRng1Max + 1)) == 0 &&
           "Correctness Assumption");
    mNumRaysPerOrg = mNumRays / (mNonRng1Max + 1) / (mNonRng1Max + 1);
  }

  ~rng_loc_nonrng_dir_source() {}

  void fill_ray(RTCRay &pRay, rng::i_rng &pRng, rng::i_rng::i_state &pRngState1,
                rng::i_rng::i_state &pRngState2,
                rng::i_rng::i_state &pRngState3,
                rng::i_rng::i_state &pRngState4) override final {
    // std::cout << "nonrngstates == " << mNonRngState1.mCurr << " "
    //           << mNonRngState2.mCurr << std::endl;
    //if (mNonRngState1.mCurr == 0 && mNonRngState2.mCurr == 0) {
      mCurrOrgn = mOrigin.get(pRng, pRngState1, pRngState2);
      //}
    auto dir = mDirection.get(mNonRng1, mNonRngState1, mNonRngState2);

    pRay.org_x = mCurrOrgn[0];
    pRay.org_y = mCurrOrgn[1];
    pRay.org_z = mCurrOrgn[2];
    pRay.dir_x = dir[0];
    pRay.dir_y = dir[1];
    pRay.dir_z = dir[2];
    if (pRay.dir_z == 0) {
      // pRay.dir_z can lead to infinitely bouncing rays.
      // We do a quick fix.
      pRay.dir_z = -1e-3;
    }
    // advance indices
    mNonRngState1.mCurr += 1;
    if (mNonRngState1.mCurr > mNonRng1.mMax) {
      mNonRngState1.mCurr = 0;
      mNonRngState2.mCurr += 1;
      if (mNonRngState2.mCurr > mNonRng1.mMax) {
        mNonRngState2.mCurr = 0;
        // assert(mNumRaysPerOrgIdx == 0 && "Correctness Assumption");
        // mNumRaysPerOrgIdx = 0;
      }
    }
    // std::cerr << "HERE END" << std::endl;
  }

private:
  // numeric_type mZval;
  // pair mC1;
  // pair mC2;
  size_t mNumRays;
  ray::i_origin<numeric_type> &mOrigin;
  ray::i_direction<numeric_type> &mDirection;

  util::triple<numeric_type> mCurrOrgn;

  size_t mNonRng1Max;
  // size_t mNonRng2Max;
  rng::non_rng mNonRng1;
  // rng::non_rng mNonRng2;
  rng::non_rng::state mNonRngState1;
  rng::non_rng::state mNonRngState2;
  size_t mNumRaysPerOrg;
  // size_t mNumRaysPerOrgIdx = 0;
};
} // namespace ray
} // namespace rti
