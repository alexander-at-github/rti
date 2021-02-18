#pragma once

#include <x86intrin.h> // vector instruction instrinsics

#include "i_direction.hpp"
#include "i_origin.hpp"
#include "i_source.hpp"
#include "rectangle_origin_z.hpp"

namespace rti::ray {
template <typename Ty>
class patched_source_aligned_dir : public ray::i_source {

public:
  patched_source_aligned_dir(Ty pZval, util::pair<Ty> pC1,
                                  util::pair<Ty> pC2, size_t pNumRays,
                                  size_t pNumXPatches, size_t pNumYPatches,
                                  ray::i_direction<Ty> &pDirection)
      : mZval(pZval), mC1(pC1), mC2(pC2), mNumRays(pNumRays),
        mNumXPatches(pNumXPatches), mNumYPatches(pNumYPatches),
        mPackSize(pNumXPatches * pNumYPatches), mDirection(pDirection),
        mPatchOrigin(0, {0, 0}, {0, 0}) /* some initianlization */
  {
    mpatchdx = ((double)mC2[0] - mC1[0]) / mNumXPatches;
    mpatchdy = ((double)mC2[1] - mC1[1]) / mNumYPatches;
    assert(mpatchdx > 0 && mpatchdy > 0 && "Precondition");
    auto patchC1 = util::pair<Ty>{0, 0};
    auto patchC2 = util::pair<Ty>{(Ty)mpatchdx, (Ty)mpatchdy};
    mPatchOrigin = rectangle_origin_z<float>{mZval, patchC1, patchC2};
    std::cout << "USING patched_source_aligned_dir" << std::endl;
    assert(mNumRays % mPackSize == 0 && "Precondition");
  }

  void fill_ray(RTCRay &pRay, rng::i_rng &pRng, rng::i_rng::i_state &pRngState1,
                rng::i_rng::i_state &pRngState2,
                rng::i_rng::i_state &pRngState3,
                rng::i_rng::i_state &pRngState4) override final {
    if (mPackIdx == 0) {
      mCurrPatchOrgn = mPatchOrigin.get(pRng, pRngState1, pRngState2);
      mCurrDir =  mDirection.get(pRng, pRngState3, pRngState4);
    }
    assert(mPackIdx < mPackSize && "Precondition");
    auto pckix = mPackIdx % mNumXPatches;
    auto pckiy = mPackIdx / mNumXPatches;
    // std::cerr << "packidx == " << packidx << " pckix == " << pckix
    //           << " pckiy == " << pckiy << std::endl;
    // std::cerr << "mpatchdx == " << mpatchdx << " mpatchdy == " << mpatchdy <<
    // std::endl;
    assert(0 <= pckix && pckix < mNumXPatches && "Correctness Assumption");
    assert(0 <= pckiy && pckiy < mNumYPatches && "Correctness Assumption");
    auto mx = mpatchdx * pckix;
    auto my = mpatchdy * pckiy;
    // std::cerr << "mx == " << mx << " my == " << my << std::endl;
    auto orgn = util::triple<Ty>{
        mC1[0] + (Ty)(mpatchdx * pckix) + mCurrPatchOrgn[0],
        mC1[1] + (Ty)(mpatchdy * pckiy) + mCurrPatchOrgn[1], mZval};
    auto constexpr tnear = 1e-4f;
    // std::cerr << "Setting ray origin == " << orgn[0] << ", " << orgn[1] << ",
    // "
    //           << orgn[2] << std::endl;
    reinterpret_cast<__m128 &>(pRay) =
        _mm_set_ps(tnear, (float)orgn[2], (float)orgn[1], (float)orgn[0]);
    auto constexpr time = 0.0f;
    auto const &dir = mCurrDir;
    reinterpret_cast<__m128 &>(pRay.dir_x) =
        _mm_set_ps(time, (float)dir[2], (float)dir[1], (float)dir[0]);
    //
    mPackIdx += 1;
    mPackIdx %= mPackSize;
    // if (mPackIdx >= mPackSize) {
    //   mPackIdx = 0;
    // }
  }

private:
  Ty mZval;
  rti::util::pair<Ty> mC1;
  rti::util::pair<Ty> mC2;
  size_t mNumRays;
  size_t mNumXPatches;
  size_t mNumYPatches;
  size_t mPackSize;
  size_t mPackIdx = 0;
  ray::i_direction<Ty> &mDirection;
  ray::rectangle_origin_z<float> mPatchOrigin;
  util::triple<Ty> mCurrPatchOrgn;
  util::triple<Ty> mCurrDir;
  double mpatchdx;
  double mpatchdy;
};
} // namespace rti::ray
