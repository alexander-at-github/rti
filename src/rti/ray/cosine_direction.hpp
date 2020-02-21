#pragma once

#include "rti/ray/cos_hemi.hpp"
#include "rti/ray/i_direction.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace ray {
  template<typename numeric_type>
  class cosine_direction : public rti::ray::i_direction<numeric_type> {

  public:
    cosine_direction(const rti::util::triple<rti::util::triple<numeric_type> > pBasis) :
      mBasis(pBasis) {
      // normalize the basis (in-place)
      rti::util::normalize(mBasis[0]);
      rti::util::normalize(mBasis[1]);
      rti::util::normalize(mBasis[2]);
    }

    static
    cosine_direction construct_in_opposite_direction_of_z_axis()
    {
      return cosine_direction<numeric_type> {
        {rti::util::triple<numeric_type> {0.f, 0.f, -1.f},
         rti::util::triple<numeric_type> {0.f, 1.f,  0.f},
         rti::util::triple<numeric_type> {1.f, 0.f,  0.f}}};
    }

    rti::util::triple<numeric_type>
    get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final
    {
      return rti::ray::cos_hemi::get(mBasis, pRng, pRngState);
    }

  private:
    rti::util::triple<rti::util::triple<numeric_type> > mBasis;
  };
}}
