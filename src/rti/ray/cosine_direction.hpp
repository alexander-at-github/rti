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

    using random_num_type = double;
  private:
    void draw_r1_r2(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState)
    {
      r1 = ((random_num_type) pRng.get(pRngState)) / pRng.max();
      r2 = ((random_num_type) pRng.get(pRngState)) / pRng.max();
      assert (0 <= r1 && r1 <= 1 && "Error in computing random number in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Error in computing random number in the interval [0, 1]");
    }

  public:
    rti::util::triple<numeric_type>
    get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final
    {
      draw_r1_r2(pRng, pRngState);
      return rti::ray::cos_hemi::get(mBasis, r1, r2);
    }

    rti::util::pair<random_num_type>
    get_last_random_pair()
    {
      return {r1, r2};
    }

    rti::util::triple<rti::util::triple<numeric_type> >
    get_basis()
    {
      return mBasis;
    }

  private:
    rti::util::triple<rti::util::triple<numeric_type> > mBasis;
    random_num_type r1 = 0;
    random_num_type r2 = 0;
  };
}}
