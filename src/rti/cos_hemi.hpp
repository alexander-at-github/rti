#pragma once

#include <boost/math/constants/constants.hpp>

#include "rti/i_rng.hpp"
#include "rti/utils.hpp"

namespace rti {
  class cos_hemi {
  public:
    template<typename T>
    static rti::triple<T> get(const rti::triple<rti::triple<T> >& pBasis,
                              rti::i_rng& pRng,
                              rti::i_rng::i_state& pRngState) {
      // Precondition: pBasis is normalized
      float epsilon = 1e-6;
      for (auto vec : pBasis.get_iterable()) {
        T length = std::sqrt(vec.frst * vec.frst + vec.scnd * vec.scnd + vec.thrd * vec.thrd);
        assert(1 - epsilon <= length && length <= 1 + epsilon && "Precondition: pBasis is normal");
      }

      double r1 = ((double) pRng.get(pRngState)) / pRng.max();
      double r2 = ((double) pRng.get(pRngState)) / pRng.max();
      assert (0 <= r1 && r1 <= 1 && "Error in computing random number in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Error in computing random number in the interval [0, 1]");

      const double two_pi = boost::math::constants::two_pi<double>();
      T cc1 = sqrt(r2);
      T cc2 = cos(two_pi * r1) * sqrt(1 - r2);
      T cc3 = sin(two_pi * r1) * sqrt(1 - r2);

      auto tt1 = pBasis.frst;
      rti::scale(tt1, cc1);
      auto tt2 = pBasis.scnd;
      rti::scale(tt2, cc2);
      auto tt3 = pBasis.thrd;
      rti::scale(tt3, cc3);

      return rti::sum(tt1, tt2, tt3);
      //assert(false && "Not implemented");
    }
  };
}
