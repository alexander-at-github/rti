#pragma once

#include <boost/math/constants/constants.hpp>

#include "rti/i_rng.hpp"
#include "rti/utils.hpp"

namespace rti {
  class cos_hemi {
  public:
    template<typename Ty>
    static rti::triple<Ty> get(const rti::triple<rti::triple<Ty> >& pBasis,
                              rti::i_rng& pRng,
                              rti::i_rng::i_state& pRngState) {
      // Precondition: pBasis is normalized
      Ty epsilon = 1e-6;
      for (auto vec : pBasis) {
        Ty length = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
        assert(1 - epsilon <= length && length <= 1 + epsilon && "Precondition: pBasis is normal");
      }

      // TODO: continue here
      Ty r1 = ((Ty) pRng.get(pRngState)) / pRng.max();
      Ty r2 = ((Ty) pRng.get(pRngState)) / pRng.max();
      assert (0 <= r1 && r1 <= 1 && "Error in computing random number in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Error in computing random number in the interval [0, 1]");

      const Ty two_pi = boost::math::constants::two_pi<Ty>();
      Ty cc1 = sqrt(r2);
      Ty cc2 = cos(two_pi * r1) * sqrt(1 - r2);
      Ty cc3 = sin(two_pi * r1) * sqrt(1 - r2);

      auto tt1 = pBasis[0];
      rti::scale(cc1, tt1);
      auto tt2 = pBasis[1];
      rti::scale(cc2, tt2);
      auto tt3 = pBasis[2];
      rti::scale(cc3, tt3);

      auto result = rti::sum(tt1, tt2, tt3);
      assert(rti::is_normalized(result) && "Postcondition");
      return result;
    }
  };
}
