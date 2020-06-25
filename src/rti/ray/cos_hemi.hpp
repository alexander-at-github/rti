#pragma once


#include "rti/rng/i_rng.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace ray {
  class cos_hemi {
  public:
    template<typename Ty>
    static rti::util::triple<Ty> get(const rti::util::triple<rti::util::triple<Ty> >& pBasis,
                              rti::rng::i_rng& pRng,
                              rti::rng::i_rng::i_state& pRngState) {
      // Precondition: pBasis is normalized
      assert(rti::util::is_normalized(pBasis[0]) &&
             rti::util::is_normalized(pBasis[1]) &&
             rti::util::is_normalized(pBasis[2]) && "Precondition");
      Ty epsilon = 1e-6;
      for (auto vec : pBasis) {
        Ty length = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
        if ( ! (1 - epsilon <= length && length <= 1 + epsilon) ) {
          std::cerr
            << "Assertion error is about to happen with length == "
            << length << std::endl;
        }
        assert(1 - epsilon <= length && length <= 1 + epsilon &&
               "Precondition: pBasis is normal");
      }

      Ty r1 = ((Ty) pRng.get(pRngState)) / pRng.max();
      Ty r2 = ((Ty) pRng.get(pRngState)) / pRng.max();
      assert (0 <= r1 && r1 <= 1 && "Error in computing random number in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Error in computing random number in the interval [0, 1]");

      constexpr Ty two_pi = 2 * rti::util::pi();
      Ty cc1 = sqrt(r2);
      Ty cc2 = cos(two_pi * r1) * sqrt(1 - r2);
      Ty cc3 = sin(two_pi * r1) * sqrt(1 - r2);

      auto tt1 = pBasis[0];
      rti::util::scale(cc1, tt1);
      auto tt2 = pBasis[1];
      rti::util::scale(cc2, tt2);
      auto tt3 = pBasis[2];
      rti::util::scale(cc3, tt3);

      auto result = rti::util::sum(tt1, tt2, tt3);
      assert(rti::util::is_normalized(result) && "Postcondition");
      return result;
    }
  };
}} // namespace
