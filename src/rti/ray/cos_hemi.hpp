#pragma once

#include <boost/math/constants/constants.hpp>

#include "rti/rng/i_rng.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace ray {
  class cos_hemi {
  public:
    template<typename numeric_type, typename random_num_type>
    static rti::util::triple<numeric_type> get(const rti::util::triple<rti::util::triple<numeric_type> >& pBasis,
                                     random_num_type& r1,
                                     random_num_type& r2) {
      assert (0 <= r1 && r1 <= 1 && "Precondition: random number must be in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Precondition: random number must be in the interval [0, 1]");
      // Precondition: pBasis is normalized
      numeric_type epsilon = 1e-6;
      for (auto vec : pBasis) {
        numeric_type length = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
        if ( ! (1 - epsilon <= length && length <= 1 + epsilon) ) {
          // std::cerr << "Assertion error is about to happen with length == " << length << std::endl;
        }
        assert(1 - epsilon <= length && length <= 1 + epsilon && "Precondition: pBasis is normal");
      }

      const numeric_type two_pi = boost::math::constants::two_pi<numeric_type>();
      numeric_type cc1 = sqrt(r2);
      numeric_type cc2 = cos(two_pi * r1) * sqrt(1 - r2);
      numeric_type cc3 = sin(two_pi * r1) * sqrt(1 - r2);

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
