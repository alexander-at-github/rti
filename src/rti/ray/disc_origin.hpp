#pragma once

#include "i_origin.hpp"

namespace rti { namespace ray {
  template<typename numeric_type>
  class disc_origin : public rti::ray::i_origin<numeric_type> {
  public:

    disc_origin(rti::util::triple<numeric_type> coords,
                rti::util::triple<numeric_type> normal,
                numeric_type radius) :
      disc_origin(coords[0], coords[1], coords[2],
                  normal[0], normal[1], normal[2],
                  radius) {}

    disc_origin(numeric_type coordx, numeric_type coordy, numeric_type coordz,
                numeric_type normx, numeric_type normy, numeric_type normz,
                numeric_type radius) :
      coordx(coordx), coordy(coordy), coordz(coordz),
      normx(normx), normy(normy), normz(normz),
      radius(radius) {
      orthonormalbasis = rti::util::get_orthonormal_basis<numeric_type>({normx, normy, normz});
    }

    rti::util::triple<numeric_type>
    get(rti::rng::i_rng& rng, rti::rng::i_rng::i_state& rngstate) const override final
    {
      using internal_numeric_type = double;
      auto r1 = (internal_numeric_type) 1.0;
      auto r2 = (internal_numeric_type) 1.0;
      do {
        r1 = (internal_numeric_type) rng.get(rngstate);
        r2 = (internal_numeric_type) rng.get(rngstate);
        r1 -= rng.max()/2;
        r2 -= rng.max()/2;
        r1 = r1 / (rng.max()/2) * radius;
        r2 = r2 / (rng.max()/2) * radius;
        assert(-radius <= r1 && r1 <= radius &&
               "Error in computin random number in the interval [-radius, +radius]");
        assert(-radius <= r2 && r2 <= radius &&
               "Error in computin random number in the interval [-radius, +radius]");
      } while ( sqrt(r1*r1 + r2*r2) > radius );
      // Remark: orthonormalbasis[0] is a scaled version of the normal vector
      auto b1 = orthonormalbasis[1];
      auto b2 = orthonormalbasis[2];
      assert(rti::util::is_normalized(b1) && rti::util::is_normalized(b2) &&
             "Correctness assumption");
      auto b1s = rti::util::scale((numeric_type) r1, b1);
      auto b2s = rti::util::scale((numeric_type) r2, b2);
      return rti::util::sum(b1s, b2s, {coordx, coordy, coordz});
    }
  private:
    numeric_type coordx, coordy, coordz;
    numeric_type normx, normy, normz;
    numeric_type radius;
    rti::util::triple<rti::util::triple<numeric_type> > orthonormalbasis;
  };
}}
