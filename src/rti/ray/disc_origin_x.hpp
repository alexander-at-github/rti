#pragma once

#include "i_origin.hpp"

namespace rti { namespace ray {
  // Ty is intended to be a numeric type
  template<typename Ty>
  class disc_origin_x : public rti::ray::i_origin<Ty> {
  public:

    disc_origin_x(Ty pX, Ty pY, Ty pZ, Ty pRadius) :
      mX(pX),
      mY(pY),
      mZ(pZ),
      mRadius(pRadius) {}

    rti::util::triple<Ty>get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final {
      Ty r1 = 1;
      Ty r2 = 1;
      do {
        r1 = (Ty) pRng.get(pRngState);
        r2 = (Ty) pRng.get(pRngState);
        r1 -= pRng.max()/2;
        r2 -= pRng.max()/2;
        r1 = r1 / (pRng.max()/2) * mRadius;
        r2 = r2 / (pRng.max()/2) * mRadius;
        assert(-mRadius <= r1 && r1 <= mRadius &&
               "Error in computin random number in the interval [-mRadius, +mRadius]");
        assert(-mRadius <= r2 && r2 <= mRadius &&
               "Error in computin random number in the interval [-mRadius, +mRadius]");
      } while ( sqrt(r1*r1 + r2*r2) > mRadius );
      return {mX, mY+r1, mZ+r2};
    }

  private:
    Ty mX;
    Ty mY;
    Ty mZ;
    Ty mRadius;
  };
}} // namespace
