#pragma once

#include "rti/i_origin.hpp"

namespace rti {
  // T is intended to be a numeric type
  template<typename T>
  class disc_origin_x : public i_origin<T> {
  public:

    disc_origin_x(T pX, T pY, T pZ, T pRadius) :
      disc_origin_x(pX, pY, pZ, pRadius, 1) {}

    disc_origin_x(T pX, T pY, T pZ, T pRadius, unsigned int pSeed) :
      mX(pX),
      mY(pY),
      mZ(pZ),
      mRadius(pRadius),
      mSeed(pSeed) {}

    // With unique_ptr one cannot return covariant types
    std::unique_ptr<i_origin<T> > clone() const override final {
      return std::make_unique<disc_origin_x<T> >(mX, mY, mZ, mRadius, mSeed);
    }

    rti::triple<T> get() override final {
      T r1 = 1;
      T r2 = 1;
      do {
        r1 = (T) rand_r(&mSeed); // stdlib.h
        r2 = (T) rand_r(&mSeed); // stdlib.h
        r1 -= RAND_MAX/2;
        r2 -= RAND_MAX/2;
        r1 = r1 / (RAND_MAX/2) * mRadius;
        r2 = r2 / (RAND_MAX/2) * mRadius;
        assert(-mRadius <= r1 && r1 <= mRadius && "Error in computin random number in the interval [-mRadius, +mRadius]");
        assert(-mRadius <= r2 && r2 <= mRadius && "Error in computin random number in the interval [-mRadius, +mRadius]");
      } while ( sqrt(r1*r1 + r2*r2) > mRadius );
      return {mX, mY+r1, mZ+r2};
    }

  private:
    T mX;
    T mY;
    T mZ;
    T mRadius;
    unsigned int mSeed;
  };
} // namespace rti
