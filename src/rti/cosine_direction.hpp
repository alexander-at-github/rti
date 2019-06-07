#pragma once

#include <boost/math/constants/constants.hpp>

#include "rti/i_direction.hpp"

namespace rti {
  template<typename T>
  class cosine_direction : public i_direction<T> {
    // This class is not thread safe.
  public:
    enum class alignment_t {X, Y, Z};

    cosine_direction() :
      cosine_direction(defaultAlignment) {}

    cosine_direction(alignment_t pAlignment) :
      cosine_direction(pAlignment, 1, 2) {} // magic numbers

    cosine_direction(unsigned int pSeed1, unsigned int pSeed2) :
      cosine_direction(defaultAlignment, pSeed1, pSeed2) {}

    cosine_direction(alignment_t pAlignment, unsigned int pSeed1, unsigned int pSeed2) :
      mAlignment(pAlignment),
      mSeed1(pSeed1),
      mSeed2(pSeed2) {}

    // Note: with unique_ptr one cannot return covariant types
    std::unique_ptr<i_direction<T> > clone() const override final {
      return std::make_unique<cosine_direction<T> >(mAlignment, mSeed1, mSeed2);
    }

    // returns a random direction
    rti::triple<T> get() override final {
      const double two_pi = boost::math::constants::two_pi<double>();
      // rand_r() returns a value in the interval [0, RAND_MAX]
      // A call to rand_r() modifies the values passed by reference.
      double r1 = ((double) rand_r(&mSeed1)) / RAND_MAX; // stdlib.h
      double r2 = ((double) rand_r(&mSeed2)) / RAND_MAX; // stdlib.h
      assert (0 <= r1 && r1 <= 1 && "Error in computing random number in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Error in computing random number in the interval [0, 1]");
      T xx = sqrt(r2);
      T yy = cos(two_pi * r1) * sqrt(1 - r2);
      T zz = sin(two_pi * r1) * sqrt(1 - r2);
      if (mAlignment == alignment_t::X) {
        return {xx, yy, zz};
      }
      if (mAlignment == alignment_t::Y) {
        // Swap coordinates
        return {yy, xx, zz};
      }
      assert(mAlignment == alignment_t::Z && "Error, unexpected direction alignment");
      // Swap coordinates
      return {zz, yy, xx};
    }

  private:
    // Axis alignment
    alignment_t mAlignment;
    constexpr static alignment_t defaultAlignment = alignment_t::X;
    // Seeds for rand_r() in stdlib.h
    unsigned int mSeed1;
    unsigned int mSeed2;
  };
} // namespace rti
