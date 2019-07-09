#pragma once

#include "rti/i_origin.hpp"

namespace rti {
  // T is intended to be a numeric type
  template<typename T>
  class constant_origin : public i_origin<T> {
  public:

    constant_origin(T pX, T pY, T pZ) :
      mX(pX),
      mY(pY),
      mZ(pZ) {}

    rti::triple<T> get(rti::i_rng& pRng, rti::i_rng::i_state& pRngState) const override final {
      return {mX, mY, mZ};
    }

  private:
    T mX;
    T mY;
    T mZ;
  };
} // namespace rti
