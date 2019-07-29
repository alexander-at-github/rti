#pragma once

#include "rti/ray/i_origin.hpp"

namespace rti { namespace ray {
  // T is intended to be a numeric type
  template<typename T>
  class constant_origin : public rti::ray::i_origin<T> {
  public:

    constant_origin(T pX, T pY, T pZ) :
      mX(pX),
      mY(pY),
      mZ(pZ) {}

    rti::util::triple<T> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final {
      return {mX, mY, mZ};
    }

  private:
    T mX;
    T mY;
    T mZ;
  };
}} // namespace rti
