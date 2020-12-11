#pragma once

#include "i_origin.hpp"

namespace rti { namespace ray {
  // Ty is intended to be a numeric type
  template<typename Ty>
  class constant_origin : public rti::ray::i_origin<Ty> {
  public:

    constant_origin(Ty pX, Ty pY, Ty pZ) :
      mX(pX),
      mY(pY),
      mZ(pZ) {}

    rti::util::triple<Ty> get(rti::rng::i_rng& pRng,
                              rti::rng::i_rng::i_state& pRngState1,
                              rti::rng::i_rng::i_state& pRngState2
                              ) const override final {
      return {mX, mY, mZ};
    }

  private:
    Ty mX;
    Ty mY;
    Ty mZ;
  };
}} // namespace rti
