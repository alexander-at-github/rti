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

    // With unique_ptr one cannot return covariant types
    std::unique_ptr<i_origin<T> > clone() const override final {
      return std::make_unique<constant_origin<T> >(mX, mY, mZ);
    }

    rti::triple<T> get() override final {
      return {mX, mY, mZ};
    }
  private:
    T mX;
    T mY;
    T mZ;
  };
} // namespace rti
