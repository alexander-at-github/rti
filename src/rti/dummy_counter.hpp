#pragma once

#include "rti/i_hit_counter.hpp"

namespace rti {
  class dummy_counter : public i_hit_counter {
  public:

    dummy_counter() {}

    dummy_counter(dummy_counter const& pDc) {}

    // dummy_counter& operator=(dummy_counter const& pOther) {
    //   if (this != &pOther) {
    //     // copy stuff from pOther to this
    //   }
    //   return *this;
    // }

    dummy_counter(dummy_counter const& pDc1, dummy_counter const& pDc2) {}

    void use(const RTCRayHit& pRayhit) override final {
      return;
    }

  };
} // namespace rti
