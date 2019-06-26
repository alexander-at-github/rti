#pragma once

#include "rti/i_rng.hpp"

namespace {
  class kiss_rng : public rti::i_rng {
  public:

    // Defines the state of this RNG
    struct state : public rti::i_rng::i_state {
      unsigned int mS1, mS2, mS3, mS4;
    };

    uint64_t get(state* pState) const override final {
      assert(false && "Not implemented");
      return 0;
    }
  }
} // namespace rti
