#pragma once

#include "rti/i_rng.hpp"

namespace rti {
  class linux_rng : public rti::i_rng {
  public:

    // Define the state for this RNG
    struct state : public rti::i_rng::i_state {
      unsigned int mSeed;
    };

    // This function simply maps to  rand_r() in stdlib.h. For information on
    // the behaviour see, e.g., the man page.
    uint64_t get(rti::i_rng::i_state* pState) const override final {
      // Precondition:
      assert(typeid(pState) == typeid(state {}) && "Error: precondition violated");

      return (uint64_t) rand_r(&(reinterpret_cast<state*>(pState)->mSeed));
    }
  };
} // namespace rti
