#pragma once

#include "rti/i_rng.hpp"

namespace {
  class kiss_rng : public rti::i_rng {
  public:

    // Defines the state of this RNG
    struct state : public rti::i_rng::i_state {

      state(unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4) :
        mS1(p1),
        mS2(p2),
        mS3(p3),
        mS4(p4) {}

      std::unique_ptr<rti::i_rng::i_state> clone() const override final {
        return std::make_unique<state>(mS1, mS2, mS3, mS4);
      }
      unsigned int mS1, mS2, mS3, mS4;
    };

    uint64_t get(rti::i_rng::i_state& pState) const override final {
      assert(false && "Not implemented");
      return 0;
    }
  }
} // namespace rti
