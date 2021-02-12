#pragma once

#include "i_rng.hpp"

namespace rti { namespace rng {
class non_rng : public i_rng {
public:
  class state : public i_state {
  public:
    std::unique_ptr<i_rng::i_state> clone() const override final {
      assert(false);
      return std::make_unique<state> ();
    }
    size_t mCurr = 0;
  };
  non_rng(size_t pMax) : mMax(pMax) {}
  uint64_t get(rti::rng::i_rng::i_state& pState) const override final {
    auto &state = *reinterpret_cast<non_rng::state*> (&pState);
    auto &curr = state.mCurr;
    if (curr > mMax) {
      std::cerr << "ERROR non_rng" << std::endl;
      assert(false);
    }
    return curr;
  }
  uint64_t min() const override final {
    return 0;
  }
  
  uint64_t max() const override final {
    return mMax;
  }
  size_t mMax;
};
  }}
