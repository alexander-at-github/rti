#pragma once

#include "rti/ray/i_direction.hpp"

namespace rti { namespace ray {
  // Since it is a dummy anyway it provides only a float implementation
  class dummy_direction : public rti::ray::i_direction<float> {
  public:

    dummy_direction() {}

    rti::util::triple<float> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {
      uint64_t rr = pRng.get(pRngState);
      uint64_t mXYidx = pRng.get(pRngState) % mXYs.size();
      float xx = ((float)rr) * 4e-7; // magic number; seems a good value for testing
      float yy = mXYs[mXYidx][0];
      float zz = mXYs[mXYidx][1];
      return {xx, yy, zz};
    }

  private:
    // 16 coordinates evenly distributed around a circle of radius 10
    std::vector<rti::util::pair<float> > mXYs
    {{10.f, 0.f},
       {9.238795325112868f, 3.826834323650898f},
       {7.0710678118654755f, 7.0710678118654755f},
       {3.826834323650898f, 9.238795325112868f},
     {0.f, 10.f},
       {-3.826834323650898f, 9.238795325112868f},
       {-7.0710678118654755f, 7.0710678118654755f},
       {-9.238795325112868f, 3.826834323650898f},
     {-10.f, 0.f},
       {-9.238795325112868f, -3.826834323650898f},
       {-7.0710678118654755f, -7.0710678118654755f},
       {-3.826834323650898f, -9.238795325112868f},
     {0.f, -10.f},
       {3.826834323650898f, -9.238795325112868f},
       {7.0710678118654755f, -7.0710678118654755f},
       {9.238795325112868f, -3.826834323650898f}};
  };
}} // namespace
