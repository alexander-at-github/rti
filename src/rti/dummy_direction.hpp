#pragma once

#include "rti/i_direction.hpp"

namespace rti {
  // Since it is a dummy anyway it provides only a float implementation
  class dummy_direction : public i_direction<float> {
    // This class is not thread safe, because of the use of rand_r(),
    // the seed as member variable, and the index variable mXYidx.
  public:

    dummy_direction() :
      mSeed(1) {} // magic number

    dummy_direction(unsigned int pSeed) :
      mSeed(pSeed) {}

    dummy_direction(unsigned int pSeed, size_t pXYidx) :
      mSeed(pSeed),
      mXYidx(pXYidx) {}

    std::unique_ptr<i_direction<float> > clone() const override final {
      return std::make_unique<dummy_direction>(mSeed, mXYidx);
    }

    rti::triple<float> get() override final {
      int rr = rand_r(&mSeed); // stdlib.h
      if (mXYidx >= mXYs.size()) {
        mXYidx = 0;
      }
      float xx = ((float)rr) * 4e-7; // magic number; seems a good value for testing
      float yy = mXYs[mXYidx].frst;
      float zz = mXYs[mXYidx].scnd;
      mXYidx += 1;
      return {xx, yy, zz};
    }

  private:
    unsigned int mSeed;
    size_t mXYidx = 0;
    // 16 coordinates evenly distributed around a circle of radius 10
    std::vector<rti::pair<float> > mXYs
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
} // namespace rti
