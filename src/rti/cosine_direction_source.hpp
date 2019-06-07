#pragma once

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <thread> // FIXME: Should be removed in the future. Need to find another solution anyway.

#include <embree3/rtcore.h>

#include "rti/i_ray_source.hpp"

namespace rti {
  class cosine_direction_source : public i_ray_source {
  public:
    // takes a function, which takes no arguments and returns a pair, as parameter and
    // sets seeds mS1 and mS2 to return values of this function
    cosine_direction_source(std::function<rti::pair<unsigned int> ()> fun) :
      mS1(fun().frst),
      mS2(fun().scnd) {}

    cosine_direction_source(unsigned int pSeed1, unsigned int pSeed2) :
      mS1(pSeed1),
      mS2(pSeed2) {}

    cosine_direction_source() :
      mS1((unsigned int) std::hash<std::thread::id>{}(std::this_thread::get_id())),
      mS2((unsigned int) (std::hash<std::thread::id>{}(std::this_thread::get_id())+1)) {
      // std::cerr
      //   << "[Constructor consine_direction_source()] "
      //   <<  std::hash<std::thread::id>{}(std::this_thread::get_id())
      //   << " mS1 == " << mS1 << " mS2 == " << mS2 << std::endl;
    }

    // Sets the direction only
    void fill_ray(RTCRay& pRay) const override final {
      //std::cerr << "[consine_direction_source::set_ray()] " << std::hash<std::thread::id>{}(std::this_thread::get_id()) << std::endl;
      //const double pi = boost::math::constants::pi<double>();
      const double two_pi = boost::math::constants::two_pi<double>();

      // rand_r() returns a value in the interval [0, RAND_MAX]
      // A call to rand_r() modifies the values passed by reference.
      double r1 = ((double) rand_r(&mS1)) / RAND_MAX; // stdlib.h
      double r2 = ((double) rand_r(&mS2)) / RAND_MAX; // stdlib.h
      assert (0 <= r1 && r1 <= 1 && "Error in computing random number in the interval [0, 1]");
      assert (0 <= r2 && r2 <= 1 && "Error in computing random number in the interval [0, 1]");

      // pRay.dir_x = cos(two_pi * r1) * sqrt(1 - r2);
      // pRay.dir_y = sin(two_pi * r1) * sqrt(1 - r2);
      // pRay.dir_z = sqrt(r2);
      pRay.dir_x = sqrt(r2);
      pRay.dir_z = cos(two_pi * r1) * sqrt(1 - r2);
      pRay.dir_y = sin(two_pi * r1) * sqrt(1 - r2);
    }

  private:
    unsigned int mS1;
    unsigned int mS2;
  };
} //namespace rti
