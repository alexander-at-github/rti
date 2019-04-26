#pragma once
#include <embree3/rtcore.h>

namespace rti {
class Tools {
  public:
  static RTCRay getRay(
      float org_x, float org_y, float org_z,
      float dir_x, float dir_y, float dir_z) {
    RTCRay ray;
    ray.org_x = org_x;
    ray.org_y = org_y;
    ray.org_z = org_z;
    ray.dir_x = dir_x;
    ray.dir_y = dir_y;
    ray.dir_z = dir_z;
    // TODO: set other necessary fields
    return ray;
  }
  static RTCRayHit getRayHit(RTCRay ray, RTCHit hit) {
    return RTCRayHit{ray, hit};
  }
};
} // namespace rti
