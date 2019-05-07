#pragma once

#include <iostream>
#include <limits>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>

namespace rti {
class utils {
  public:
  static RTCRay constructRay(
      float org_x, float org_y, float org_z,
      float dir_x, float dir_y, float dir_z) {
    RTCRay ray;
    ray.org_x = org_x;
    ray.org_y = org_y;
    ray.org_z = org_z;
    ray.dir_x = dir_x;
    ray.dir_y = dir_y;
    ray.dir_z = dir_z;
    // start of ray
    ray.tnear = 0;
    // Maximum length of ray
    ray.tfar = std::numeric_limits<float>::max();
    // TODO: set other necessary fields
    return ray;
    // "When no intersection is found, the ray/hit data is not updated. When an intersection is
    // found, the hit distance is written into the tfar member of the ray and all hit data is set,
    // such as unnormalized geometry normal in object space (Ng hit member), local hit coordinates
    // (u, v hit member), instance ID (instID hit member), geometry ID (geomID hit member), and
    // primitive ID (primID hit member). See Section RTCHit for the hit layout description."
    // Source: https://www.embree.org/api.html#rtcintersect1
  }
  static RTCHit constructHit() {
    return RTCHit{};
  }
  static RTCRayHit constructRayHit(RTCRay ray, RTCHit hit) {
    return RTCRayHit{ray, hit};
  }
};
} // namespace rti
