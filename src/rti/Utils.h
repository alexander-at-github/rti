#pragma once
#include <embree3/rtcore.h>

#include <iostream>
#include <limits>

#include <assert.h>

namespace rti {
class Utils {
  public:
  struct Vertex {
    float xx, yy, zz; // No padding here!
    // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
    // single precision x, y, z floating point coordinates
    // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred
    // from the size of that buffer. "
    // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
  };
  // Embree needs aligned memory:
  //struct alignas(4) Triangle {
  struct Triangle {
    uint32_t v0, v1, v2;
    //int v0, v1, v2;
    // "RTC_GEOMETRY_TYPE_TRIANGLE: The index buffer contains an array of three
    // 32-bit indices per triangle (RTC_FORMAT_UINT format)"
    // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
  };
  static void invertnormals(Triangle* ptriangles, size_t numtriangles) {
    for (size_t idx = 0; idx < numtriangles; ++idx) {
      std::swap<uint32_t>(ptriangles[idx].v1, ptriangles[idx].v2);
    }
  }
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
