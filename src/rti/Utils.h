#pragma once
#include <embree3/rtcore.h>

#include <iostream>

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
    // "RTC_GEOMETRY_TYPE_TRIANGLE: The index buffer contains an array of three
    // 32-bit indices per triangle (RTC_FORMAT_UINT format)"
    // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
  };
  static void invertnormals(Triangle* ptriangles, size_t numtriangles) {
    for (size_t idx = 0; idx < numtriangles; ++idx) {
      std::swap<uint32_t>(ptriangles[idx].v1, ptriangles[idx].v2);
    }
  }
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
