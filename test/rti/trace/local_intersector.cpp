#include <gtest/gtest.h>

#include <embree3/rtcore.h>
#include <immintrin.h> // AVX
#include <pmmintrin.h>
#include <xmmintrin.h>

#include "rti/trace/local_intersector.hpp"
#include "rti/rng/mt64_rng.hpp"

using namespace rti;
using numeric_type = float;

TEST(local_intersector_test, simple_hit_front) {
  auto ray = RTCRay {};
  ray.org_x = 0; ray.org_y = 0; ray.org_z = 0;
  ray.dir_x = 0; ray.dir_y = 1; ray.dir_z = 0;
  auto disc = util::quadruple<numeric_type> {1, 10, 0, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {0, -1, 0};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_TRUE(intersects);
}

TEST(local_intersector_test, simple_hit_back) {
  auto ray = RTCRay {};
  ray.org_x = 0; ray.org_y = 0; ray.org_z = 0;
  ray.dir_x = 0; ray.dir_y = 1; ray.dir_z = 0;
  auto disc = util::quadruple<numeric_type> {1, 10, 0, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {0, 1, 0};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, simple_miss) {
  auto ray = RTCRay {};
  ray.org_x = 0; ray.org_y = 0; ray.org_z = 0;
  ray.dir_x = 0; ray.dir_y = 1; ray.dir_z = 1;
  auto disc = util::quadruple<numeric_type> {1, 10, 0, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {0, -1, 0};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, simple_parallel_miss) {
  auto ray = RTCRay {};
  ray.org_x = 0; ray.org_y = 0; ray.org_z = 0;
  ray.dir_x = 0; ray.dir_y = 1; ray.dir_z = 0;
  auto disc = util::quadruple<numeric_type> {1, 10, 0, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {-1, 0, 0};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, almost_parallel_miss) {
  auto ray = RTCRay {};
  ray.org_x = 0; ray.org_y = 0; ray.org_z = 0;
  ray.dir_x = 0; ray.dir_y = 1; ray.dir_z = 0;
  auto disc = util::quadruple<numeric_type> {1, 10, 0, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {-1, -0.01, 0};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, some_hit) {
  auto ray = RTCRay {};
  ray.org_x = 2; ray.org_y = 10; ray.org_z = 1.5;
  ray.dir_x = 0.133; ray.dir_y = -1; ray.dir_z = 0.05;
  auto disc = util::quadruple<numeric_type> {3, 5, 2, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {-1, 1.5, 0.5};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_TRUE(intersects);
}

TEST(local_intersector_test, miss_by_shooting_away) {
  auto ray = RTCRay {};
  ray.org_x = 2; ray.org_y = 10; ray.org_z = 1.5;
  ray.dir_x = -0.2; ray.dir_y = 1; ray.dir_z = -0.21;
  auto disc = util::quadruple<numeric_type> {3, 5, 2, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {-1, 1.5, 0.5};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, tight_miss) {
  auto ray = RTCRay {};
  ray.org_x = 2; ray.org_y = 10; ray.org_z = 1.5;
  ray.dir_x = -0.125; ray.dir_y = -1; ray.dir_z = 0.05;
  auto disc = util::quadruple<numeric_type> {3, 5, 2, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {-1, 1.5, 0.5};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, miss_by_starting_at_the_wrong_side_and_shooting_away) {
  auto ray = RTCRay {};
  ray.org_x = 1; ray.org_y = 1; ray.org_z = 1;
  ray.dir_x = 0; ray.dir_y = -1; ray.dir_z = 0;
  auto disc = util::quadruple<numeric_type> {3, 10, 2, 2}; // x, y, z, radius
  auto dnormal = util::triple<numeric_type> {0, 1, 0};
  static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
  util::normalize(*reinterpret_cast<util::triple<numeric_type>*> (&ray.dir_x));
  util::normalize(dnormal);
  auto intersects = trace::local_intersector::intersect(ray, disc, dnormal);
  ASSERT_FALSE(intersects);
}

TEST(local_intersector_test, compare_with_embree) {
  // ATTENTION: Our implementation treats hits from the back differently than Embree!
  auto discs = std::vector<util::quadruple<float> > {
    {-2, -2, 0, 3},
    { 3,  5, 0, 2.76}
  };
  auto normals = std::vector<util::triple<float> > {
    {0.5, 0.2, 1},
    {0.1, -0.155, 1}
  };
  // auto normals = std::vector<util::triple<float> > {
  //   {0, 0, 1},
  //   {0, 0, 1}
  // };
  for (auto& normal : normals) {
    util::normalize(normal);
  }
  auto const& numpoints = discs.size();
  
  struct point_4f_t {
    float xx, yy, zz, radius;
  };
  struct normal_vec_3f_t {
    float xx, yy, zz;
  };
  {  // if multi-threaded, then this block needs to be parallelized.
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  }
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  auto geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);
  auto eediscs = (point_4f_t*)
    rtcSetNewGeometryBuffer(geometry,
                            RTC_BUFFER_TYPE_VERTEX,
                            0, // slot
                            RTC_FORMAT_FLOAT4,
                            sizeof(point_4f_t),
                            numpoints);
  auto eenormals = (normal_vec_3f_t*)
    rtcSetNewGeometryBuffer(geometry,
                            RTC_BUFFER_TYPE_NORMAL,
                            0, // slot
                            RTC_FORMAT_FLOAT3,
                            sizeof(normal_vec_3f_t),
                            numpoints);
  for (size_t idx = 0; idx < numpoints; ++idx) {
    eediscs[idx].xx = discs[idx][0];
    eediscs[idx].yy = discs[idx][1];
    eediscs[idx].zz = discs[idx][2];
    eediscs[idx].radius = discs[idx][3];
    eenormals[idx].xx = normals[idx][0];
    eenormals[idx].yy = normals[idx][1];
    eenormals[idx].zz = normals[idx][2];
  }
  // rtcSetGeometryOccludedFilterFunction(geometry, &occluded_filter);
  // rtcSetGeometryIntersectFilterFunction(geometry, &intersect_filter);
  rtcCommitGeometry(geometry);
  auto scene = rtcNewScene(device);
  auto bbquality = RTC_BUILD_QUALITY_HIGH;
  rtcSetSceneBuildQuality(scene, bbquality);
  //rtcSetSceneFlags(scene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
  auto geometryID = rtcAttachGeometry(scene, geometry);
  rtcCommitScene(scene);

  auto rng = rng::mt64_rng {};
  auto rngstate = rng::mt64_rng::state {};
  auto rtccontext = RTCIntersectContext {};
  rtcInitIntersectContext(&rtccontext);
  auto numsamples = 20 * 1024 * 1024;
  auto numhits = 0;
  for (size_t idx = 0; idx < numsamples; ++idx) {
    auto p1 = util::triple<float> {};
    p1[0] = ((float) rng.get(rngstate) / rng.max()) * 600 - 10; // x 
    p1[1] = ((float) rng.get(rngstate) / rng.max()) * 20 - 10; // y
    p1[2] = 10; // z
    auto p2 = util::triple<float> {};
    p2[0] = ((float) rng.get(rngstate) / rng.max()) * 20 - 10; // x 
    p2[1] = ((float) rng.get(rngstate) / rng.max()) * 20 - 10; // y
    p2[2] = -1; // z
    
    auto rayhit = RTCRayHit {0.f};
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.ray.tnear = 0;
    rayhit.ray.org_x = p1[0];
    rayhit.ray.org_y = p1[1];
    rayhit.ray.org_z = p1[2];
    auto direction = util::diff(p2, p1);
    util::normalize(direction);
    rayhit.ray.dir_x = direction[0];
    rayhit.ray.dir_y = direction[1];
    rayhit.ray.dir_z = direction[2];

    rtcIntersect1(scene, &rtccontext, &rayhit);
    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
      // an Embree hit
      numhits += 1;
      auto& ray = rayhit.ray;
      auto intersectsD0 = trace::local_intersector::intersect(ray, discs[0], normals[0]);
      auto intersectsD1 = trace::local_intersector::intersect(ray, discs[1], normals[1]);
      if (rayhit.hit.primID == 0) {
        // if ( ! intersectsD0 ) {
        //   std::cout << ray.org_x << " " << ray.org_y << " " << ray.org_z << " "
        //             << ray.dir_x << " " << ray.dir_y << " " << ray.dir_z << std::endl;
        //   std::cout << "hitpoint: "
        //             << ray.org_x + ray.dir_x * ray.tfar << " "
        //             << ray.org_y + ray.dir_y * ray.tfar << " "
        //             << ray.org_z + ray.dir_z * ray.tfar << std::endl;
        //   //auto intersectsD0_2 = trace::local_intersector::intersect(ray, discs[0], normals[0]);
        // }
        // ASSERT_TRUE (intersectsD0);
        // ASSERT_FALSE(intersectsD1);
      } else if (rayhit.hit.primID == 1) {
        // if ( ! intersectsD1 ) {
        //   std::cout << "ray: "
        //             << ray.org_x << " " << ray.org_y << " " << ray.org_z << " "
        //             << ray.dir_x << " " << ray.dir_y << " " << ray.dir_z << std::endl;
        //   std::cout << "tfar: " << ray.tfar << std::endl;
        //   std::cout << "hitpoint: "
        //             << ray.org_x + ray.dir_x * ray.tfar << " "
        //             << ray.org_y + ray.dir_y * ray.tfar << " "
        //             << ray.org_z + ray.dir_z * ray.tfar << std::endl;
        //   auto intersectsD1_2 = trace::local_intersector::intersect_debug(ray, discs[1], normals[1]);
        // }
        ASSERT_FALSE(intersectsD0);
        ASSERT_TRUE (intersectsD1);
      } else {
        ASSERT_TRUE(false);
      }
    }
  }
  rtcReleaseGeometry(geometry);
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);
  //std::cout << numhits << " hits out of " << numsamples << " sampled particles" << std::endl;
}
