#include <gtest/gtest.h>

#include "rti/trace/local_intersector.hpp"

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
