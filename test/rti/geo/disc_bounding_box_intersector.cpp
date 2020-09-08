#include <gtest/gtest.h>

#include "rti/geo/disc_bounding_box_intersector.hpp"
#include "rti/util/utils.hpp"

TEST(disc_bounding_box_intersector, disc_fully_inside) {
  auto bounding_box = rti::util::pair<rti::util::pair<float> >
    {rti::util::pair<float> {-1, -1}, rti::util::pair<float> {1, 1}};
  auto dbi = rti::geo::disc_bounding_box_intersector {bounding_box};
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, (float) rti::util::pi()) << "Computed area incorrect; disc completely inside the box";
}


// TEST(test_suite_name, test_name) {
// }
