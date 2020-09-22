#include <gtest/gtest.h>

#include "rti/geo/disc_bounding_box_intersector.hpp"
#include "rti/util/utils.hpp"

class disc_bounding_box_intersector_test_1 : public ::testing::Test {
protected:
  rti::util::pair<rti::util::pair<float> > boundingbox
    {rti::util::pair<float> {-1, -1}, rti::util::pair<float> {1, 1}};
  rti::geo::disc_bounding_box_intersector dbbi {boundingbox};
};

class disc_bounding_box_intersector_test_2 : public ::testing::Test {
protected:
  rti::util::pair<rti::util::pair<float> > boundingbox
    {rti::util::pair<float> {-2, -4}, rti::util::pair<float> {8, 1}};
  rti::geo::disc_bounding_box_intersector dbbi {boundingbox};
};

TEST_F(disc_bounding_box_intersector_test_1, disc_fully_inside_r_1) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_fully_inside_r_0_5) {
  auto radius = 0.5f;
  auto disc = rti::util::quadruple<float> {0.2, 0.2, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_fully_outside_x_direction) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {2, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0);
}

TEST_F(disc_bounding_box_intersector_test_1, disc_fully_outside_y_direction) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0, 2, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0);
}

TEST_F(disc_bounding_box_intersector_test_1, disc_inside_parallel_to_bounding_box_x_direction) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0.5, 0, 0, radius};
  auto normal = rti::util::triple<float> {1, 0, 0};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_outside_parallel_to_bounding_box_x_direction) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {1.5, 0, 0, radius};
  auto normal = rti::util::triple<float> {1, 0, 0};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0);
}

TEST_F(disc_bounding_box_intersector_test_1, disc_fully_inside_tilted_close_to_boundary) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0.2, 0, 0, radius};
  auto normal = rti::util::triple<float> {1, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_half_inside_x_direction_right) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {1, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5f * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_a_bit_outside_x_direction_right) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0.5, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  auto alpha = 2 * std::acos(((double) 1.0 /* bounding box */ - disc[0]) / radius);
  auto correct_value =
    (double) radius * radius * rti::util::pi() -
    radius * radius / 2 * (alpha - std::sin(alpha));
  ASSERT_TRUE(std::abs(area_inside - correct_value) < 1e-7);
}

TEST_F(disc_bounding_box_intersector_test_1, disc_almost_entirely_outside_x_direction_right) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {1.9, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  //std::cerr << "area_inside == " << area_inside << std::endl;
  auto alpha = 2 * std::acos(((double) 1.0 /* bounding box */ - disc[0]) / radius);
  //std::cerr << "alpha == " << alpha << std::endl;
  auto correct_value =
    (double) radius * radius * rti::util::pi() -
    radius * radius / 2 * (alpha - std::sin(alpha));
  // std::cerr << "computed value == " << area_inside << std::endl;
  // std::cerr << "correct value == " << correct_value << std::endl;
  ASSERT_TRUE(std::abs(area_inside - correct_value) < 1e-6);
}

TEST_F(disc_bounding_box_intersector_test_1, disc_half_inside_x_direction_left) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {-1, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5f * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_almost_entirely_outside_x_direction_left) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {-1.9, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  //std::cerr << "area_inside == " << area_inside << std::endl;
  auto alpha = 2 * std::acos(((double) 1.0 /* bounding box */ + disc[0]) / radius);
  //std::cerr << "alpha == " << alpha << std::endl;
  auto correct_value =
    (double) radius * radius * rti::util::pi() -
    radius * radius / 2 * (alpha - std::sin(alpha));
  // std::cerr << "computed value == " << area_inside << std::endl;
  // std::cerr << "correct value == " << correct_value << std::endl;
  ASSERT_TRUE(std::abs(area_inside - correct_value) < 1e-6);
}

TEST_F(disc_bounding_box_intersector_test_1, disc_half_inside_y_direction_top) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0, 1, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5f * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_half_inside_y_direction_bottom) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {1, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5f * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_half_inside_tilted_y_direction_bottom) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {1, 0, 0, radius};
  auto normal = rti::util::triple<float> {0, 1, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5f * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_1, disc_inside_parallel_to_bounding_box_y_direction) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {0, 0.8, 0, radius};
  auto normal = rti::util::triple<float> {0, 1, 0};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, radius * radius * (float) rti::util::pi());
}


TEST_F(disc_bounding_box_intersector_test_2, disc_half_inside_x_direction_right) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {8, -1, 0, radius};
  auto normal = rti::util::triple<float> {0, 1, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5 * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_2, disc_fully_inside_x_direction_right) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {7.18, -1, 0, radius}; // 7.18 should be quite close
  auto normal = rti::util::triple<float> {1, 1, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_2, disc_half_inside_x_direction_left) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {-2, -1, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5 * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_2, disc_half_inside_y_direction_bottom) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {3, -4, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5 * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_2, disc_half_inside_y_direction_top) {
  auto radius = 1.0f;
  auto disc = rti::util::quadruple<float> {3, 1, 0, radius};
  auto normal = rti::util::triple<float> {0, 0, 1};
  auto area_inside = dbbi.area_inside(disc, normal);
  ASSERT_EQ(area_inside, 0.5 * radius * radius * (float) rti::util::pi());
}

TEST_F(disc_bounding_box_intersector_test_2, disc_mostly_inside_tilted_y_direction_bottom) {
  auto radius = 1.0f;
  //
  auto distance = 0.7f;
  auto disc = rti::util::quadruple<float> {3, (-4+distance) , 0, radius};
  // the normal is transformed from {0, -1, 1} to turn around the y-axis.
  auto normal = rti::util::triple<float> {0.6, -1, 0.8};
  auto area_inside = dbbi.area_inside(disc, normal);

  auto alpha = 2 * std::acos(std::sqrt(distance * distance * 2) / radius);
  alpha = 2 * rti::util::pi() - alpha;
  //std::cerr << "alpha == " << alpha << std::endl;
  auto correct_value =
    radius * radius / 2 * (alpha - std::sin(alpha));
  // std::cerr << "computed value == " << area_inside << std::endl;
  // std::cerr << "correct value == " << correct_value << std::endl;
  ASSERT_TRUE(std::abs(area_inside - correct_value) < 1e-6);
}

// TEST(test_suite_name, test_name) {
// }