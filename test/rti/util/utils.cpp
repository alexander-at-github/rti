#include<gtest/gtest.h>

#include "rti/util/utils.hpp"

TEST(utils, standard_normal_pdf_2d)
{
  auto epsilon = 1e-6;
  auto value_0_0 = rti::util::standard_normal_pdf_2d<double>({0.0, 0.0});
  //std::cerr << "value_0_0 == " << value_0_0 << std::endl;
  ASSERT_TRUE(value_0_0 - 1.0/(2*rti::util::pi) < epsilon) << "value of standard normal density is incorrect";
  auto value_16_m24 = rti::util::standard_normal_pdf_2d<double>({1.6, -2.4});
  //std::cerr << "value_16_m24 == " << value_16_m24 << std::endl;
  ASSERT_TRUE(value_16_m24 - 0.0024840199925583277 < epsilon) << "value of standard normal density is incorrect";
}

TEST(utils, normal_pdf_2d)
{
  auto epsilon = 1e-6;

  auto lastvalue = rti::util::normal_pdf_2d<double>
    ({0.0, 0.0}, // x and y
     {0.0, 0.0}, // mu(x) and mu(y)
     {1.0, 1.0}); // var(x) and var(y)
  ASSERT_TRUE(lastvalue - 1.0/(2*rti::util::pi) < epsilon) << "value of standard normal density is incorrect";

  auto m1 = 5.0;
  auto m2 = 6.0;
  lastvalue = rti::util::normal_pdf_2d<double>
    ({m1, m2},
     {m1, m2},
     {1.0, 1.0});
  ASSERT_TRUE(lastvalue - 1.0/(2*rti::util::pi) < epsilon) << "value of standard normal density is incorrect";

  auto s1 = 2.0;
  auto s2 = 3.0;
  lastvalue = rti::util::normal_pdf_2d<double>
    ({m1, m2},
     {m1, m2},
     {s1, s2});
  // std::cerr << "lastvalue == " << lastvalue << " nvalue == " << 1.0/(2 * rti::util::pi * std::sqrt(s1 * s2)) << std::endl;
  ASSERT_TRUE(lastvalue - 1.0/(2 * rti::util::pi * std::sqrt(s1 * s2)) < epsilon)
    << "value of standard normal density is incorrect";
  // should be around 0.06497473 (using dmvnorm() in R)
  ASSERT_TRUE(lastvalue - 0.06497473 < epsilon)
    << "value of standard normal density is incorrect";

  lastvalue = rti::util::normal_pdf_2d<double>
    ({3.0, 4.2},
     {-0.4, 1.1},
     {3.3, 2.6});
  // reference value calculated with dmvnorm() in R
  ASSERT_TRUE(lastvalue - 0.001485229 < epsilon) << "value of standard normal density is incorrect";
}
