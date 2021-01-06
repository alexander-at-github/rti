#include <array>
#include <cassert>
#include <sstream>

#include <gtest/gtest.h>
#include <gnuplot-iostream.h>

#include "rti/ray/rectangle_origin_z.hpp"
#include "rti/rng/cstdlib_rng.hpp"
#include "rti/util/utils.hpp"

TEST(rectangle_origin_z_double, histogram_w_cstdlib_rng) {
  using numeric_type = double;
  auto zval = 2.0; // double
  auto c1 = rti::util::pair<numeric_type> {2.0, 2.0};
  auto c2 = rti::util::pair<numeric_type> {12.0, 24.0};
  // Test assumptions
  assert(c1[0] <= c2[0] && "Error in test specification; Test assumption fails");
  assert(c1[1] <= c2[1] && "Error in test specification; Test assumption fails");
  auto org = rti::ray::rectangle_origin_z<numeric_type> {zval, c1, c2};
  auto rng = rti::rng::cstdlib_rng {};
  auto rngstate1 = rti::rng::cstdlib_rng::state {1234567890};
  auto rngstate2 = rti::rng::cstdlib_rng::state { 987654321};

  constexpr auto numxbins = 100u;
  constexpr auto numybins = 100u;
  auto hist = std::array<std::array<size_t, numybins>, numxbins> {};

  auto numiter = (size_t) 1e6;
  for (size_t idx = 0; idx < numiter; ++idx) {
    auto sample = org.get(rng, rngstate1, rngstate2);
    // TODO: write a test that checks that state is changed
    ASSERT_EQ(zval, sample[2]) << "Sampled z-value is incorrect";
    assert(sample[0] - c1[0] >= 0 && "Error in test specification");
    // intentional truncation of number to nearest integer not greate than sample[0]
    auto xbin = size_t ((sample[0] - c1[0]) / (c2[0] - c1[0]) * numxbins);
    if (xbin == numxbins) xbin -= 1;
    // assert(0 <= xbin); // always true
    assert(xbin < numxbins && "Sanity check fails");
    assert(sample[1] - c1[1] >= 0 && "Error in test specification");
    // intentional truncation of number to nearest integer not greate than sample[1]
    auto ybin = size_t ((sample[1] - c1[1]) / (c2[1] - c1[1]) * numybins);
    //assert(0 <= ybin); // always true
    assert(ybin < numybins && "Sanity check fails");
    hist[xbin][ybin]++;
  }

  auto pm3dstr = std::stringstream {};
  auto mapstr = std::stringstream {};

  pm3dstr
    //<< "set term postscript enhanced color \n"
    << "set title 'Histogram of uniform samples on a rectangle' \n"
    << "set pm3d\n"
    << "set key off\n"
    << "splot '-' with pm3d \n";

  mapstr
    << "set title 'Histogram of uniform samples on a rectangle' \n"
    << "set pm3d map\n"
    << "set key off\n"
    << "splot '-' with pm3d \n";

  for (size_t ix1 = 0; ix1 < hist.size(); ++ix1) {
    auto mm = hist[ix1];
    for (size_t ix2 = 0; ix2 < mm.size(); ++ix2) {
      pm3dstr << ix1 << " " << ix2 << " " << mm[ix2] << "\n";
      mapstr << ix1 << " " << ix2 << " " << mm[ix2] << "\n";
    }
    pm3dstr << "\n";
    mapstr << "\n";
  }

  //auto gnuplot = Gnuplot{}; // this syntax does not work with Gnuplot
  Gnuplot gnuplot3d {"gnuplot -persist"};
  Gnuplot gnuplotmap {"gnuplot -persist"};
  gnuplot3d << pm3dstr.str();
  gnuplotmap << mapstr.str();
}
