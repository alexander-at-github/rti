#include <array>
#include <cassert>
#include <sstream>

#include <gtest/gtest.h>
#include <gnuplot-iostream.h>

#include "rti/ray/rectangle_origin_z.hpp"
#include "rti/rng/cstdlib_rng.hpp"
#include "rti/rng/mt64_rng.hpp"
#include "rti/util/utils.hpp"

using namespace rti;

// Test Fixture
class rectangle_origin_z_double : public ::testing::Test {
protected:
  using nt = double;

  nt zval = 2.0;
  // Corners of testing domain
  util::pair<nt> c1 {2.0, 2.0}; // x, y
  util::pair<nt> c2 {12.0, 24.0};
  // util::pair<nt> c1 {0, 0}; // x, y
  // util::pair<nt> c2 {38.0, 8.0};
  // Sytem under test
  ray::rectangle_origin_z<nt> org {zval, c1, c2};
  
  constexpr static size_t numxbins = 100;
  constexpr static size_t numybins = 100;
  // constexpr static size_t numxbins = 38;
  // constexpr static size_t numybins = 8;
  std::array<std::array<size_t, numybins>, numxbins> histogram {{0u}};

  constexpr static size_t numiter = 1e7;
  
  rectangle_origin_z_double()
  {
    // I am not sure if the histogram values are initialized correctly when
    // the histogram variables is declared. (Execution shows that they are not,
    // at least with my compiler.) Hence, we explicitely initialize it
    // here.
    initialize_histogram();
  }
  
  void run_with_i_rng
  (rng::i_rng& rng,
   rng::i_rng::i_state& rngstate1,
   rng::i_rng::i_state& rngstate2)
  {
    for (size_t idx = 0; idx < numiter; ++idx) {
      auto sample = org.get(rng, rngstate1, rngstate2);
      //std::cout << "sample == " << sample[0] << " " << sample[1] << " " << sample[2] << std::endl;
      ASSERT_EQ(zval, sample[2]) << "sampled z-value must be correct";
      assert(sample[0] - c1[0] >= 0 && "test specification must be correct");
      auto xbin = (size_t) ((sample[0] - c1[0]) / (c2[0] - c1[0]) * numxbins);
      if (xbin == numxbins) xbin -= 1;
      assert(xbin < numxbins && "sanity assertion");
      assert(sample[1] - c1[1] >= 0 && "test specification must be correct");
      auto ybin = (size_t) ((sample[1] - c1[1]) / (c2[1] - c1[1]) * numybins);
      assert(ybin < numybins && "sanity assertion");
      histogram[xbin][ybin] += 1;
    }
  }

  template<typename generator>
  void run_with_uniform_real_distribution
  (generator& rng1,
   std::uniform_real_distribution<>& dist1,
   generator& rng2,
   std::uniform_real_distribution<>& dist2)
  {
    for (size_t idx = 0; idx < numiter; ++idx) {
      auto xx = dist1(rng1);
      auto yy = dist2(rng2);
      auto xbin = (size_t) ((xx - c1[0]) / (c2[0] - c1[0]) * numxbins);
      assert(xbin < numxbins && "sanity assertion");
      auto ybin = (size_t) ((yy - c1[1]) / (c2[1] - c1[1]) * numybins);
      assert(ybin < numybins && "sanity assertion");
      histogram[xbin][ybin] += 1;
    }
  }

  static
  std::string with_gnuplot_header_3d(std::string titlemsg)
  {
    return 
      std::string("set title '") + titlemsg + "'\n" +
      "set pm3d\n" +
      "set key off\n" +
      "splot '-' with pm3d \n";
  }

  static
  std::string with_gnuplot_header_map(std::string titlemsg)
  {
    return
      std::string("set title '") + titlemsg + "'\n" +
      "set pm3d map\n" +
      "set key off\n" +
      "splot '-' with pm3d \n";
  }
  
  void plot(std::string gnuplot_header)
  {
    Gnuplot {"gnuplot -persist"} << gnuplot_header << histogram_to_str();
  }

  util::pair<nt> histogram_min_max()
  {
    auto min = std::numeric_limits<nt>::max();
    auto max = std::numeric_limits<nt>::lowest();
    for (auto const& vv : histogram) {
      for (auto const& ee : vv) {
        if (ee > max) max = ee;
        if (ee < min) min = ee;
      }
    }
    return {min, max};
  }

private:

  void initialize_histogram()
  {
    for (auto& vv : histogram) {
      for (auto& ee : vv) {
        ee = 0u;
      }
    }
  }
  
  std::string histogram_to_str()
  {
    auto str = std::stringstream {};
    for (size_t ix1 = 0; ix1 < histogram.size(); ++ix1) {
      auto mm = histogram[ix1];
      for (size_t ix2 = 0; ix2 < mm.size(); ++ix2) {
        str << ix1 << " " << ix2 << " " << mm[ix2] << "\n";
      }
      str << "\n";
    }
    return str.str();
  }
};

TEST_F(rectangle_origin_z_double, histogram_w_cstdlib_rng) {
  auto rng = rng::cstdlib_rng {};
  auto rngstate1 = rng::cstdlib_rng::state {1234567890};
  auto rngstate2 = rng::cstdlib_rng::state { 987654321};
  run_with_i_rng(rng, rngstate1, rngstate2);
  auto minmax = histogram_min_max();
  auto title = std::string("cstdlib\\_rng; min: ") +
    std::to_string(minmax[0]) + " max: " + std::to_string(minmax[1]);
  plot(with_gnuplot_header_3d(title));
  plot(with_gnuplot_header_map(title));
}

TEST_F(rectangle_origin_z_double, histogram_w_mt64_rng) {
  auto rng = rng::mt64_rng {};
  auto rngstate1 = rng::mt64_rng::state {1234567890};
  auto rngstate2 = rng::mt64_rng::state { 987654321};
  run_with_i_rng(rng, rngstate1, rngstate2);
  auto minmax = histogram_min_max();
  auto title = std::string("mt64\\_rng; min: ") +
    std::to_string(minmax[0]) + " max: " + std::to_string(minmax[1]);
  plot(with_gnuplot_header_3d(title));
  plot(with_gnuplot_header_map(title));
}

TEST_F(rectangle_origin_z_double, histogram_w_mt64_stdlib_uniform_distribution) {
  auto rng1 = std::mt19937_64 {std::mt19937_64::default_seed};
  auto dist1 = std::uniform_real_distribution<nt> {c1[0], c2[0]}; // xmin, xmax
  auto rng2 = std::mt19937_64 {std::mt19937_64::default_seed / 2};
  auto dist2 = std::uniform_real_distribution<nt> {c1[1], c2[1]}; // ymin, ymax
  run_with_uniform_real_distribution(rng1, dist1, rng2, dist2);
  auto minmax = histogram_min_max();
  auto title = std::string("mt64 with std::uniform\\_real\\_distribution; min: ") +
    std::to_string(minmax[0]) + " max: " + std::to_string(minmax[1]);
  plot(with_gnuplot_header_3d(title));
  plot(with_gnuplot_header_map(title));
}

TEST_F(rectangle_origin_z_double, histogram_w_ranlux48_stdlib_uniform_distribution) {
  auto rng1 = std::ranlux48 {(int64_t) 1234567890};
  auto dist1 = std::uniform_real_distribution<nt> {c1[0], c2[0]}; // xmin, xmax
  auto rng2 = std::ranlux48 {(int64_t)  987654321};
  auto dist2 = std::uniform_real_distribution<nt> {c1[1], c2[1]}; // ymin, ymax
  run_with_uniform_real_distribution(rng1, dist1, rng2, dist2);
  auto minmax = histogram_min_max();
  auto title = std::string("ranlux48 with std::uniform\\_real\\_distribution; min: ") +
    std::to_string(minmax[0]) + " max: " + std::to_string(minmax[1]);
  plot(with_gnuplot_header_3d(title));
  plot(with_gnuplot_header_map(title));
}
