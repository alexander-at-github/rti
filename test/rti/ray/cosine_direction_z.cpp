#include <array>
#include <sstream>

#include <gtest/gtest.h>
#include <gnuplot-iostream.h>

#include "rti/ray/cosine_direction_z.hpp"
#include "rti/rng/cstdlib_rng.hpp"
#include "rti/util/utils.hpp"

TEST(cosine_direction_z_float, plot_w_cstdlib_rng) {
  using numeric_type = float;
  using namespace rti;

  auto cdir = ray::cosine_direction_z<numeric_type> {};
  auto rng = rng::cstdlib_rng {};
  auto rngstate1 = rti::rng::cstdlib_rng::state {1234567890};
  auto rngstate2 = rti::rng::cstdlib_rng::state { 987654321};

  constexpr auto numiter = (size_t) 1e5;

  auto samples = std::array<util::triple<numeric_type>, numiter> {};

  for (size_t idx = 0 ; idx < numiter; ++idx) {
    auto sample = cdir.get(rng, rngstate1, rngstate2);
    samples[idx] = sample;
  }

  auto gpstr = std::stringstream {};
  gpstr
    << "set title 'Samples on a hemisphere weighted by the cosine"
    << " of the zenith angle' \n"
    //<< "set dgrid3d 128, 128, 1 \n"
    //<< "set pm3d \n"
    << "unset hidden3d \n"
    //<< set ticslevel 0.5
    //<< "set view 0,0 \n"
    //<< "set view 90,0 \n"
    << "set view 60,30 \n"
    << "set autoscale \n"
    << "set parametric \n"
    << "set palette color \n"
    //<< "set style data lines \n"
    << "set style data points \n"
    << "set contour \n"
    << "set xlabel 'x-axis' \n"
    << "set ylabel 'y-axis' \n"
    << "set zlabel 'z-axis' \n"
    << "set key off \n"
    << "splot '-' with dots palette title '' \n";
  for (auto const& sample : samples)
    gpstr << sample[0] << " " << sample[1] << " " << sample[2] << "\n";
  Gnuplot gnup {};
  gnup << gpstr.str();
}
