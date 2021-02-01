#include <iostream>

#include "../device.hpp"
#include "../particle/i_particle.hpp"

class particle : public rti::particle::i_particle<numeric_type> {
public:
  numeric_type
  get_sticking_probability
  (RTCRay& rayin,
   RTCHit& hitin,
   rti::geo::meta_geometry<numeric_type>& geometry,
   rti::rng::i_rng& rng,
   rti::rng::i_rng::i_state& rngstate) override final
  {
    return 0.8;
  }

  void init_new() override final {}
};

int main(int argc, char** argv)
{
  auto points = std::vector<std::array<float, 3> > {{0,0,0}, {1,1,0}};
  auto normals = std::vector<std::array<float, 3> > {{0,0,1}, {0,0,1}};
  auto spacing = std::vector<float> {std::sqrt(4), std::sqrt(4)};

  auto rtidevice = rti::device<float, particle> {};
  rtidevice.set_points(points);
  rtidevice.set_normals(normals);
  rtidevice.set_grid_spacing(spacing);
  rtidevice.set_number_of_rays(1 * 1024 * 1024);
  rtidevice.run();
  auto mcestimates = rtidevice.get_mc_estimates();

  auto seperator = ",";
  auto sep = "";
  for(size_t idx = 0; idx < mcestimates.size(); ++idx) {
    std::cout << sep << mcestimates[idx];
    sep = seperator;
  }
  std::cout << std::endl;

  exit(EXIT_SUCCESS);
}
