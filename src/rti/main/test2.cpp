#include <iostream>

#include "../device.hpp"
#include "../particle/i_particle.hpp"
#include "../particle/i_particle_factory.hpp"

template<typename numeric_type>
class particle : public rti::particle::i_particle<numeric_type> {
public:
  numeric_type process_hit(size_t primID, std::array<numeric_type, 3> direction) override final
  {
    return 0.8;
  }

  void init_new() override final
  {
    return;
  }
};

template<typename numeric_type>
class particle_factory : public rti::particle::i_particle_factory<numeric_type> {
public:
  std::unique_ptr<rti::particle::i_particle<numeric_type> > create() override final
  {
    return std::make_unique<particle<numeric_type>>();
  }
};

int main(int argc, char** argv)
{
  auto points = std::vector<std::array<float, 3> > {{0,0,0}, {1,1,0}};
  auto normals = std::vector<std::array<float, 3> > {{0,0,1}, {0,0,1}};
  auto spacing = std::vector<float> {std::sqrt(4), std::sqrt(4)};

  auto rtidevice = rti::device<float> {};
  auto particlefactory = std::make_unique<particle_factory<float> > ();
  rtidevice.set_points(points);
  rtidevice.set_normals(normals);
  rtidevice.set_grid_spacing(spacing);
  rtidevice.set_number_of_rays(1 * 1024 * 1024);
  rtidevice.register_particle_factory(std::move(particlefactory));
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
