#pragma once

#include <array>
#include <memory>
#include <vector>

#include "particle/i_particle.hpp"
#include "particle/i_particle_factory.hpp"

namespace rti {
  // forward declaration of implementation specific class
  template<typename numeric_type>
  class deviceImpl;

  // The main interface of rtilib
  template<typename numeric_type>
  class device final {

  public:
    device();
    // Because if the Pimpl idiom we need explicite destructor, move constructor, and move assignment.
    // For details see book Effective Modern C++ section Item 22: When using the Pimpl Idiom, define special
    // member functions in the implementation file.
    ~device();
    device(device<numeric_type>&&);
    device<numeric_type>& operator=(device<numeric_type>&&);
    void set_points(std::vector<std::array<numeric_type, 3> >);
    void set_normals(std::vector<std::array<numeric_type, 3> >);
    void set_grid_spacing(std::vector<numeric_type>);
    //void set_sticking_coefficient(numeric_type);
    void set_number_of_rays(size_t);
    // The register_particle_factory() function takes ownership of the particle.
    // That is, you cannot use the unique_ptr after passing it into this function.
    // Call this function in one of two ways:
    //   (a) auto factory = std::make_unique<concrete::factory_class> (\*constructor-arguments*\)
    //       instance.register_particel_factory(std::move(factory))
    //   (b) instance.register_particle_factory(std::unique_ptr<rti::particle::i_particle_factory>
    //         (new concrete::factory_class(\*constructor-arguments*\)))
    void register_particle_factory(std::unique_ptr<rti::particle::i_particle_factory<numeric_type> >);
    void run();
    std::vector<numeric_type> get_mc_estimates();

    std::vector<size_t> get_hit_cnts();

    // non-exposed implementation details
  private:
    std::unique_ptr<rti::deviceImpl<numeric_type> > pimpl;
  };

  // Instantiations of the device class template
  template class device<float>;
  template class device<double>;
}
