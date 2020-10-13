#pragma once

#include <memory>

#include "i_particle.hpp"

namespace rti { namespace particle {
  template<typename numeric_type>
  class i_particle_factory {
  public:
    virtual ~i_particle_factory() {}
    virtual std::unique_ptr<rti::particle::i_particle<numeric_type> > create() = 0;
  };
}}
