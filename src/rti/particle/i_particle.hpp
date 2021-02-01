#pragma once

#include "../geo/meta_geometry.hpp"
#include "../rng/i_rng.hpp"

namespace rti { namespace particle {
  // An interface for particles. The user of rtilib has to provide an implementation.
  template<typename numeric_type>
  class i_particle {
  public:
    virtual ~i_particle() {}
    // The code in an implementation of process_hit() will be called by multiple threads in parallel.
    // One has to make sure the implementation is thread-safe.
    virtual numeric_type get_sticking_probability(RTCRay& rayin, RTCHit& hitin, rti::geo::meta_geometry<numeric_type>& geometry, rti::rng::i_rng& rng, rti::rng::i_rng::i_state& rngstate) = 0;
    virtual void init_new() = 0;
  };
}}
