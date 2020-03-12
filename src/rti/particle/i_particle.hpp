#pragma once

namespace rti { namespace particle {
  // An interface for particles. The user of rtilib has to provide an implementation.
  template<typename numeric_type>
  class i_particle {
  public:
    virtual ~i_particle() {}
    // The code in an implementation of process_hit() will be called by multiple threads in parallel.
    // One has to make sure the implementation is thread-safe.
    virtual numeric_type process_hit(size_t primID) = 0;
    virtual void init_new() = 0;
  };
}}
