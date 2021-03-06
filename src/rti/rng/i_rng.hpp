#pragma once

#include <memory>

namespace rti { namespace rng {
  // Realization of the i_rng interface is intended to be a random number
  // generator together with the definition of a struct which holds all the
  // state (e.g. seeds) which the random number generator needs. That is,
  // this interface defines how a stateless random number generator relates
  // to a state which is held by the user of this interface.
  class i_rng {
    // A pure virtual class / interface
  public:
    virtual ~i_rng() {}

    // A definition of the interface of a state
    struct i_state {
      virtual ~i_state() {}
      virtual std::unique_ptr<i_state> clone() const = 0;
    };

    // A definition of this function will most likely alter the content of its
    // argument.
    virtual uint64_t get(i_state& pState) const = 0;

    //constexpr
    virtual uint64_t min() const = 0;
    //constexpr
    virtual uint64_t max() const = 0;
  };
}} // namespace
