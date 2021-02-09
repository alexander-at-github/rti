#pragma once

#include <random>

#include "i_rng.hpp"

namespace rti { namespace rng {
  class mt64_rng_v2 : public rti::rng::i_rng {
  public:
    // define state of this RNG
    struct state : public rti::rng::i_rng::i_state {
      state() : state(std::mt19937_64::default_seed) {}
      state(unsigned int pSeed) : mMT(pSeed) {}
    public:
      state(std::mt19937_64 pMT) : mMT(pMT) {}
    public:
      std::unique_ptr<rti::rng::i_rng::i_state> clone() const override final {
        // I am not sure what state the new MT-rand-generator will have
        return std::make_unique<state> (mMT); // works without "new" keyword, works only with public constructor
        //return std::unique_ptr<state> (new state (mMT)); // works with private constructor, works only with "new" keyword
      }
      std::mt19937_64 mMT;
    };

    // define the functions of the RNG
    uint64_t get(rti::rng::i_rng::i_state& pState) const override final {
      // Precondition:
      // The parameter pState needs to be of type rti::rng::mt64_rng_v2::state.
      // This sentence is verified in the following assertion.
      auto stateObjectForAssertion = state {};
      //std::cout << "*pState: " << typeid(*pState).name() << std::endl << std::endl;
      //std::cout
      // << "stateObjectForAssertion: "
      // << typeid(stateObjectForAssertion).name()
      // << std::endl << std::endl;
      assert(typeid(pState) == typeid(stateObjectForAssertion) && "Error: precondition violated");
      auto castStatePtr = reinterpret_cast<state*>(&pState);
      return (uint64_t) castStatePtr->mMT(); // call operator()() function on the mersenne twister
    }

    // constexpr
    uint64_t min() const override final {
      return std::mt19937_64::min();
    }

    // constexpr
    uint64_t max() const override final {
      return std::mt19937_64::max();
    }
  };
}} // namespace
