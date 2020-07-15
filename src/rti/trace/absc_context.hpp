#pragma once

#include <embree3/rtcore.h>

#include "rti/rng/i_rng.hpp"
#include "rti/util/utils.hpp"

// This class needs to be used according to the following protocol. If one does not
// follow the protocol, then the behaviour is undefined.
//
// (0) Call the static function rti::trace::point_cloud_context<numeric_type>::register_intersect_filter_funs().
//     It will register the intersection filter functions for you in Embree.
//     This needs to be done only once (in one thread).
//     The call of this function needs to be done before the rtcCommitScene() call.
//     (See https://www.embree.org/api.html#scene-object)
// (2) Call init() on each rti::trace::absc_context<numeric_type> intstance at least once to initialize
//     the Embree data structures contained in the instance.
// (3) Set the ray weight to a start value by calling init_ray_weight() whenever you need it.
// (4) For each tracing operation proceed as follows
//   (4.1) Call absc_context.intersect1() to trace a ray.
//   (4.2) absc_context.tfar contains a reliable value in any case.
//   (4.3) absc_context.reflect is a boolean value specifying if a reflection should happen.
//         If a reflection should happen, then absc_context.rayout contains the new ray.

namespace rti { namespace trace {
  template<typename numeric_type>
  class absc_context {

  public:

    virtual ~absc_context() {}

    absc_context(bool pReflect, rti::util::pair<rti::util::triple<numeric_type> >& pRayout, numeric_type pTfar, rti::rng::i_rng& rng, rti::rng::i_rng::i_state& rngstate) :
      reflect(pReflect),
      rayout(pRayout),
      tfar(pTfar),
      rng(rng),
      rngstate(rngstate) {
      set_initial_ray_weight(1); // set default weight of a ray
      init_ray_weight();
    }

  protected:

    struct context_c_wrapper {
      // wraping the Embree rtc context and the rti context into a C-style struct
      //
      // See https://www.embree.org/api.html#rtcinitintersectcontext
      // Data layout:
      // The first member HAS TO BE the RTC context, that is, context from the Embree library.
      RTCIntersectContext mRtcContext; // not a reference or pointer but full data stored here
      rti::trace::absc_context<numeric_type>& mAbscContext;

      context_c_wrapper(rti::trace::absc_context<numeric_type>& pAC) :
        mAbscContext(pAC) {}
    };
    context_c_wrapper mContextCWrapper {*this};

  public:

    virtual void intersect1(RTCScene&, RTCRayHit&) = 0;
    virtual void init() = 0;
    virtual bool compute_exposed_areas_by_sampling() = 0;

    void init_ray_weight()
    {
      rayWeight = initialrayweight;
    }

    // Public data members
    double rayWeight; // initialize to some value
    bool reflect; // initialize to some value
    rti::util::pair<rti::util::triple<numeric_type> >& rayout;
    numeric_type tfar; // initialize to some value

  private:

    double initialrayweight;

  public:

    void set_initial_ray_weight(double iw)
    {
      initialrayweight = iw;
    }

  public:

    double get_initial_ray_weight()
    {
      return initialrayweight;
    }

    double get_ray_weight_lower_threshold()
    {
      // =================================================================
      // choosing a good value for the weight lower threshold is important
      // =================================================================
      return 0.1 * initialrayweight;
    }

    double get_ray_renew_weight()
    {
      return 0.3 * initialrayweight;
    }

  protected:

    rti::rng::i_rng& rng;
    rti::rng::i_rng::i_state& rngstate;

  protected:

    void rejection_control_check_weight_reweight_or_kill()
    {
      assert(get_ray_weight_lower_threshold() <= get_ray_renew_weight() && "Precondition");
      // We do what is sometimes called Roulette in MC literatur.
      // Jun Liu calls it "rejection controll" in his book.
      // If the weight of the ray is above a certain threshold, we always reflect.
      // If the weight of the ray is below the threshold, we randomly decide to either kill the
      // ray or increase its weight (in an unbiased way).
      this->reflect = true;
      if (this->rayWeight < get_ray_weight_lower_threshold()) {
        RLOG_DEBUG << "in rejection_control_check_weight_reweight_or_kill() "
                   << "(rayWeight < get_ray_weight_lower_threshold()) holds" << std::endl;
        // We want to set the weight of (the reflection of) the ray to RAY_NEW_WEIGHT.
        // In order to stay  unbiased we kill the reflection with a probability of
        // (1 - current.rayWeight / get_ray_renew_weight()).
        auto rndm = this->rng.get(this->rngstate);
        assert(this->rayWeight <= get_ray_renew_weight() && "Assumption");
        auto killProbability = 1.0f - this->rayWeight / get_ray_renew_weight();
        if (rndm < (killProbability * this->rng.max())) {
          kill_the_ray();
        } else {
          set_ray_weight_to_new_weight(get_ray_renew_weight());
        }
      }
      // If the ray has enough weight, then we do not kill in any case.
    }

  private:

    void kill_the_ray()
    {
      this->reflect = false;
      // RLOG_TRACE << "K";
    }

    void set_ray_weight_to_new_weight(double renewweight)
    {
      this->rayWeight = renewweight;
      // RLOG_TRACE << "S";
    }
  };
}}
