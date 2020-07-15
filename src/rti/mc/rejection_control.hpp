#pragma once

namespace rti { namespace mc {

  template<typename numeric_type>
  class rejection_control {

  public:
    rejection_control
    (rti::rng::i_rng& rng, rti::rng::i_rng::i_state& rngstate) :
      rng(rng),
      rngstate(rngstate) {
      assert(false && "deprecated"):
    }

    template<typename numeric_type_o>
    void
    check_weight_reweight_or_kill
    (rti::trace::absc_context<numeric_type>& context, numeric_type_o weightlow, numeric_type_o renewweight)
    {
      assert(weightlow <= renewweight && "Precondition");
      // We do what is sometimes called Roulette in MC literatur.
      // Jun Liu calls it "rejection controll" in his book.
      // If the weight of the ray is above a certain threshold, we always reflect.
      // If the weight of the ray is below the threshold, we randomly decide to either kill the
      // ray or increase its weight (in an unbiased way).
      context.reflect = true;
      if (context.rayWeight < weightlow) {
        RLOG_DEBUG << "in rti::mc::rejection_control (rayWeight < weightlow) holds" << std::endl;
        // We want to set the weight of (the reflection of) the ray to RAY_NEW_WEIGHT.
        // In order to stay  unbiased we kill the reflection with a probability of
        // (1 - current.rayWeight / renewweight).
        auto rndm = this->rng.get(this->rngstate);
        assert(context.rayWeight <= renewweight && "Assumption");
        auto killProbability = 1.0f - context.rayWeight / renewweight;
        if (rndm < (killProbability * this->rng.max())) {
          kill_the_ray(context);
        } else {
          set_ray_weight_to_new_weight(context, renewweight);
        }
      }
      // If the ray has enough weight, then we do not kill in any case.
    }

  private:
    void kill_the_ray(rti::trace::absc_context<numeric_type>& context)
    {
      context.reflect = false;
      RLOG_TRACE << "K";
    }

    template<typename numeric_type_o>
    void set_ray_weight_to_new_weight(rti::trace::absc_context<numeric_type>& context, numeric_type_o renewweight)
    {
      context.rayWeight = renewweight;
      RLOG_TRACE << "S";
    }

    rti::rng::i_rng& rng;
    rti::rng::i_rng::i_state& rngstate;
  };
}}
