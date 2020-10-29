#pragma once

namespace rti { namespace mc {

  template<typename numeric_type>
  class rejection_control {
    
  private:

    static
    double get_ray_weight_lower_threshold(numeric_type const& initweight)
    {
      // =================================================================
      // choosing a good value for the weight lower threshold is important
      // =================================================================
      return 0.1 * initweight;
    }

    static
    double get_ray_renew_weight(numeric_type const& initweight)
    {
      return 0.3 * initweight;
    }
    
  public:

    static
    bool check_weight_reweight_or_kill
    (numeric_type& rayweight,
     numeric_type const& initweight,
     rng::i_rng& rng,
     rng::i_rng::i_state& rngstate)
    {
      auto lowerthshld = get_ray_weight_lower_threshold(initweight);
      auto renewweight = get_ray_renew_weight(initweight);
      assert(lowerthshld <= renewweight && "Precondition");
      // We do what is sometimes called Roulette in MC literatur.
      // Jun Liu calls it "rejection controll" in his book.
      // If the weight of the ray is above a certain threshold, we always reflect.
      // If the weight of the ray is below the threshold, we randomly decide to either kill the
      // ray or increase its weight (in an unbiased way).
      if (rayweight >= lowerthshld) {
        // continue the ray without any modification
        return true;
      }
      assert (rayweight < lowerthshld && "Correcntess Assertion");
      // We want to set the weight of (the reflection of) the ray to the value of get_ray_renew_weight().
      // In order to stay  unbiased we kill the reflection with a probability of
      // (1 - current.rayWeight / get_ray_renew_weight()).
      auto rndm = rng.get(rngstate);
      assert (rayweight <= renewweight && "Assumption");
      auto killProbability = 1.0f - rayweight / renewweight;
      if (rndm < (killProbability * rng.max())) {
        // kill the ray
        return false;
      }
      // set_ray_weight_to_new_weight
      rayweight = renewweight;
      // continue ray
      return true;
    }
  };
}}
