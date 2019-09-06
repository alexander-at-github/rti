#pragma once

#include <embree3/rtcore.h>

#include <rti/trace/context.hpp>

namespace rti { namespace surface {
  template<typename Ty>
  class model { // : public rti::surface::i_surface_model
  public:
    model(rti::reflection::i_reflection_model<Ty>& pReflectionModel,
          rti::trace::i_hit_counter& pHitCounter) :
      mReflectionModel(pReflectionModel),
      mHitCounter(pHitCounter) {
      std::cerr << "Warning: This class uses a dummy epsilon value" << std::endl;
    }
    // // A member function returning a pointer to a RTCFilterFunctionN
    // RTCFilterFunctionN get_filter_callback() {
    //   // TODO: Use a lambda!
    //   return &this->filter_fun;
    // }
    //private:
    // Callback filters
    void filter_fun(RTCFilterFunctionNArguments const* args) {
      float epsilon = 0.1; // TODO: FIX
      assert(args->N == 1 && "Assumption: This code can handle only single rays");
      // Question: Can one change the filter function in the RTCObject within the filter?
      // Even if it works, it is probably not a good idea.
      //
      auto valid = args->valid; // a pointer
      valid[0] = 0; // set invalid; causes continuation of ray
      // Cast to our project specific context
      auto context = static_cast<rti::trace::context<Ty>*> (args->context);
      if (context->firstHit) {
        context->firstHitTFar = context.ray[0].tfar;
        context->tfarMax = context.ray[0].tfar + epsilon; // TODO: FIX to proper epsilon value
        return;
      }
      // else
      if (context->ray[0].tfar > context->tfarMax) {
        valid[0] = -1; // do not continue tracing this ray
        return;
      }
      // TODO: decide on reflection
      //context->reflect = true;
    }
  private:
    rti::reflection::i_reflection_model<Ty>& mReflectionModel;
    rti::trace::i_hit_counter& mHitCounter;
  };
}} // namespace
