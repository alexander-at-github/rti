#pragma once

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
// (3) Set absc_context.rayWeight to initial value by calling init_ray_weight() whenever you need it.
// (4) For each tracing operation proceed as follows
//   (4.1) Call absc_context.intersect1() to trace a ray.
//   (4.2) absc_context.tfar contains a reliable value in any case.
//   (4.3) absc_context.reflect is a boolean value specifying if a reflection should happen.
//         If a reflection should happen, then absc_context.rayout contains the new ray.

namespace rti { namespace trace {
  template<typename numeric_type> // intended to be a numeric type
  class absc_context {
  public:

    virtual ~absc_context() {}
    // Constructor to initialize member variables
    absc_context(float pRayWeight, bool pReflect, rti::util::pair<rti::util::triple<numeric_type> >& pRayout, numeric_type pTfar) :
      rayWeight(pRayWeight),
      reflect(pReflect),
      rayout(pRayout),
      tfar(pTfar) {}
  protected:
    struct context_c_wrapper {
      // wraping the Embree rtc context and the rti context into a C-style struct
      //
      // See https://www.embree.org/api.html#rtcinitintersectcontext
      // Data layout:
      // The first member HAS TO BE the RTC context, that is, context from the Embree library.
      RTCIntersectContext mRtcContext; // not a reference or pointer but full data stored here
      rti::trace::absc_context<numeric_type>& mAbscContext;
      // c-tor
      context_c_wrapper(rti::trace::absc_context<numeric_type>& pAC) :
        mAbscContext(pAC) {}
    };
    context_c_wrapper mContextCWrapper {*this};
  public:

    // virtual void register_intersect_filter_funs(rti::geo::i_geometry<numeric_type>& pGeometry,
    //                                             rti::geo::i_boundary<numeric_type>& pBoundary) = 0;
    /* register_intersect_filter_function cannot be a member of a class (it can only be static).
       One needs to call it in the factory function.
    */
    virtual void intersect1(RTCScene&, RTCRayHit&) = 0;
    virtual void init() = 0;
    virtual void init_ray_weight() = 0;

    virtual bool compute_exposed_areas_by_sampling() = 0;

    // Public data members
    double rayWeight; // initialize to some value
    bool reflect; // initialize to some value
    rti::util::pair<rti::util::triple<numeric_type> >& rayout;
    numeric_type tfar; // initialize to some value
  };
}}
