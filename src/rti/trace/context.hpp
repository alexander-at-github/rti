#pragma once

#include <vector>

#include <embree3/rtcore.h>

#include "rti/reflection/i_reflection_model.hpp"
#include "rti/reflection/specular.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/i_hit_accumulator.hpp"

// This class needs to be used according to the following protocol. If one does not
// follow the protocol, then the behaviour is undefined.
//
// (0) Call the static function rti::trace::context<Ty>::register_intersect_filter_funs().
//     It will register the intersection filter functions for you in Embree.
//     This needs to be done only once (in one thread).
//     The call of this function needs to be done before the geometries are attached to the
//     scene though. (Is that correct?)
// (1) Construct one rti::trace::context<Ty> instance in each thread.
//     Rationale: The functions in rti::trace::context<Ty> are not thread-safe.
// (2) Call init() on each rti::trace::context<Ty> intstance at least once to initialize
//     the Embree data structures contained in the instance.
// (3) Set context.rayWeight = context.INITIAL_RAY_WEIGHT whenever you need it.
// (4) For each tracing operation proceed as follows
//   (4.1) Call context.intersect1() to trace a ray.
//   (4.2) context.tfar contains a reliable value in any case.
//   (4.3) context.reflect is a boolean value specifying if a reflection should happen.
//         If a reflection should happen, then context.rayout contains the new ray.

namespace rti { namespace trace {
  template<typename Ty> // intended to be a numeric type
  class context {
  public:
    // a wrapper aroung struct RTCIntersectContext which attaches additional,
    // rti-project specific data to the context.
    // See https://www.embree.org/api.html#rtcinitintersectcontext
    //
    // Data layout:
    // The first member is the RTC context, that is, context from the Embree library.
    RTCIntersectContext mRtcContext; // not a reference or pointer but full data stored here
    //
    static constexpr float INITIAL_RAY_WEIGHT = 1.0f;
    // =================================================================
    // CHOOSING A GOOD VALUE FOR THE WEIGHT LOWER THRESHOLD IS IMPORTANT
    // =================================================================
    static constexpr float RAY_WEIGHT_LOWER_THRESHOLD = 0.1f;
    static constexpr float RAY_RENEW_WEIGHT = 0.3f;
    //
    // additional data
    // Here we violate the naming convention (to name a member variable with a string
    // starting with the character 'm') on purpose. Rationale: The context provides
    // (contextual) local data to the outside.
    // TODO: change the naming?
    //
    // geometry related data
  private:
    bool geoNotIntersected = true;
    Ty geoFirstHitTFar = 0; // the type will probably be float since Embree uses float for its RTCRay.tfar
    Ty geoTFarMax = 0;
    rti::util::pair<rti::util::triple<Ty> > geoRayout;
    // boundary related data
    bool boundNotIntersected = true;
    Ty boundFirstHitTFar = 0;
    rti::util::pair<rti::util::triple<Ty> > boundRayout;
    // other data
  public:
    float rayWeight = INITIAL_RAY_WEIGHT;
    bool reflect = false;
    rti::util::pair<rti::util::triple<Ty> >& rayout = geoRayout; // initialize to some value
    Ty tfar = 0; // initialize to some value
    //
    // Class members which will be used in a conventional way.
  private:
    unsigned int mGeometryID;
    rti::geo::i_geometry<Ty>& mGeometry;
    rti::reflection::i_reflection_model<Ty>& mReflectionModel;
    rti::trace::i_hit_accumulator<Ty>& mHitAccumulator;
    unsigned int mBoundaryID;
    rti::geo::i_boundary<Ty>& mBoundary;
    rti::reflection::i_reflection_model<Ty>& mBoundaryReflectionModel;

    rti::rng::i_rng& mRng;
    rti::rng::i_rng::i_state& mRngState;

    // struct hit_t {
    //   unsigned int geomID;
    //   unsigned int primID;
    // };
    // A vector of primitive IDs collected through the filter function filter_fun_geometry() which we
    // will then post process in the post_process_intersection() function.
    // Initialize to some reasonable size. The vector may grow, if needed.
    // Weird initialization necessary
    std::vector<unsigned int> mGeoHitPrimIDs;


    // Constructor
  public:
    context(unsigned int pGeometryID,
            rti::geo::i_geometry<Ty>& pGeometry,
            rti::reflection::i_reflection_model<Ty>& pReflectionModel,
            rti::trace::i_hit_accumulator<Ty>& pHitAccumulator,
            unsigned int pBoundaryID,
            rti::geo::i_boundary<Ty>& pBoundary,
            rti::reflection::i_reflection_model<Ty>& pBoundaryReflectionModel,
            rti::rng::i_rng& pRng,
            rti::rng::i_rng::i_state& pRngState) :
      mGeometryID(pGeometryID),
      mGeometry(pGeometry),
      mReflectionModel(pReflectionModel),
      mHitAccumulator(pHitAccumulator),
      mBoundaryID(pBoundaryID),
      mBoundary(pBoundary),
      mBoundaryReflectionModel(pBoundaryReflectionModel),
      mRng(pRng),
      mRngState(pRngState) {
      std::cerr << "rti::trace::context::context()" << std::endl;
      std::cerr << "Warning: This class uses a dummy epsilon value" << std::endl;
      //
      mGeoHitPrimIDs.reserve(16); // Reserve some reasonable number of hit elements for one ray
      mGeoHitPrimIDs.clear();
    }

  public:
    static
    void register_intersect_filter_funs(rti::geo::i_geometry<Ty>& pGeometry,
                                        rti::geo::i_boundary<Ty>& pBoundary) {
      RLOG_DEBUG << "register_intersect_filter_funs()" << std::endl;
      rtcSetGeometryIntersectFilterFunction(pGeometry.get_rtc_geometry(), &filter_fun_geometry);
      rtcSetGeometryIntersectFilterFunction(pBoundary.get_rtc_geometry(), &filter_fun_boundary);
      rtcCommitGeometry(pGeometry.get_rtc_geometry()); // TODO: needed?
      rtcCommitGeometry(pBoundary.get_rtc_geometry());
    }

  private:
    // Callback filters
    static // static cause C++ allows to pass a function (take the adress of a function)
    // only on static functions. The context object will be passed through the
    // argument.
    void filter_fun_geometry(RTCFilterFunctionNArguments const* args) {
      RLOG_DEBUG << "filter_fun_geometry()" << std::endl;
      assert(args->N == 1 && "Precondition");
      // This function gets a pointer to a context object in args.context
      auto cc = args->context;
      // The following cast also characterizes a precondition to this function
      auto  rticonstcontextptr = reinterpret_cast<rti::trace::context<Ty> const*> (cc);
      auto  rticontextptr = const_cast<rti::trace::context<Ty>*> (rticonstcontextptr);
      //auto& rtccontext = rticontextptr->mRtcContext;
      assert(args->N == 1 && "Precondition: for the cast");
      auto  rayptr = reinterpret_cast<RTCRay*> (args->ray); // ATTENTION: this cast works only if N == 1 holds
      auto  hitptr = reinterpret_cast<RTCHit*> (args->hit); // ATTENTION: this cast works only if N == 1 holds
      auto  validptr = args->valid;
      // Reference all the data in the context with local variables
      // auto& reflect =    rticontextptr->reflect;
      auto& geometryID = rticontextptr->mGeometryID;
      auto& rng =        rticontextptr->mRng;
      auto& rngState =   rticontextptr->mRngState;
      auto& geometry =   rticontextptr->mGeometry;
      auto& geoRayout =  rticontextptr->geoRayout;

      assert(hitptr->geomID == geometryID && "Assumption");

      auto epsilon = 0.02f;                                         // TODO: FIX

      // Debug
      if(rticontextptr->geoNotIntersected) {
        RLOG_DEBUG << "filter_fun_geometry(): FIRST-HIT IS SET ";
      } else {
        RLOG_DEBUG << "filter_fun_geometry(): FIRST-HIT IS *NOT* SET ";
      }
      RLOG_DEBUG << "hitptr->geomID == " << hitptr->geomID << " primID == " << hitptr->primID << std::endl;

      if(rticontextptr->geoNotIntersected) {
        rticontextptr->mGeoHitPrimIDs.clear();
      }
      if ( (! rticontextptr->geoNotIntersected) && rayptr->tfar > rticontextptr->geoTFarMax) {
        // hit is outside of the range which we are interested in
        // stop this ray
        //validptr[0] = -1;
        return;
      }
      if (rticontextptr->geoNotIntersected) {
        geoRayout = rticontextptr->mReflectionModel.use(*rayptr, *hitptr, geometry, rng, rngState);
        // set tfar
        rticontextptr->geoFirstHitTFar = rayptr->tfar;
        rticontextptr->geoTFarMax = rayptr->tfar + epsilon;
        RLOG_DEBUG << "filter_fun_geometry(): geoTFarMax set to " << rticontextptr->geoTFarMax << std::endl;
        rticontextptr->geoNotIntersected = false;
      }
      // Add the hit to the hit-vector
      rticontextptr->mGeoHitPrimIDs.push_back(hitptr->primID);
      validptr[0] = 0; // continue this ray
      return;
    }

    static
    void filter_fun_boundary(RTCFilterFunctionNArguments const* args) {
      RLOG_DEBUG << "filter_fun_boundary()" << std::endl;
      assert(args->N == 1 && "Precondition");
      // This function gets a pointer to a context object in args.context
      auto cc = args->context;
      // The following cast also characterizes a precondition to this function
      auto  rticonstcontextptr = reinterpret_cast<rti::trace::context<Ty> const*> (cc);
      auto  rticontextptr = const_cast<rti::trace::context<Ty>*> (rticonstcontextptr);
      // auto& rtccontext = rticontextptr->mRtcContext;
      assert(args->N == 1 && "Precondition: for the cast");
      auto  rayptr = reinterpret_cast<RTCRay*> (args->ray); // ATTENTION: this cast works only if N == 1 holds
      auto  hitptr = reinterpret_cast<RTCHit*> (args->hit); // ATTENTION: this cast works only if N == 1 holds
      auto  validptr = args->valid;
      // Reference all the data in the context with local variables
      auto& boundaryID = rticontextptr->mBoundaryID;
      auto& rng =        rticontextptr->mRng;
      auto& rngState =   rticontextptr->mRngState;
      auto& boundary =   rticontextptr->mBoundary;
      auto& boundRayout = rticontextptr->boundRayout;

      assert(hitptr->geomID == boundaryID && "Assumption");
      // Debug
      if(rticontextptr->boundNotIntersected) {
        RLOG_DEBUG << "filter_fun_boundary(): FIRST-HIT IS SET ";
      } else {
        RLOG_DEBUG << "filter_fun_boundary(): FIRST-HIT IS *NOT* SET ";
      }
      RLOG_DEBUG << "hitptr->geomID == " << hitptr->geomID << " primID == " << hitptr->primID << std::endl;
      // Non-Debug
      if(rticontextptr->boundNotIntersected) {
        boundRayout = rticontextptr->mBoundaryReflectionModel.use(*rayptr, *hitptr, boundary, rng, rngState);
        rticontextptr->boundFirstHitTFar = rayptr->tfar;
        rticontextptr->boundNotIntersected = false;
      }
      // We always accept boundary hits. That is, we do not change args->valid[0]
    }

    // static // static cause C++ allows to pass a function (take the adress of a function)
    //        // only on static functions. The context object will be passed through the
    //        // argument.
    // void filter_fun(RTCFilterFunctionNArguments const* args) {
    //   assert(args->N == 1 && "Precondition");
    //   // This function gets a pointer to a context object in args.context
    //   auto cc = args->context;

    //   // The following cast also characterizes a precondition to this function
    //   auto  rticonstcontextptr = reinterpret_cast<rti::trace::context<Ty> const*> (cc);
    //   auto  rticontextptr = const_cast<rti::trace::context<Ty>*> (rticonstcontextptr);
    //   // auto& rtccontext = rticontextptr->mRtcContext;
    //   assert(args->N == 1 && "Precondition: for the cast");
    //   auto  rayptr = reinterpret_cast<RTCRay*> (args->ray); // ATTENTION: this cast works only if N == 1 holds
    //   auto  hitptr = reinterpret_cast<RTCHit*> (args->hit); // ATTENTION: this cast works only if N == 1 holds
    //   auto  validptr = args->valid;

    //   // Reference all the data in the context with local variables
    //   //auto& reflect =    rticontextptr->reflect;
    //   auto& rayout =     rticontextptr->rayout;
    //   auto& geometryID = rticontextptr->mGeometryID;
    //   auto& boundaryID = rticontextptr->mBoundaryID;
    //   auto& rng =        rticontextptr->mRng;
    //   auto& rngState =   rticontextptr->mRngState;
    //   auto& geometry =   rticontextptr->mGeometry;
    //   auto& boundary =   rticontextptr->mBoundary;

    //   auto epsilon = 0.02f;                                         // TODO: FIX

    //   if(rticontextptr->geoNotIntersected) {
    //     RLOG_DEBUG << "filter_fun(): FIRST-HIT IS SET ";
    //     //reflect = false;
    //     rticontextptr->mGeoHitPrimIDs.clear();
    //   } else {
    //     // Debug
    //     RLOG_DEBUG << "filter_fun(): FIRST-HIT IS *NOT* SET ";
    //   }
    //   RLOG_DEBUG << "geomID == " << hitptr->geomID << " primID == " << hitptr->primID << std::endl;

    //   //auto& geomID = RTCHitN_geomID(hitptr, args->N, 0);
    //   auto& geomID = hitptr->geomID;
    //   assert((geomID == RTC_INVALID_GEOMETRY_ID || geomID == boundaryID || geomID == geometryID) &&
    //          "Precondition");
    //   if (geomID == RTC_INVALID_GEOMETRY_ID) {
    //     // No hit
    //     // the intersection function will return automatically when the member "validptr" is not changed
    //     return;
    //   }
    //   // ray hit

    //   // Four cases: (1.1) first boundary hit, (1.2) non-first boundary hit
    //   //             (2.1) first geometry hit, (2.2) non-first geometry hit
    //   if (geomID == boundaryID) {
    //     RLOG_DEBUG << "filter_fun(): geomID == boundaryID holds" << std::endl;
    //     // stop this ray
    //     //validptr[0] = -1;
    //     if (rticontextptr->geoNotIntersected) {
    //       // Whether or not to reflect is decided on the first hit
    //       rticontextptr->geoNotIntersected = false;
    //       rayout = rticontextptr->mBoundaryReflectionModel.use(*rayptr, *hitptr, boundary, rng, rngState);
    //       //reflect = true;
    //       //
    //       // Also set geoTFarMax for correct ray-logging with rti::util::RAYLOG
    //       rticontextptr->geoTFarMax = rayptr->tfar + epsilon;

    //     }
    //     return;
    //   }
    //   // else
    //   assert(geomID == geometryID && "Assumption");
    //   if ( (! rticontextptr->geoNotIntersected) && rayptr->tfar > rticontextptr->geoTFarMax) {
    //     // hit is outside of the range which we are interested in
    //     // stop this ray
    //     //validptr[0] = -1;
    //     return;
    //   }
    //   if (rticontextptr->geoNotIntersected) {
    //     rayout = rticontextptr->mReflectionModel.use(*rayptr, *hitptr, geometry, rng, rngState);
    //     //reflect = true;
    //     // We will recheck the reflection decission in the post_process_intersect() function

    //     // set tfar
    //     rticontextptr->geoFirstHitTFar = rayptr->tfar;
    //     rticontextptr->geoTFarMax = rayptr->tfar + epsilon;
    //     RLOG_DEBUG << "filter_fun(): geoTFarMax set to " << rticontextptr->geoTFarMax << std::endl;

    //     rticontextptr->geoNotIntersected = false;
    //   }
    //   //mHitcounter.use(pRayhit);
    //   // Add the hit to the hit-vector
    //   rticontextptr->mGeoHitPrimIDs.push_back({hitptr->geomID, hitptr->primID});

    //   validptr[0] = 0; // continue this ray
    //   return;
    // }

    private:
    void post_process_intersect(RTCRayHit& pRayHit) {
      RLOG_DEBUG
        << "post_process_intersect(); mGeoHitPrimIDs.size() == " << this->mGeoHitPrimIDs.size() << "; ";
      for(auto const& ee : this->mGeoHitPrimIDs)
        RLOG_DEBUG << ee << " ";
      RLOG_DEBUG
        << std::endl;
      RLOG_DEBUG << "this->rayWeight == " << this->rayWeight << std::endl;

      assert((pRayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) || this->boundNotIntersected &&
             "Violation of Assumption: \
              If the Embree intersectin test returns RTC_INVALID_GEOMETRY, then no element of the \
              bounding box (boundary) has been hit and this->boundFistHit still equals true.");
      // Note that (pRayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID && ! this->mGeoHitPrimIDs.empty())
      // can hold since we do not always clear the vector mGeoHitPrimIDs. (We clear it only when we
      // have a first new hit.) The right way to check the ray did not hit anything is to check
      // the geomID and this->geoNotIntersected.
      if (pRayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID && this->geoNotIntersected) {
        // The ray did not hit anything
        this->reflect = false;
        this->tfar = pRayHit.ray.tfar;
        return;
      }
      this->reflect = true;
      if ( !this->boundNotIntersected ) {
        if ( this->geoNotIntersected || this->boundFirstHitTFar < this->geoFirstHitTFar) {
          // We hit the bounding box
          this->rayout = this->boundRayout;
          this->tfar = this->boundFirstHitTFar;
          return;
        }
      }
      // else
      assert( ! this->geoNotIntersected && ! this->mGeoHitPrimIDs.empty() && "Assertion");
      // Normal hit on the surface
      this->rayout = this->geoRayout;
      //this->tfar = this->geoFirstHitTFar;
      this->tfar = this->geoTFarMax; // we show tfarMax
      // deliver energy onto the surface
      auto weightfraction = this->rayWeight / this->mGeoHitPrimIDs.size();
      for (size_t idx = 0; idx < this->mGeoHitPrimIDs.size(); ++idx) {
        auto hitprimid = this->mGeoHitPrimIDs[idx];
        auto sticking = mGeometry.get_sticking_coefficient(); // TODO
        assert(0 <= sticking && sticking <= 1 && "Assumption");
        auto valuetodrop = weightfraction * sticking;
        this->mHitAccumulator.use(hitprimid, valuetodrop);
        this->rayWeight -= valuetodrop;
      }
      // We do what is called Roulette in MC literatur.
      // If the weight of the ray is above a certain threshold, we always reflect.
      // If the weight of the ray is below the threshold, we randomly decide to either kill the
      // ray or increase its weight (in an unbiased way).
      if (this->rayWeight < RAY_WEIGHT_LOWER_THRESHOLD) {
        RLOG_DEBUG << "in post_process_intersect() (rayWeight < RAY_WEIGHT_LOWER_THRESHOLD) holds" << std::endl;
        // We want to set the weight of (the reflection of) the ray to RAY_NEW_WEIGHT.
        // In order to stay  unbiased we kill the reflection with a probability of
        // (1 - rticontextptr->rayWeight / RAY_RENEW_WEIGHT).
        auto rndm = this->mRng.get(this->mRngState);
        assert(this->rayWeight < RAY_RENEW_WEIGHT && "Assumption");
        auto killProbability = 1.0f - this->rayWeight / RAY_RENEW_WEIGHT;
        if (rndm < (killProbability * this->mRng.max())) {
          // kill the ray
          this->reflect = false;
        } else {
          // The ray survived the roulette: set the ray weight to the new weight
          //reflect = true;
          this->rayWeight = RAY_RENEW_WEIGHT;
        }
      }
      // If the ray has enough weight, then we reflect it in any case.
    }
  public:
    void intersect1(RTCScene& pScene, RTCRayHit& pRayHit) {
      // prepare
      pRayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
      pRayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
      this->geoNotIntersected = true;
      this->boundNotIntersected = true;
      // Embree intersect
      // This call uses our intersection functions which must have been registered with
      // register_intersect_filter_funs() .
      auto rtcContextPtr = reinterpret_cast<RTCIntersectContext*> (this);
      rtcIntersect1(pScene, rtcContextPtr, &pRayHit);
      // post process
      this->post_process_intersect(pRayHit);
    }

    void init() {
      // // RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT flag uses an optimized traversal
      // // algorithm for incoherent rays (default)
      rtcInitIntersectContext(&this->mRtcContext);
    }
  };
}} // namespace
