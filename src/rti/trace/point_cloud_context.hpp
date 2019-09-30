#pragma once

#include <vector>
#include <set>

#include <embree3/rtcore.h>

#include "rti/reflection/i_reflection_model.hpp"
#include "rti/reflection/specular.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/i_hit_accumulator.hpp"

// This class needs to be used according to the following protocol. If one does not
// follow the protocol, then the behaviour is undefined.
//
// (0) Call the static function rti::trace::point_cloud_context<Ty>::register_intersect_filter_funs().
//     It will register the intersection filter functions for you in Embree.
//     This needs to be done only once (in one thread).
//     The call of this function needs to be done before the rtcCommitScene() call.
//     (See https://www.embree.org/api.html#scene-object)
// (1) Construct one rti::trace::point_cloud_context<Ty> instance in each thread.
//     Rationale: The functions in rti::trace::point_cloud_context<Ty> are not thread-safe.
// (2) Call init() on each rti::trace::point_cloud_context<Ty> intstance at least once to initialize
//     the Embree data structures contained in the instance.
// (3) Set point_cloud_context.rayWeight = point_cloud_context.INITIAL_RAY_WEIGHT whenever you need it.
// (4) For each tracing operation proceed as follows
//   (4.1) Call point_cloud_context.intersect1() to trace a ray.
//   (4.2) point_cloud_context.tfar contains a reliable value in any case.
//   (4.3) point_cloud_context.reflect is a boolean value specifying if a reflection should happen.
//         If a reflection should happen, then point_cloud_context.rayout contains the new ray.

namespace rti { namespace trace {
  template<typename Ty> // intended to be a numeric type
  class point_cloud_context {
  public:
    // a wrapper aroung struct RTCIntersectContext which attaches additional,
    // rti-project specific data to the context.
    // See https://www.embree.org/api.html#rtcinitintersectcontext
    //
    // Data layout:
    // The first member HAS TO BE the RTC context, that is, context from the Embree library.
    RTCIntersectContext mRtcContext; // not a reference or pointer but full data stored here
    //
    static constexpr float INITIAL_RAY_WEIGHT = 1.0f;
    // =================================================================
    // CHOOSING A GOOD VALUE FOR THE WEIGHT LOWER THRESHOLD IS IMPORTANT
    // =================================================================
    static constexpr float RAY_WEIGHT_LOWER_THRESHOLD = 0.1f;
    static constexpr float RAY_RENEW_WEIGHT = 3 * RAY_WEIGHT_LOWER_THRESHOLD;
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
    rti::geo::absc_point_cloud_geometry<Ty>& mGeometry;
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
    point_cloud_context(unsigned int pGeometryID,
            rti::geo::i_geometry<Ty>& pGeometry,
            rti::reflection::i_reflection_model<Ty>& pReflectionModel,
            rti::trace::i_hit_accumulator<Ty>& pHitAccumulator,
            unsigned int pBoundaryID,
            rti::geo::i_boundary<Ty>& pBoundary,
            rti::reflection::i_reflection_model<Ty>& pBoundaryReflectionModel,
            rti::rng::i_rng& pRng,
            rti::rng::i_rng::i_state& pRngState) :
      mGeometryID(pGeometryID),
      // TODO: FIX the cast
      mGeometry(*dynamic_cast<rti::geo::absc_point_cloud_geometry<Ty>*>(&pGeometry)),
      mReflectionModel(pReflectionModel),
      mHitAccumulator(pHitAccumulator),
      mBoundaryID(pBoundaryID),
      mBoundary(pBoundary),
      mBoundaryReflectionModel(pBoundaryReflectionModel),
      mRng(pRng),
      mRngState(pRngState) {
      //
      std::cerr << "rti::trace::point_cloud_context::point_cloud_context()" << std::endl;
      std::cerr << "Warning: This class uses a dummy epsilon value" << std::endl;
      //
      mGeoHitPrimIDs.reserve(32); // magic number // Reserve some reasonable number of hit elements for one ray
      mGeoHitPrimIDs.clear();
    }

  public:
    static
    void register_intersect_filter_funs(rti::geo::i_geometry<Ty>& pGeometry,
                                        rti::geo::i_boundary<Ty>& pBoundary) {
      // The following cast characterizes a precondition to this function
      auto pPCGeoPointer = dynamic_cast<rti::geo::absc_point_cloud_geometry<Ty>*> (&pGeometry);
      auto& pPCGeo = *pPCGeoPointer;
      RLOG_DEBUG << "register_intersect_filter_funs()" << std::endl;
      rtcSetGeometryIntersectFilterFunction(pPCGeo.get_rtc_geometry(), &filter_fun_geometry);
      rtcSetGeometryIntersectFilterFunction(pBoundary.get_rtc_geometry(), &filter_fun_boundary);
      rtcCommitGeometry(pPCGeo.get_rtc_geometry()); // TODO: needed?
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
      auto  rticonstcontextptr = reinterpret_cast<rti::trace::point_cloud_context<Ty> const*> (cc);
      auto  rticontextptr = const_cast<rti::trace::point_cloud_context<Ty>*> (rticonstcontextptr);
      // The rticontextptr now serves an equal function as the this pointer in a conventional
      // (non-static) member function.
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

      // Debug
      if(rticontextptr->geoNotIntersected) {
        RLOG_DEBUG << "filter_fun_geometry(): FIRST-HIT IS SET ";
      } else {
        RLOG_DEBUG << "filter_fun_geometry(): FIRST-HIT IS *NOT* SET ";
      }
      RLOG_DEBUG << "hitptr->geomID == " << hitptr->geomID << " primID == " << hitptr->primID << std::endl;
      // Non-Debug
      if(rticontextptr->geoNotIntersected) {
        rticontextptr->mGeoHitPrimIDs.clear();
      }
      if ( (! rticontextptr->geoNotIntersected) && rayptr->tfar > rticontextptr->geoTFarMax) {
        // hit is outside of the range which we are interested in
        // stop this ray
        //validptr[0] = -1;
        return;
      }

      // Check whether we hit from the front or from the back of the disc.
      // (It seems Embrees backface culling does not work.)
      //
      // If the dot product of the ray direction and the surface normal is greater than zero, then
      // we hit the back face of the disc.
      if (rti::util::dot_product(rti::util::triple<Ty> {rayptr->dir_x, rayptr->dir_y, rayptr->dir_z},
                                 geometry.get_normal(hitptr->primID)) > 0) {
        if (rticontextptr->geoNotIntersected) {
          validptr[0] = 0; // continue this ray
        } else {
          validptr[0] = -1; // do not continue this ray
        }
        return;
      }

      if (rticontextptr->geoNotIntersected) {
        geoRayout = rticontextptr->mReflectionModel.use(*rayptr, *hitptr, geometry, rng, rngState);
        // set tfar
        rticontextptr->geoFirstHitTFar = rayptr->tfar;
        // determine epsilon:
        // Set epsilon equal to three times the radius of the first disc which is hit by the ray.
        // the fourth element of the primitive is the radius
        auto epsilon = 2 * geometry.get_prim(hitptr->primID)[3];
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
      auto  rticonstcontextptr = reinterpret_cast<rti::trace::point_cloud_context<Ty> const*> (cc);
      auto  rticontextptr = const_cast<rti::trace::point_cloud_context<Ty>*> (rticonstcontextptr);
      // The rticontextptr now serves an equal function as the this pointer in a conventional
      // (non-static) member function.
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

      // //
      // // It may happen that a disc protrudes the boundary. If one would stop searching for intersections
      // // at any boundary hit (that is, at this location in the code), then one could miss discs (the ones
      // // which pertrude) and one could produce an incorrect result on these discs.
      // // Additionally, boundary hits can happen before geometry hits eventhough their sparial ordering
      // // is the other way around (that is, eventhough the geometry is closer to the ray source than the
      // // boundary). Hence, we always want to continue tracing after a boundary hit.
      // validptr[0] = 0; // continue this ray

      /* old code */
      if ( (! rticontextptr->geoNotIntersected) &&
          rayptr->tfar < rticontextptr->geoTFarMax) {
        args->valid[0] = 0; // continue this ray
      }
      // Accept boundary hit.
      args->valid[0] = -1; // accept / stop this rayptr
    }

    private:
    void post_process_intersect(RTCRayHit& pRayHit) {
      RLOG_DEBUG
        << "post_process_intersect(); mGeoHitPrimIDs.size() == " << this->mGeoHitPrimIDs.size() << "; ";
      for(auto const& ee : this->mGeoHitPrimIDs)
        RLOG_DEBUG << ee << " ";
      RLOG_DEBUG
        << std::endl;
      RLOG_DEBUG << "this->rayWeight == " << this->rayWeight << std::endl;

      // // // assert((pRayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) || this->boundNotIntersected &&
      // // //        "Violation of Assumption: \
      // // //         If the Embree intersectin test returns RTC_INVALID_GEOMETRY, then no element of the \
      // // //         bounding box (boundary) has been hit and this->boundNotIntersected still equals true.");
      // // // // Note that (pRayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID && ! this->mGeoHitPrimIDs.empty())
      // // // // may hold since we do not always clear the vector mGeoHitPrimIDs. (We clear it only when we
      // // // // have a first new hit.) The right way to check that the ray did not hit anything is to check
      // // // // the geomID and this->geoNotIntersected.

      // The right way to check that the ray did not hit anything is to check this->geoNotIntersected and
      // this->boundNotIntersected.

      if (this->geoNotIntersected && this->boundNotIntersected) {
        // The ray did not hit anything
        this->reflect = false;
        this->tfar = pRayHit.ray.tfar;
        return;
      }
      this->reflect = true;

      // if ( !this->boundNotIntersected ) {
      //   if ( this->geoNotIntersected ||
      //        ( ! this->geoNotIntersected && // Here we consider that a disc which pertrudes
      //                                       // the bounding box might still be relevant.
      //          this->boundFirstHitTFar < (this->geoFirstHitTFar - (this->geoTFarMax - this->geoFirstHitTFar)))) {
      //     // We hit the bounding box
      //     this->rayout = this->boundRayout;
      //     this->tfar = this->boundFirstHitTFar;
      //     return;
      //   }
      // }
      if ( ! this->boundNotIntersected ) {
        if ( this->geoNotIntersected ||
             this->boundFirstHitTFar < this->geoFirstHitTFar) {
          // We hit bounding box
          this->rayout = this->boundRayout;
          this->tfar = this->boundFirstHitTFar;
          return;
        }
      }

      // else
      assert( ! this->geoNotIntersected && ! this->mGeoHitPrimIDs.empty() && "Assertion");
      // Normal hit on the surface
      this->rayout = this->geoRayout;
      this->tfar = this->geoTFarMax; // we show tfarMax
      // deliver energy onto the surface
      // Note if Embree uses RTC_BUILD_QUALITY_HIGH, then there might exist duplicate
      // hits in this->mGeoHitPrimIDs .
      auto set = std::set<unsigned int> {};
      //auto weightfraction = this->rayWeight / this->mGeoHitPrimIDs.size();
      auto sumvaluedroped = 0.0; // double
      for (size_t idx = 0; idx < this->mGeoHitPrimIDs.size(); ++idx) {
        auto hitprimid = this->mGeoHitPrimIDs[idx];
        if (set.count(hitprimid) != 0) {
          // Duplicate
          continue;
        }
        set.insert(hitprimid);
        auto sticking = mGeometry.get_sticking_coefficient();
        assert(0 <= sticking && sticking <= 1 && "Assumption");
        auto valuetodrop = this->rayWeight * sticking;
        this->mHitAccumulator.use(hitprimid, valuetodrop);
        sumvaluedroped += valuetodrop;
      }
      auto hitprimcount = set.size();
      assert(this->rayWeight >= sumvaluedroped / hitprimcount && "Assertion");
      this->rayWeight -= sumvaluedroped / hitprimcount;
      // We do what is sometimes called Roulette in MC literatur.
      // Jun Liu calls it "rejection controll" in his book.
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
