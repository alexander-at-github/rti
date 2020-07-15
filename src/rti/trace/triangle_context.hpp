#pragma once

#include <vector>
#include <set>

#include <embree3/rtcore.h>

#include "rti/reflection/i_reflection_model.hpp"
#include "rti/reflection/specular.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/absc_context.hpp"
#include "rti/trace/i_hit_accumulator.hpp"

// This class needs to be used according to a protocol! See base class rti::trace::absc_context

namespace rti { namespace trace {
  template<typename numeric_type> // intended to be a numeric type
  class triangle_context : public rti::trace::absc_context<numeric_type> {

  private:

    // geometry related data
    bool geoNotIntersected = true;
    numeric_type geoFirstHitTFar = 0; // the type will probably be float since Embree uses float for its RTCRay.tfar
    numeric_type geoTFarMax = 0;
    rti::util::pair<rti::util::triple<numeric_type> > geoRayout {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // boundary related data
    bool boundNotIntersected = true;
    numeric_type boundFirstHitTFar = 0;
    rti::util::pair<rti::util::triple<numeric_type> > boundRayout {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  private:

    unsigned int mGeometryID = RTC_INVALID_GEOMETRY_ID; // initialize to some useful value
    rti::geo::triangle_geometry<numeric_type>& mGeometry;
    rti::reflection::i_reflection_model<numeric_type>& mReflectionModel;
    rti::trace::i_hit_accumulator<numeric_type>& mHitAccumulator;
    unsigned int mBoundaryID = RTC_INVALID_GEOMETRY_ID; // initialize to some useful value
    rti::geo::i_boundary<numeric_type>& mBoundary;
    rti::reflection::i_reflection_model<numeric_type>& mBoundaryReflectionModel;

    // A vector of primitive IDs collected through the filter function filter_fun_geometry() which we
    // will then post process in the post_process_intersection() function.
    // Initialize to some reasonable size. The vector may grow, if needed.
    // Weird initialization necessary
    std::vector<unsigned int> mGeoHitPrimIDs {};

  public:

    triangle_context(unsigned int pGeometryID,
            rti::geo::triangle_geometry<numeric_type>& pGeometry,
            rti::reflection::i_reflection_model<numeric_type>& pReflectionModel,
            rti::trace::i_hit_accumulator<numeric_type>& pHitAccumulator,
            unsigned int pBoundaryID,
            rti::geo::i_boundary<numeric_type>& pBoundary,
            rti::reflection::i_reflection_model<numeric_type>& pBoundaryReflectionModel,
            rti::rng::i_rng& pRng,
            rti::rng::i_rng::i_state& pRngState) :
      // initialize members of virtual base class
      rti::trace::absc_context<numeric_type>(false, geoRayout, 0, pRng, pRngState), // initialize to some values
      mGeometryID(pGeometryID),
      mGeometry(pGeometry),
      mReflectionModel(pReflectionModel),
      mHitAccumulator(pHitAccumulator),
      mBoundaryID(pBoundaryID),
      mBoundary(pBoundary),
      mBoundaryReflectionModel(pBoundaryReflectionModel) {
      mGeoHitPrimIDs.reserve(32); // magic number // Reserve some reasonable number of hit elements for one ray
      mGeoHitPrimIDs.clear();
    }

  public:

    static
    void register_intersect_filter_funs(rti::geo::i_geometry<numeric_type>& pGeometry,
                                        rti::geo::i_boundary<numeric_type>& pBoundary)
    {
      // The following cast characterizes a precondition to this function
      auto pPCGeoPointer = dynamic_cast<rti::geo::triangle_geometry<numeric_type>*> (&pGeometry);
      auto& pPCGeo = *pPCGeoPointer;
      RLOG_DEBUG << "register_intersect_filter_funs()" << std::endl;
      rtcSetGeometryIntersectFilterFunction(pPCGeo.get_rtc_geometry(), &filter_fun_geometry);
      rtcSetGeometryIntersectFilterFunction(pBoundary.get_rtc_geometry(), &filter_fun_boundary);
      rtcCommitGeometry(pPCGeo.get_rtc_geometry());
      rtcCommitGeometry(pBoundary.get_rtc_geometry());
    }

  private:

    // Callback filters
    static // static cause C++ allows to pass a function (take the adress of a function)
    // only on static functions. The context object will be passed through the
    // argument.
    void filter_fun_geometry(RTCFilterFunctionNArguments const* args)
    {
      RLOG_DEBUG << "filter_fun_geometry()" << std::endl;
      assert(args->N == 1 && "Precondition");
      // This function gets a pointer to a context object in args.context
      auto cc = args->context;

      //std::cerr << "filter_fun_geometry(): the address of the rtc context: " << cc << std::endl;
      auto ccnonconst = const_cast<RTCIntersectContext*>(cc);
      auto rtiabscontextptr = &reinterpret_cast<typename rti::trace::absc_context<numeric_type>::context_c_wrapper*> (ccnonconst)->mAbscContext;
      auto rticontextptr = reinterpret_cast<rti::trace::triangle_context<numeric_type>*> (rtiabscontextptr);
      //std::cerr << "filter_fun_geometry(): address of the rti context: " << rticontextptr << std::endl;

      // The rticontextptr now serves an equal function as the this pointer in a conventional
      // (non-static) member function.
      assert(args->N == 1 && "Precondition: for the cast");
      auto  rayptr = reinterpret_cast<RTCRay*> (args->ray); // ATTENTION: this cast works only if N == 1 holds
      auto  hitptr = reinterpret_cast<RTCHit*> (args->hit); // ATTENTION: this cast works only if N == 1 holds
      auto  validptr = args->valid;
      // Reference all the data in the context with local variables
      // auto& reflect =    rticontextptr->reflect;
      auto& geometryID = rticontextptr->mGeometryID;
      auto& rng =        rticontextptr->rng;
      auto& rngState =   rticontextptr->rngstate;
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
      geoRayout = rticontextptr->mReflectionModel.use(*rayptr, *hitptr, geometry, rng, rngState);
      rticontextptr->geoFirstHitTFar = rayptr->tfar;
      rticontextptr->geoNotIntersected = false;
      // In a triangle context, always accept hits.
      return;
    }

    static
    void filter_fun_boundary(RTCFilterFunctionNArguments const* args)
    {
      RLOG_DEBUG << "filter_fun_boundary()" << std::endl;
      assert(args->N == 1 && "Precondition");
      // This function gets a pointer to a context object in args.context
      auto cc = args->context;

      //std::cerr << "filter_fun_boundary(): the address of the rtc context: " << cc << std::endl;
      auto ccnonconst = const_cast<RTCIntersectContext*>(cc);
      auto rtiabscontextptr = &reinterpret_cast<typename rti::trace::absc_context<numeric_type>::context_c_wrapper*> (ccnonconst)->mAbscContext;
      auto rticontextptr = reinterpret_cast<rti::trace::triangle_context<numeric_type>*> (rtiabscontextptr);
      //std::cerr << "filter_fun_boundary(): address of the rti context: " << rticontextptr << std::endl;

      // The rticontextptr now serves an equal function as the this pointer in a conventional
      // (non-static) member function.
      assert(args->N == 1 && "Precondition: for the cast");
      auto  rayptr = reinterpret_cast<RTCRay*> (args->ray); // ATTENTION: this cast works only if N == 1 holds
      auto  hitptr = reinterpret_cast<RTCHit*> (args->hit); // ATTENTION: this cast works only if N == 1 holds
      auto  validptr = args->valid;
      // Reference all the data in the context with local variables
      auto& boundaryID = rticontextptr->mBoundaryID;
      auto& rng =        rticontextptr->rng;
      auto& rngState =   rticontextptr->rngstate;
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
      boundRayout = rticontextptr->mBoundaryReflectionModel.use(*rayptr, *hitptr, boundary, rng, rngState);
      rticontextptr->boundFirstHitTFar = rayptr->tfar;
      rticontextptr->boundNotIntersected = false;
      return;
    }

  private:

    void post_process_intersect(RTCRayHit& pRayHit)
    {
      RLOG_DEBUG
        << "post_process_intersect(); mGeoHitPrimIDs.size() == " << this->mGeoHitPrimIDs.size() << "; ";
      for(auto const& ee : this->mGeoHitPrimIDs)
        RLOG_DEBUG << ee << " ";
      RLOG_DEBUG
        << std::endl;
      RLOG_DEBUG << "this->rayWeight == " << this->rayWeight << std::endl;

      if (this->geoNotIntersected && this->boundNotIntersected) {
        // The ray did not hit anything
        this->reflect = false;
        this->tfar = pRayHit.ray.tfar;
        return;
      }

      this->reflect = true;

      if ( ( ! this->boundNotIntersected && this->geoNotIntersected) ||
           ( ! this->boundNotIntersected && this->boundFirstHitTFar < this->geoFirstHitTFar)) {
        // /the geometry has not been hit. Nevertheless, the boundary has been hit
        this->rayout = this->boundRayout;
        this->tfar = this->boundFirstHitTFar;
        return;
      }
      // A hit on the geometry
      assert ( ! this->geoNotIntersected && "Assumption");
      this->rayout = this->geoRayout;
      this->tfar = this->geoFirstHitTFar;
      auto& hitprimId = pRayHit.hit.primID;
      auto valuetodrop = this->rayWeight * this->mGeometry.get_sticking_coefficient();
      this->mHitAccumulator.use(hitprimId, valuetodrop);
      this->rayWeight -= valuetodrop;

      this->rejection_control_check_weight_reweight_or_kill();
    }

  public:

    void intersect1(RTCScene& pScene, RTCRayHit& pRayHit) override final
    {
      // prepare
      pRayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
      pRayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
      this->geoNotIntersected = true;
      this->boundNotIntersected = true;
      // Embree intersect
      // This call uses our intersection functions which must have been registered with
      // register_intersect_filter_funs() .

      //std::cerr << "intersect1(): the address of the rtc context: "
      //          << &this->mContextCWrapper.mRtcContext << std::endl;
      //std::cerr << "intersect1(): address of the rti context: " << this << std::endl;

      // performing ray queries in a scene is thread-safe
      rtcIntersect1(pScene, &this->mContextCWrapper.mRtcContext, &pRayHit);

      // post process
      this->post_process_intersect(pRayHit);
    }

    void init() override final
    {
      // // RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT flag uses an optimized traversal
      // // algorithm for incoherent rays (default)
      rtcInitIntersectContext(&this->mContextCWrapper.mRtcContext);
      // this->mRtcContext.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
      // this->mRtcContext.filter = nullptr;
      // assert(RTC_MAX_INSTANCE_LEVEL_COUNT == 1 && "Assumption");
      // this->mRtcContext.instID[0] = RTC_INVALID_GEOMETRY_ID; // initialize to some value
    }

    bool compute_exposed_areas_by_sampling() override final
    {
      return false;
    }
  };
}}
