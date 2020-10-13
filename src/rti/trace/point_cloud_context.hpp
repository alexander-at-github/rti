#pragma once

#include <vector>
#include <set>

#include <embree3/rtcore.h>

#include "../reflection/i_reflection_model.hpp"
#include "../reflection/specular.hpp"
#include "absc_context.hpp"
#include "dummy_counter.hpp"
#include "i_hit_accumulator.hpp"

// This class needs to be used according to a protocol! See base class rti::trace::absc_context

namespace rti { namespace trace {
  template<typename numeric_type>
  class point_cloud_context : public rti::trace::absc_context<numeric_type> {

  private:
    // geometry related data
    bool geoNotIntersected = true;
    numeric_type geoFirstHitTFar = 0;
    numeric_type geoTFarMax = 0;
    rti::util::pair<rti::util::triple<numeric_type> > geoRayout {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // boundary related data
    bool boundNotIntersected = true;
    numeric_type boundFirstHitTFar = 0;
    rti::util::pair<rti::util::triple<numeric_type> > boundRayout {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // other data
  private:
    unsigned int mGeometryID = RTC_INVALID_GEOMETRY_ID;
    rti::geo::absc_point_cloud_geometry<numeric_type>& mGeometry;
    rti::reflection::i_reflection_model<numeric_type>& mReflectionModel;
    rti::trace::i_hit_accumulator<numeric_type>& mHitAccumulator;
    unsigned int mBoundaryID = RTC_INVALID_GEOMETRY_ID;
    rti::geo::i_boundary<numeric_type>& mBoundary;
    rti::reflection::i_reflection_model<numeric_type>& mBoundaryReflectionModel;

    rti::particle::i_particle<numeric_type>& particle;

    // A vector of primitive IDs collected through the filter function filter_fun_geometry() which we
    // will then post process in the post_process_intersection() function.
    // Initialize to some reasonable size. The vector may grow, if needed.
    // Weird initialization necessary
    std::vector<unsigned int> mGeoHitPrimIDs {};


  public:
    point_cloud_context(unsigned int pGeometryID,
            rti::geo::absc_point_cloud_geometry<numeric_type>& pGeometry,
            rti::reflection::i_reflection_model<numeric_type>& pReflectionModel,
            rti::trace::i_hit_accumulator<numeric_type>& pHitAccumulator,
            unsigned int pBoundaryID,
            rti::geo::i_boundary<numeric_type>& pBoundary,
            rti::reflection::i_reflection_model<numeric_type>& pBoundaryReflectionModel,
            rti::rng::i_rng& pRng,
            rti::rng::i_rng::i_state& pRngState,
            rti::particle::i_particle<numeric_type>& particle) :
      rti::trace::absc_context<numeric_type>(false, geoRayout, 0, pRng, pRngState), // initialize to some values
      mGeometryID(pGeometryID),
      mGeometry(pGeometry),
      mReflectionModel(pReflectionModel),
      mHitAccumulator(pHitAccumulator),
      mBoundaryID(pBoundaryID),
      mBoundary(pBoundary),
      mBoundaryReflectionModel(pBoundaryReflectionModel),
      particle(particle) {
      mGeoHitPrimIDs.reserve(32); // magic number // Reserve some reasonable number of hit elements for one ray
      mGeoHitPrimIDs.clear();
    }

  public:
    static
    void register_intersect_filter_funs(rti::geo::i_geometry<numeric_type>& pGeometry,
                                        rti::geo::i_boundary<numeric_type>& pBoundary)
    {
      // The following cast characterizes a precondition to this function
      auto pPCGeoPointer = dynamic_cast<rti::geo::absc_point_cloud_geometry<numeric_type>*> (&pGeometry);
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

      auto ccnonconst = const_cast<RTCIntersectContext*>(cc);
      auto rtiabscontextptr = &reinterpret_cast
        <typename rti::trace::absc_context<numeric_type>::context_c_wrapper*> (ccnonconst)->mAbscContext;
      auto rticontextptr = reinterpret_cast<rti::trace::point_cloud_context<numeric_type>*> (rtiabscontextptr);

      // The rticontextptr now serves an equal function as the this pointer in a conventional (non-static) member function.
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

      // // Debug
      // if(rticontextptr->geoNotIntersected) {
      //   RLOG_DEBUG << "filter_fun_geometry(): FIRST-HIT IS SET ";
      // } else {
      //   RLOG_DEBUG << "filter_fun_geometry(): FIRST-HIT IS *NOT* SET ";
      // }
      // RLOG_DEBUG << "hitptr->geomID == " << hitptr->geomID << " primID == " << hitptr->primID << std::endl;
      // // Non-Debug
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
      if (rti::util::dot_product(rti::util::triple<numeric_type> {rayptr->dir_x, rayptr->dir_y, rayptr->dir_z},
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
        // Set epsilon equal to 0.5 times the radius of the first disc which is hit by the ray.
        // the fourth element of the primitive is the radius
        auto epsilon = 0.5 * geometry.get_prim(hitptr->primID)[3];
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
    void filter_fun_boundary(RTCFilterFunctionNArguments const* args)
    {
      RLOG_DEBUG << "filter_fun_boundary()" << std::endl;
      assert(args->N == 1 && "Precondition");
      // This function gets a pointer to a context object in args.context
      auto cc = args->context;

      auto ccnonconst = const_cast<RTCIntersectContext*>(cc);
      auto rtiabscontextptr = &reinterpret_cast
        <typename rti::trace::absc_context<numeric_type>::context_c_wrapper*> (ccnonconst)->mAbscContext;
      auto rticontextptr = reinterpret_cast<rti::trace::point_cloud_context<numeric_type>*> (rtiabscontextptr);

      // The rticontextptr now serves an equal function as the this pointer in a conventional (non-static) member function.
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
      // // Debug
      // if(rticontextptr->boundNotIntersected) {
      //   RLOG_DEBUG << "filter_fun_boundary(): FIRST-HIT IS SET ";
      // } else {
      //   RLOG_DEBUG << "filter_fun_boundary(): FIRST-HIT IS *NOT* SET ";
      // }
      // RLOG_DEBUG << "hitptr->geomID == " << hitptr->geomID << " primID == " << hitptr->primID << std::endl;
      // // Non-Debug
      if(rticontextptr->boundNotIntersected) {
        boundRayout = rticontextptr->mBoundaryReflectionModel.use(*rayptr, *hitptr, boundary, rng, rngState);
        rticontextptr->boundFirstHitTFar = rayptr->tfar;
        rticontextptr->boundNotIntersected = false;
      }

      if ( (! rticontextptr->geoNotIntersected) &&
          rayptr->tfar < rticontextptr->geoTFarMax) {
        args->valid[0] = 0; // continue this ray
      }
      // Accept boundary hit.
      args->valid[0] = -1; // accept / stop this rayptr
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

      // Note that (pRayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID && ! this->mGeoHitPrimIDs.empty())
      // may hold since we do not always clear the vector mGeoHitPrimIDs.
      // The right way to check that the ray did not hit anything is to check this->geoNotIntersected and
      // this->boundNotIntersected.

      if (this->geoNotIntersected && this->boundNotIntersected) {
        // The ray did not hit anything
        this->reflect = false;
        this->tfar = pRayHit.ray.tfar;
        return;
      }
      this->reflect = true;

      if ( ! this->boundNotIntersected ) {
        if ( this->geoNotIntersected ||
             this->boundFirstHitTFar < this->geoFirstHitTFar) {
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
        auto const ray = pRayHit.ray;
        auto sticking = particle.process_hit(pRayHit.hit.primID, {ray.dir_x, ray.dir_y, ray.dir_z});
        assert(0 <= sticking && sticking <= 1 && "Assumption");
        auto valuetodrop = this->rayWeight * sticking;
        this->mHitAccumulator.use(hitprimid, valuetodrop);
        sumvaluedroped += valuetodrop;
      }
      auto hitprimcount = set.size();
      assert(this->rayWeight >= sumvaluedroped / hitprimcount && "Assertion");
      this->rayWeight -= sumvaluedroped / hitprimcount;

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
      rtcIntersect1(pScene, &this->mContextCWrapper.mRtcContext, &pRayHit);
      this->post_process_intersect(pRayHit);
    }

    void init() override final
    {
      // RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT flag uses an optimized traversal
      // algorithm for incoherent rays (default).
      // this->mRtcContext.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
      rtcInitIntersectContext(&this->mContextCWrapper.mRtcContext);
    }

    bool compute_exposed_areas_by_sampling() override final
    {
      return false;
    }
  };
}}
