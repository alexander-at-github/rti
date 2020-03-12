#pragma once

#include <vector>
#include <set>

#include <embree3/rtcore.h>

#include "rti/mc/rejection_control.hpp"
#include "rti/particle/i_particle.hpp"
#include "rti/reflection/i_reflection_model.hpp"
#include "rti/reflection/specular.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/absc_context.hpp"
#include "rti/trace/i_hit_accumulator.hpp"

// This class needs to be used according to a protocol! See base class rti::trace::absc_context

namespace rti { namespace trace {
  template<typename numeric_type> // intended to be a numeric type
  class point_cloud_context_simplified : public rti::trace::absc_context<numeric_type> {

  public:
    // managment for ray weights
    static constexpr float INITIAL_RAY_WEIGHT = 1.0f;
    // =================================================================
    // CHOOSING A GOOD VALUE FOR THE WEIGHT LOWER THRESHOLD IS IMPORTANT
    // =================================================================
    static constexpr float RAY_WEIGHT_LOWER_THRESHOLD = 0.005f;
    static constexpr float RAY_RENEW_WEIGHT = 20 * RAY_WEIGHT_LOWER_THRESHOLD; // magic number

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

    rti::rng::i_rng& mRng;
    rti::rng::i_rng::i_state& mRngState;

    rti::mc::rejection_control<numeric_type> rejectioncontrol;
    rti::particle::i_particle<numeric_type>& particle;


    // A vector of primitive IDs collected through the filter function filter_fun_geometry() which we
    // will then post process in the post_process_intersection() function.
    // Initialize to some reasonable size. The vector may grow, if needed.
    // Weird initialization necessary
    std::vector<unsigned int> mGeoHitPrimIDs {};

  public:
    point_cloud_context_simplified(unsigned int pGeometryID,
            rti::geo::absc_point_cloud_geometry<numeric_type>& pGeometry,
            rti::reflection::i_reflection_model<numeric_type>& pReflectionModel,
            rti::trace::i_hit_accumulator<numeric_type>& pHitAccumulator,
            unsigned int pBoundaryID,
            rti::geo::i_boundary<numeric_type>& pBoundary,
            rti::reflection::i_reflection_model<numeric_type>& pBoundaryReflectionModel,
            rti::rng::i_rng& pRng,
            rti::rng::i_rng::i_state& pRngState,
            rti::particle::i_particle<numeric_type>& particle) :
      rti::trace::absc_context<numeric_type>(INITIAL_RAY_WEIGHT, false, geoRayout, 0),// initialize to some values
      mGeometryID(pGeometryID),
      mGeometry(pGeometry),
      mReflectionModel(pReflectionModel),
      mHitAccumulator(pHitAccumulator),
      mBoundaryID(pBoundaryID),
      mBoundary(pBoundary),
      mBoundaryReflectionModel(pBoundaryReflectionModel),
      mRng(pRng),
      mRngState(pRngState),
      rejectioncontrol(mRng, mRngState),
      particle(particle) {
      mGeoHitPrimIDs.reserve(32); // magic number // Reserve a reasonable number of hit elements for one ray
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
      // Don't do anything
    }

  public:
    void intersect1(RTCScene& pScene, RTCRayHit& pRayHit) override final
    {
      // prepare
      pRayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
      pRayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

      // Embree intersect
      // performing ray queries in a scene is thread-safe
      rtcIntersect1(pScene, &(this->mContextCWrapper.mRtcContext), &pRayHit);

      this->tfar = pRayHit.ray.tfar;

      if (pRayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
        // no hit
        this->reflect = false;
        RLOG_TRACE << "-h";
        return;
      }
      // else
      // a hit
      if (pRayHit.hit.geomID == this->mBoundaryID) {
        RLOG_DEBUG << "Ray hit boundary" << std::endl;
        this->boundRayout = this->mBoundaryReflectionModel.use(
          pRayHit.ray, pRayHit.hit, this->mBoundary, this->mRng, this->mRngState);
        this->rayout = this->boundRayout;
        this->reflect = true;
        RLOG_TRACE << "b";
        return;
      }
      // else
      if (pRayHit.hit.geomID == this->mGeometryID) {
        // ray hit the geometry
        RLOG_DEBUG << "Ray hit primitive with ID " << pRayHit.hit.primID << std::endl;

        // If the dot product of the ray direction and the surface normal is greater than zero, then
        // we hit the back face of the disc.
        auto const& ray = pRayHit.ray;
        auto const& hit = pRayHit.hit;
        if (rti::util::dot_product(rti::util::triple<numeric_type> {ray.dir_x, ray.dir_y, ray.dir_z},
                                   this->mGeometry.get_normal(hit.primID)) > 0) {
          // Hit from the back
          this->reflect = false;
          RLOG_TRACE << "a";
          return;
          // TODO: One could consider to let hits with very small values in tfar through.
        }

        auto sticking = particle.process_hit(pRayHit.hit.primID);
        this->geoRayout = this->mReflectionModel.use(
          pRayHit.ray, pRayHit.hit, this->mGeometry, this->mRng, this->mRngState);
        this->rayout = this->geoRayout;
        auto valuetodrop = this->rayWeight * sticking;
        this->mHitAccumulator.use(pRayHit.hit.primID, valuetodrop);
        this->rayWeight -= valuetodrop;
      } else {
        assert(false && "Assumption");
      }

      rejectioncontrol.check_weight_reweight_or_kill
        (*this, this->RAY_WEIGHT_LOWER_THRESHOLD, this->RAY_RENEW_WEIGHT);
      if (this->reflect) {
        RLOG_TRACE << "r";
      } else {
        RLOG_TRACE << "-r";
      }
    }

    void init() override final
    {
      // RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT flag uses an optimized traversal
      // algorithm for incoherent rays (default).
      // this->mRtcContext.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
      rtcInitIntersectContext(&this->mContextCWrapper.mRtcContext);
    }

    void init_ray_weight() override final
    {
      this->rayWeight = this->INITIAL_RAY_WEIGHT;
      RLOG_TRACE << "I";
    }
  };

}}
