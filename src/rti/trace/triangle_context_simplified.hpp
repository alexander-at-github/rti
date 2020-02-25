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
  template<typename Ty> // intended to be a numeric type
  class triangle_context_simplified : public rti::trace::absc_context<Ty> {
  private:
    // managment for ray weights
    static constexpr float INITIAL_RAY_WEIGHT = 1.0f;
    // =================================================================
    // CHOOSING A GOOD VALUE FOR THE WEIGHT LOWER THRESHOLD IS IMPORTANT
    // =================================================================
    static constexpr float RAY_WEIGHT_LOWER_THRESHOLD = 0.1f;
    static constexpr float RAY_RENEW_WEIGHT = 3 * RAY_WEIGHT_LOWER_THRESHOLD; // magic number
  private:
    // geometry related data
    rti::util::pair<rti::util::triple<Ty> > geoRayout {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // boundary related data
    rti::util::pair<rti::util::triple<Ty> > boundRayout {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // other data
  private:
    unsigned int mGeometryID = RTC_INVALID_GEOMETRY_ID; // initialize to some useful value
    rti::geo::triangle_geometry<Ty>& mGeometry;
    rti::reflection::i_reflection_model<Ty>& mReflectionModel;
    rti::trace::i_hit_accumulator<Ty>& mHitAccumulator;
    unsigned int mBoundaryID = RTC_INVALID_GEOMETRY_ID; // initialize to some useful value
    rti::geo::i_boundary<Ty>& mBoundary;
    rti::reflection::i_reflection_model<Ty>& mBoundaryReflectionModel;

    rti::rng::i_rng& mRng;
    rti::rng::i_rng::i_state& mRngState;

    // A vector of primitive IDs collected through the filter function filter_fun_geometry() which we
    // will then post process in the post_process_intersection() function.
    // Initialize to some reasonable size. The vector may grow, if needed.
    // Weird initialization necessary
    std::vector<unsigned int> mGeoHitPrimIDs {};

    // Constructor
  public:
    triangle_context_simplified(unsigned int pGeometryID,
            rti::geo::triangle_geometry<Ty>& pGeometry,
            rti::reflection::i_reflection_model<Ty>& pReflectionModel,
            rti::trace::i_hit_accumulator<Ty>& pHitAccumulator,
            unsigned int pBoundaryID,
            rti::geo::i_boundary<Ty>& pBoundary,
            rti::reflection::i_reflection_model<Ty>& pBoundaryReflectionModel,
            rti::rng::i_rng& pRng,
            rti::rng::i_rng::i_state& pRngState) :
      // initialize members of virtual base class
      rti::trace::absc_context<Ty>(INITIAL_RAY_WEIGHT, false, geoRayout, 0), // initialize to some values
      mGeometryID(pGeometryID),
      mGeometry(pGeometry),
      mReflectionModel(pReflectionModel),
      mHitAccumulator(pHitAccumulator),
      mBoundaryID(pBoundaryID),
      mBoundary(pBoundary),
      mBoundaryReflectionModel(pBoundaryReflectionModel),
      mRng(pRng),
      mRngState(pRngState) {
      mGeoHitPrimIDs.reserve(32); // magic number // Reserve some reasonable number of hit elements for one ray
      mGeoHitPrimIDs.clear();
    }

  public:
    static
    void register_intersect_filter_funs(rti::geo::i_geometry<Ty>& pGeometry,
                                        rti::geo::i_boundary<Ty>& pBoundary) {
      // The following cast characterizes a precondition to this function
      auto pPCGeoPointer = dynamic_cast<rti::geo::triangle_geometry<Ty>*> (&pGeometry);
      auto& pPCGeo = *pPCGeoPointer;
      RLOG_DEBUG << "register_intersect_filter_funs()" << std::endl;
      // Don't do anything.
    }

  public:
    void intersect1(RTCScene& pScene, RTCRayHit& pRayHit) override final {
      // In this class, triangle_context_intersect, the intersection is a straight forward
      // Embree intersection.
      // prepare
      pRayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
      pRayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

      // Embree intersect
      // performing ray queries in a scene is thread-safe
      rtcIntersect1(pScene, &this->mContextCWrapper.mRtcContext, &pRayHit);

      this->tfar = pRayHit.ray.tfar;

      if (pRayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
        // no hit
        this->reflect = false;
        return;
      }
      // else
      // a hit
      if (pRayHit.hit.geomID == this->mBoundaryID) {
        // ray hit the boundary
        this->boundRayout = this->mBoundaryReflectionModel.use(
          pRayHit.ray, pRayHit.hit, this->mBoundary, this->mRng, this->mRngState);
        this->rayout = this->boundRayout;
        this->reflect = true;
        return;
      } else if (pRayHit.hit.geomID == this->mGeometryID) {
        // ray hit the geometry
        this->geoRayout = this->mReflectionModel.use(
          pRayHit.ray, pRayHit.hit, this->mGeometry, this->mRng, this->mRngState);
        this->rayout = this->geoRayout;
        auto valuetodrop = this->rayWeight * this->mGeometry.get_sticking_coefficient();
        this->mHitAccumulator.use(pRayHit.hit.primID, valuetodrop);
        this->rayWeight -= valuetodrop;
      } else {
        assert(false && "Assumption");
      }

      weight_check_reweight_kill();
    }
  private:
    void weight_check_reweight_kill() {
      // We do what is sometimes called Roulette in MC literatur.
      // Jun Liu calls it "rejection controll" in his book.
      // If the weight of the ray is above a certain threshold, we always reflect.
      // If the weight of the ray is below the threshold, we randomly decide to either kill the
      // ray or increase its weight (in an unbiased way).
      this->reflect = true;
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
    void init() override final {
      // // RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT flag uses an optimized traversal
      // // algorithm for incoherent rays (default)
      rtcInitIntersectContext(&this->mContextCWrapper.mRtcContext);
      // this->mRtcContext.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
      // this->mRtcContext.filter = nullptr;
      // assert(RTC_MAX_INSTANCE_LEVEL_COUNT == 1 && "Assumption");
      // this->mRtcContext.instID[0] = RTC_INVALID_GEOMETRY_ID; // initialize to some value
    }

    void init_ray_weight() override final {
      this->rayWeight = this->INITIAL_RAY_WEIGHT;
    }
  };
}}
