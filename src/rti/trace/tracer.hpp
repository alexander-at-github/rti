#pragma once

#include <boost/core/demangle.hpp>

#include <cmath>
#include <chrono>
#include <omp.h>

#include "rti/geo/i_boundary.hpp"
#include "rti/geo/i_geometry.hpp"
#include "rti/ray/i_source.hpp"
#include "rti/reflection/diffuse.hpp"
#include "rti/reflection/specular.hpp"
#include "rti/trace/counter.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/result.hpp"
#include "rti/util/timer.hpp"

namespace rti { namespace trace {
  template<typename Ty>
  class tracer {
  public:

    tracer(rti::geo::i_geometry<Ty>& pGeo, rti::geo::i_boundary<Ty>& pBoundary, rti::ray::i_source& pSource) :
      mGeo(pGeo),
      mBoundary(pBoundary),
      mSource(pSource) {
      assert (pGeo.get_rtc_device() == pBoundary.get_rtc_device() &&
        "the geometry and the boundary need to refer to the same Embree (rtc) device");
    }

    rti::trace::result run() {
      // Prepare a data structure for the result.
      auto result = rti::trace::result {};
      result.inputFilePath = mGeo.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(mGeo).name());

      // Prepare Embree
      auto device = mGeo.get_rtc_device();
      auto scene = rtcNewScene(device);
      // No scene flags.
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      // Selecting higher build quality results in better rendering performance but slower3
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      auto bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);
      auto geometry = mGeo.get_rtc_geometry();
      // rtcCommitGeometry(geometry); // Removed; should be done in the implementations of i_geometry
      auto boundary = mBoundary.get_rtc_geometry();

      // rtcAttachGeometry() is thread safe
      auto geometryID = rtcAttachGeometry(scene, geometry);
      auto boundaryID = rtcAttachGeometry(scene, boundary);

      // Use openMP for parallelization
      #pragma omp parallel
      {
        rtcJoinCommitScene(scene);
        // TODO: move to the other parallel region at the bottom
      }

      // The geometry object will be destructed when the scene object is
      // destructed, because the geometry object is attached to the scene
      // object.
      rtcReleaseGeometry(geometry);
      //
      auto context = RTCIntersectContext {};
      rtcInitIntersectContext(&context);

      // *Ray queries*
      //size_t nrexp = 27;
      auto nrexp = 20; // int
      auto numRays = std::pow(2.0, nrexp); // returns a double
      result.numRays = numRays; // Save the number of rays also to the test result

      //rti::reflection::specular reflectionModel;
      //rti::reflection::diffuse reflectionModel(0.015625);
      auto reflectionModel = rti::reflection::diffuse<Ty> {0.01};
      auto boundaryReflection = rti::reflection::specular<Ty> {};

      auto geohitc = 0ull; // unsigned long long int
      auto nongeohitc = 0ull;
      auto hitCounter = rti::trace::counter {mGeo.get_num_primitives()};
      //rti::trace::counter hitCounter; // TODO: switch
      #pragma omp declare \
        reduction(hit_counter_combine : \
                  rti::trace::counter : \
                  omp_out = rti::trace::counter(omp_out, omp_in)) \
        initializer(omp_priv = rti::trace::counter(omp_orig))
        //initializer(omp_priv = rti::trace::counter(omp_orig))

      // omp_set_dynamic(false);

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      auto rng = std::make_unique<rti::rng::cstdlib_rng>();

      // Start timing
      auto timer = rti::util::timer {};

      #pragma omp parallel \
        reduction(+ : geohitc, nongeohitc) \
        reduction(hit_counter_combine : hitCounter)
      {
        // Thread local data goes here, if it is not needed anymore after the
        // execution of the parallel region.
        auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        // REMARK: All data which is modified in the parallel loop should be
        // handled here explicitely.
        auto seed = (unsigned int) ((omp_get_thread_num() + 1) * 29); // multiply by magic number (prime)
        // It seems really important to use two separate seeds / states for
        // sampling the source and sampling reflections. When we use only one
        // state for both, then the variance is very high.
        auto rngSeed1 = std::make_unique<rti::rng::cstdlib_rng::state>(seed);
        auto rngSeed2 = std::make_unique<rti::rng::cstdlib_rng::state>(seed+2);
        // TODO: move this initialization to, e.g., the constructor

        // A dummy counter for the boundary
        auto boundaryCntr = rti::trace::dummy_counter {};

        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) numRays; ++idx) {

          rayhit.ray.tnear = 0;
          rayhit.ray.time = 0;
          rayhit.ray.tfar = std::numeric_limits<float>::max();

          RLOG_DEBUG << "Preparing new ray from source" << std::endl;
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed1);

          bool reflect;
          do {
            // { // Debug
            //   #pragma omp critical
            //   {
            //     std::cout << "thread " << omp_get_thread_num()
            //               << " rng-state " << rngSeed->mSeed << std::endl;
            //   }
            // } // Debug

            rayhit.ray.tfar = std::numeric_limits<float>::max();
            rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

            RLOG_DEBUG
              << "preparing ray == ("
              << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z
              << ") ("
              << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z
              << ")" << std::endl;
            // performing ray queries in a scene is thread-safe
            rtcIntersect1(scene, &context, &rayhit);

            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
              // No  hit
              nongeohitc += 1;
              break; // break do-while loop
            }
            // else
            // A hit
            if (rayhit.hit.geomID == boundaryID) {
              // Ray hit the boundary
              boundaryReflection.use(rayhit, *rng, *rngSeed2, this->mBoundary, boundaryCntr);
              continue;
            }
            geohitc += 1;

            RLOG_DEBUG << "rayhit.hit.primID == " << rayhit.hit.primID << std::endl;
            RLOG_DEBUG << "prim == " << mGeo.prim_to_string(rayhit.hit.primID) << std::endl;

            reflect = reflectionModel.use(rayhit, *rng, *rngSeed2, this->mGeo, hitCounter);
          } while (reflect);
        }
      }

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitCounter = std::make_unique<rti::trace::counter>(hitCounter);
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      // // Turn dynamic adjustment of number of threads on again
      // omp_set_dynamic(true);

      return result;
    }
  private:
    rti::geo::i_geometry<Ty>& mGeo;
    rti::geo::i_boundary<Ty>& mBoundary;
    rti::ray::i_source& mSource;
  };
}} // namespace
