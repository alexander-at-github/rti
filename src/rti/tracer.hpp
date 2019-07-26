#pragma once

#include <boost/core/demangle.hpp>

#include <cmath>
#include <chrono>
#include <omp.h>

#include "rti/bucket_counter.hpp"
#include "rti/dummy_counter.hpp"
#include "rti/i_boundary.hpp"
#include "rti/i_geometry.hpp"
#include "rti/i_ray_source.hpp"
#include "rti/diffuse_reflection.hpp"
#include "rti/specular_reflection.hpp"
#include "rti/timer.hpp"
#include "rti/trace_result.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp" // only debug

namespace rti {
  template<typename Ty>
  class tracer {
  public:

    tracer(i_geometry<Ty>& pGeo, i_boundary<Ty>& pBoundary, i_ray_source& pSource) :
      mGeo(pGeo),
      mBoundary(pBoundary),
      mSource(pSource) {
      assert (pGeo.get_rtc_device() == pBoundary.get_rtc_device() &&
        "the geometry and the boundary need to refer to the same Embree (rtc) device");
    }

    rti::trace_result run() {
      // Prepare a data structure for the result.
      auto result = trace_result {};
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

      //rti::specular_reflection reflectionModel;
      //rti::diffuse_reflection reflectionModel(0.015625);
      auto reflectionModel = rti::diffuse_reflection<Ty> {0.01};
      auto boundaryReflection = rti::specular_reflection<Ty> {};

      auto geohitc = 0ull; // unsigned long long int
      auto nongeohitc = 0ull;
      auto bucketCounter = rti::bucket_counter {45, 101};
      #pragma omp declare \
        reduction(bucket_counter_combine : \
                  rti::bucket_counter : \
                  omp_out = rti::bucket_counter::combine(omp_out, omp_in)) \
        initializer(omp_priv = rti::bucket_counter(omp_orig))
        //initializer(omp_priv = bucketCounter(45, 101))
      // TODO: IS THAT CORRECT? (see also the copy constructor of bucket_counter)

      // omp_set_dynamic(false);

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      auto rng = std::make_unique<rti::cstdlib_rng>();

      // Start timing
      auto timer = rti::timer {};

      #pragma omp parallel \
        reduction(+ : geohitc, nongeohitc) \
        reduction(bucket_counter_combine : bucketCounter)
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
        auto rngSeed1 = std::make_unique<rti::cstdlib_rng::state>(seed);
        auto rngSeed2 = std::make_unique<rti::cstdlib_rng::state>(seed+2);
        // TODO: move this initialization to, e.g., the constructor

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
              // TODO
            }
            geohitc += 1;

            RLOG_DEBUG << "rayhit.hit.primID == " << rayhit.hit.primID << std::endl;
            RLOG_DEBUG << "prim == " << mGeo.prim_to_string(rayhit.hit.primID) << std::endl;

            reflect = reflectionModel.use(rayhit, *rng, *rngSeed2, this->mGeo, bucketCounter);
          } while (reflect);
        }
      }

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      // // Turn dynamic adjustment of number of threads on again
      // omp_set_dynamic(true);

      //std::cout << bucketCounter << std::endl;

      return result;
    }
  private:
    i_geometry<Ty>& mGeo;
    i_boundary<Ty>& mBoundary;
    i_ray_source& mSource;
  };
}
