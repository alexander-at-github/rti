#pragma once

#include <boost/core/demangle.hpp>

#include <cmath>
#include <chrono>
#include <omp.h>

#include "rti/bucket_counter.hpp"
#include "rti/i_geometry.hpp"
#include "rti/i_ray_source.hpp"
#include "rti/diffuse_reflection.hpp"
//#include "rti/specular_reflection.hpp"
#include "rti/timer.hpp"
#include "rti/trace_result.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp" // only debug

namespace rti {
  class tracer {
  public:
    tracer(i_geometry& pGeo, i_ray_source& pSource) :
      mGeo(pGeo),
      mSource(pSource) {}

    rti::trace_result run() {
      // Prepare a data structure for the result.
      trace_result result {};
      result.inputFilePath = mGeo.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(mGeo).name());

      // Prepare Embree
      RTCDevice device = mGeo.get_rtc_device();
      RTCScene scene = rtcNewScene(device);
      // No scene flags.
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      // Selecting higher build quality results in better rendering performance but slower3
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      RTCBuildQuality bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);
      RTCGeometry geometry = mGeo.get_rtc_geometry();
      rtcCommitGeometry(geometry);
      (void) rtcAttachGeometry(scene, geometry);

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
      RTCIntersectContext context;
      rtcInitIntersectContext(&context);

      // *Ray queries*
      //size_t nrexp = 27;
      size_t nrexp = 24;
      size_t numRays = std::pow(2, nrexp);
      result.numRays = numRays; // Save the number of rays also to the test result

      //rti::specular_reflection reflectionModel;
      //rti::diffuse_reflection reflectionModel(0.015625);
      rti::diffuse_reflection reflectionModel(0.1);

      size_t hitc = 0;
      size_t nonhitc = 0;
      rti::bucket_counter bucketCounter(45, 101);
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
      rti::timer timer;

      #pragma omp parallel \
        reduction(+ : hitc, nonhitc) \
        reduction(bucket_counter_combine : bucketCounter)
      {
        // Thread local data goes here, if it is not needed anymore after the
        // execution of the parallel region.
        RTCRayHit rayhit = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        // REMARK: All data which is modified in the parallel loop should be
        // handled here explicitely.
        unsigned int seed = omp_get_thread_num();
        auto rngSeed = std::make_unique<rti::cstdlib_rng::state>(seed);
        // TODO: move this initialization to, e.g., the constructor

        #pragma omp for
        for (size_t idx = 0; idx < numRays; ++idx) {

          rayhit.ray.tnear = 0;
          rayhit.ray.time = 0;
          rayhit.ray.tfar = std::numeric_limits<float>::max();

          RLOG_DEBUG << "Preparing new ray from source" << std::endl;
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed);

          bool reflect;
          do {
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
              nonhitc += 1;
              break; // break do-while loop
            }
            // else
            // A hit
            hitc += 1;

            RLOG_DEBUG << "rayhit.hit.primID == " << rayhit.hit.primID << std::endl;
            RLOG_DEBUG << "prim == " << mGeo.prim_to_string(rayhit.hit.primID) << std::endl;

            reflect = reflectionModel.use(rayhit, *rng, *rngSeed, this->mGeo, bucketCounter);
          } while (reflect);
        }
      }

      result.timeNanoseconds = timer.elapsed_nanoseconds();

      // // Turn dynamic adjustment of number of threads on again
      // omp_set_dynamic(true);

      std::cout << bucketCounter << std::endl;

      return result;
    }
  private:
    i_geometry& mGeo;
    i_ray_source& mSource;
  };
}
