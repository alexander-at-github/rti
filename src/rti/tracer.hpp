#pragma once

#include <boost/core/demangle.hpp>

#include <cmath>
#include <chrono>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/tbb.h"

#include "rti/bucket_counter.hpp"
#include "rti/i_geometry.hpp"
#include "rti/i_ray_source.hpp"
#include "rti/lambertian_reflection.hpp"
#include "rti/specular_reflection.hpp"
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
      rtcCommitScene(scene);
      // The geometry object will be destructed when the scene object is
      // destructed, because the geometry object is attached to the scene
      // object.
      rtcReleaseGeometry(geometry);
      //
      RTCIntersectContext context;
      rtcInitIntersectContext(&context);

      // *Ray queries*
      //size_t nrexp = 27;
      size_t nrexp = 26;
      size_t numRays = std::pow(2,nrexp);
      result.numRays = numRays; // Save the number of rays also to the test result

      // We use both, tbb::enumerable_thread_specific and c++ thread_local variables (further
      // down).
      tbb::enumerable_thread_specific<RTCRayHit> rayhitGrp(
        RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
      tbb::enumerable_thread_specific<size_t> hitcGrp(0);
      tbb::enumerable_thread_specific<size_t> nonhitcGrp(0);
      //tbb::enumerable_thread_specific<std::vector<rti::triple<float> > > hitpointgrp;

      rti::bucket_counter bucketCounterPrototype(45, 101); // lenght 45 and 100 buckets
      tbb::enumerable_thread_specific<rti::bucket_counter> bucketCounterGrp(bucketCounterPrototype);

      // Initializing rays here does not work, cause TBB creates the members of an
      // enumerable_thread_specific object only when a thread requests it. That is,
      // at this position all enumerable_thread_specific containers are empty.

      //rti::specular_reflection reflectionModel;
      //rti::lambertian_reflection reflectionModel(0.015625);
      rti::lambertian_reflection reflectionModel(0.1);

      // Start timing
      rti::timer timer;

      tbb::parallel_for(
        // magic number: number of elements in one blocked range
        tbb::blocked_range<size_t>(0, numRays, 64),
        // capture by refernce
        [&scene, &hitcGrp, &nonhitcGrp, &rayhitGrp, &bucketCounterGrp,
         &context, &reflectionModel,
         // The only way to capture member variables is to capture the
         // this-reference.
         this]
        (const tbb::blocked_range<size_t>& range) {
          // Here it is important to use the refernce type explicitly, otherwise
          // we cannot modify the content.
          auto& rayhit = rayhitGrp.local();
          auto& hitc = hitcGrp.local();
          auto& nonhitc = nonhitcGrp.local();
          auto& bucketCounter = bucketCounterGrp.local();

          thread_local auto source = mSource.clone();

          rayhit.ray.tnear = 0;
          rayhit.ray.time = 0;
          rayhit.ray.tfar = std::numeric_limits<float>::max();

          // The range may contain only one single element. That is, the
          // tbb::blocked_range() function really is the loop.
          for(size_t idx = range.begin(); idx < range.end(); ++idx) {

            RLOG_DEBUG << "Preparing new ray from source" << std::endl;
            source->fill_ray(rayhit.ray);

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

              reflect = reflectionModel.use(rayhit, this->mGeo, bucketCounter);
            } while (reflect);
          }
        });

      result.timeNanoseconds = timer.elapsed_nanoseconds();

      result.hitc = hitcGrp.combine([](size_t xx, size_t yy){ return xx + yy; });
      result.nonhitc = nonhitcGrp.combine([](size_t xx, size_t yy){ return xx + yy; });
      rti::bucket_counter bucketResult =
        bucketCounterGrp.combine(
          [](rti::bucket_counter b1, rti::bucket_counter b2) {
            return rti::bucket_counter::combine(b1, b2);
          });
      std::cout << bucketResult << std::endl;

      return result;
    }
  private:
    i_geometry& mGeo;
    i_ray_source& mSource;
  };
}
