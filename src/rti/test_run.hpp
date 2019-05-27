#pragma once

#include <boost/core/demangle.hpp>
#include <boost/timer/timer.hpp>

#include <cmath>
#include <chrono>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/tbb.h"

#include "rti/i_geometry.hpp"
#include "rti/i_ray_source.hpp"
#include "rti/test_result.hpp"

namespace rti {
  class test_run {
  public:
    test_run(i_geometry& pGeo, i_ray_source& pRaySource) :
      mGeo(pGeo),
      // mRaySource is not used in the moment!
      mRaySource(pRaySource) {};

    rti::test_result run() {
      test_result result {};
      result.inputFilePath = mGeo.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(mGeo).name());

      RTCDevice device = mGeo.get_rtc_device();
      RTCScene scene = rtcNewScene(device);

      // No scene flags for now.
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      // Selecting higher build quality results in better rendering performance but slower
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      RTCBuildQuality bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);

      RTCGeometry geometry = mGeo.get_rtc_geometry();
      rtcCommitGeometry(geometry);

      (void) rtcAttachGeometry(scene, geometry);
      rtcCommitScene(scene);
      // Release geomtery.
      // The geometry object will be destructed when the scene object is
      // destructed, because the geometry object is attached to the scene
      // object.
      rtcReleaseGeometry(geometry);

      //////////////////////////////
      // Ray queries
      //////////////////////////////
      // The Embree tutorial uses one context per ray.
      // I think that is not necessary, though.
      RTCIntersectContext context;
      rtcInitIntersectContext(&context);

      // size_t nrexp = 26;
      // size_t nrexp = 24;
      size_t nrexp = 26;
      size_t numRays = std::pow(2,nrexp);
      result.numRays = numRays; // Save the number of rays also to the test result

      // std::cout << "Running 2**" << nrexp << " ray traces on "
      //           << boost::core::demangle(typeid(mGeo).name()) << std::endl << std::flush;

      const std::vector<rti::pair<float>> ddxys
      {{-3.f, -3.f}, {-1.f, -3.f}, { 1.f, -3.f}, {3.f, -3.f},
       {-3.f, -1.f}, {-1.f, -1.f}, { 1.f, -1.f}, {3.f, -1.f},
       {-3.f,  1.f}, {-1.f,  1.f}, { 1.f,  1.f}, {3.f,  1.f},
       {-3.f,  3.f}, {-1.f,  3.f}, { 1.f,  3.f}, {3.f,  3.f}};

      // std::vector<RTCRay> raypack(ddxys.size());
      // float magic = 0.5f; // magic number

      // for (size_t idx = 0; idx < ddxys.size(); ++idx) {
      //   // Origin:
      //   raypack[idx].org_x = 0.1;
      //   raypack[idx].org_y = 0;
      //   raypack[idx].org_z = 0;
      //   raypack[idx].dir_x = magic;
      //   raypack[idx].dir_y = ddxys[idx].frst;
      //   raypack[idx].dir_z = ddxys[idx].scnd;
      //   // start of ray
      //   raypack[idx].tnear = 0;
      //   // Maximum length of ray
      //   raypack[idx].tfar = std::numeric_limits<float>::max();
      //   raypack[idx].flags = 0;
      // }

      // Thread local variables
      tbb::enumerable_thread_specific<unsigned int> randGrp(1234); // some number
      tbb::enumerable_thread_specific<RTCRayHit> rayhitGrp(RTCRayHit {RTCRay {}, RTCHit {}});
      tbb::enumerable_thread_specific<size_t> rpidxGrp(0);

      // Initialize rays
      for (auto& rayhit : rayhitGrp) {
        rayhit.ray.org_x = 0.1;
        rayhit.ray.org_y = 0;
        rayhit.ray.org_z = 0;
        rayhit.ray.dir_x = 0.5f; // magic number
        rayhit.ray.dir_y = 0;
        rayhit.ray.dir_z = 0;
        rayhit.ray.tnear = 0;
        rayhit.ray.tfar = std::numeric_limits<float>::max();
        rayhit.ray.flags = 0;
      }

      // Start timing until end of variable scope.
      //boost::timer::auto_cpu_timer timerOld;

      boost::timer::cpu_timer timer; // also calls start

      result.startTime = std::chrono::high_resolution_clock::now();

      tbb::parallel_for(
        tbb::blocked_range<size_t>(0, numRays),
        //[raypack = raypack] // C++14; use generalized lambda capture to copy raypack
        //[&] // Capture all by refernce can easily cause errors
        [&scene, &context, &ddxys, &randGrp, &rayhitGrp, &rpidxGrp] // capture by refernce
        (const tbb::blocked_range<size_t>& range) {
          // tbb::enumerable_thread_specific<unsigned int>::reference rand = randGrp.local();
          // tbb::enumerable_thread_specific<RTCRayHit>::reference rayhit = rayhitGrp.local();
          // tbb::enumerable_thread_specific<size_t>::reference rpidx = rpidxGrp.local();
          auto rand = randGrp.local();
          auto rayhit = rayhitGrp.local();
          auto rpidx = rpidxGrp.local();

          // The range often contains only one single element. That is, the tbb::blocked_range() function
          // really is the loop.
          for(size_t idx = range.begin(); idx < range.end(); ++idx) {

            if (rpidx >= ddxys.size()) {
              rpidx = 0;
              rand = fastrand(rand);
            }

            rayhit.ray.dir_x = static_cast<float>(rand);
            rayhit.ray.dir_y = ddxys[rpidx].frst;
            rayhit.ray.dir_z = ddxys[rpidx].scnd;

            rayhit.ray.tfar = std::numeric_limits<float>::max();
            rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

            rtcIntersect1(scene, &context, &rayhit);

            // if (rayhit.ray.tfar * raypack[rpidx].dir_x > 1) {
            // if (rayhit.hit.primID == 0) {
            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
              std::cerr << ">>>>> ERROR: ray did not hit anything"
                        << " rayhit.ray.tfar == " << rayhit.ray.tfar 
                        << " raypack[rpidx].dir_x == " << rayhit.ray.dir_x
                        << " raypack[rpidx].dir_y == " << rayhit.ray.dir_y
                        << " raypack[rpidx].dir_z == " << rayhit.ray.dir_z
                        << " primID == " << rayhit.hit.primID << std::endl;
            } else {
              // std::cerr << ">>>>>> no error in this iteration"
              //           << " rayhit.ray.tfar == " << rayhit.ray.tfar 
              //           << " raypack[rpidx].dir_x == " << raypack[rpidx].dir_x
              //           << " raypack[rpidx].dir_y == " << raypack[rpidx].dir_y
              //           << " raypack[rpidx].dir_z == " << raypack[rpidx].dir_z
              //           << " primID == " << rayhit.hit.primID << std::endl;
            }
          }
        });

        result.endTime = std::chrono::high_resolution_clock::now();

        result.timeNanoseconds = timer.elapsed().wall; // Using boost timer

        return result;
    }
  private:
    i_geometry& mGeo;
    i_ray_source& mRaySource;

    //unsigned int g_seed = 1234;
    inline static
    int fastrand(unsigned int pSeed) {
      pSeed = (214013*pSeed+2531011);
      return (pSeed>>16)&0x7FFF; // bs
    }
  };
}
