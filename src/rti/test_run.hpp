#pragma once

#include<boost/timer/timer.hpp>

#include <cmath>

#include "tbb/tbb.h"

#include "rti/i_geometry_from_gmsh.hpp"
#include "rti/i_ray_source.hpp"
#include "rti/test_result.hpp"

namespace rti {
  class test_run {
  public:
    test_run(i_geometry_from_gmsh& pGeo, i_ray_source& pRaySource) :
      mGeo(pGeo),
      mRaySource(pRaySource) {};

    test_result run() {
      RTCDevice device = mGeo.get_rtc_device();
      RTCScene scene = rtcNewScene(device);

      // No scene flags for now.
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      // Selecting higher build quality results in better rendering performance but slower
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      RTCBuildQuality bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);

      RTCGeometry geometry = mGeo.get_rtc_geometry();
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Using " << typeid(mGeo).name();
      // BOOST_LOG_SEV(rti::mRLogger, blt::debug) << mGeo.to_string();
      rtcCommitGeometry(geometry);

      (void) rtcAttachGeometry(scene, geometry);
      rtcCommitScene(scene);
      // Release geomtery.
      // The geometry object will be destructed when the scene object is
      // destructed, because the geometry object is attached to the scene
      // objec = 0t
      rtcReleaseGeometry(geometry);

      //////////////////////////////
      // Ray queries
      //////////////////////////////
      // The Embree tutorial uses one context per ray.
      // I think that is not necessary, though.
      RTCIntersectContext context;
      rtcInitIntersectContext(&context);

      size_t nrexp = 26;
      size_t numRays = std::pow(2,nrexp);

      std::cout << "Running 2**" << nrexp << " ray traces on "
                << boost::core::demangle(typeid(mGeo).name()) << std::endl << std::flush;

      const std::vector<std::pair<float, float>> ddxys
      {{-3.f, -3.f}, {-1.f, -3.f}, { 1.f, -3.f}, {3.f, -3.f},
       {-3.f, -1.f}, {-1.f, -1.f}, { 1.f, -1.f}, {3.f, -1.f},
       {-3.f,  1.f}, {-1.f,  1.f}, { 1.f,  1.f}, {3.f,  1.f},
       {-3.f,  3.f}, {-1.f,  3.f}, { 1.f,  3.f}, {3.f,  3.f}};
      std::vector<RTCRay> raypack(ddxys.size());
      float magic = 0.5f; // magic number

      for (size_t idx = 0; idx < ddxys.size(); ++idx) {
        // Origin:
        raypack[idx].org_x = 1;
        raypack[idx].org_y = 0;
        raypack[idx].org_z = 0;
        raypack[idx].dir_x = magic;
        raypack[idx].dir_y = ddxys[idx].first;
        raypack[idx].dir_z = ddxys[idx].second;
        // start of ray
        raypack[idx].tnear = 0;
        // Maximum length of ray
        raypack[idx].tfar = std::numeric_limits<float>::max();
      }

      RTCHit hit{};
      RTCRayHit rayhit {raypack[0], hit}; // initialize arbitrary

      // Start timing until end of variable scope.
      boost::timer::auto_cpu_timer timer;

      tbb::parallel_for(
        tbb::blocked_range<size_t>(0,numRays),
        [&](const tbb::blocked_range<size_t>& range) {
          unsigned int rand = range.begin(); // some value
          for(size_t idx = range.begin(); idx < range.end(); ++idx) {
            size_t rpidx = idx % raypack.size();
            //raypack[rpidx].dir_x = raypack[rpidx].dir_x < 1024.f ? raypack[rpidx].dir_x * 1.5f : magic; // increase z axis of ray direction
            // raypack[rpidx].tfar = std::numeric_limits<float>::max();
            // If we copy the ray, then this is not overwritten.

            // RTCRayHit rayhit {raypack[rpidx], hit};

            if (rpidx == 0) {
              rand = this->fastrand(rand);
            }
            raypack[rpidx].dir_x = static_cast<float>(rand);
            rayhit.ray = raypack[rpidx];

            rayhit.ray.tfar = std::numeric_limits<float>::max();

            rtcIntersect1(scene, &context, &rayhit);

            // if (rayhit.ray.tfar > 100 /* magic number*/) {
            //   // No hit. That's a program error.
            //   std::cout << "ERROR: ray did not hit anything" << std::endl;
            //}
          }
        });



      // for (size_t idx = 0; idx < numRays; ++idx) { // Do some rays for now.

      //   size_t rpidx = idx % raypack.size();
      //   //raypack[rpidx].dir_x = raypack[rpidx].dir_x < 1024.f ? raypack[rpidx].dir_x * 1.5f : magic; // increase z axis of ray direction
      //   // raypack[rpidx].tfar = std::numeric_limits<float>::max();
      //   // If we copy the ray, then this is not overwritten.

      //   // RTCRayHit rayhit {raypack[rpidx], hit};

      //   if (rpidx == 0) {
      //     rand = this->fastrand();
      //   }
      //   raypack[rpidx].dir_x = static_cast<float>(rand);
      //   rayhit.ray = raypack[rpidx];

      //   rayhit.ray.tfar = std::numeric_limits<float>::max();

      //   // BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "Before rtcIntersect1()";
      //   // BOOST_LOG_SEV(rti::mRLogger, blt::trace) <<  "rayhit.ray.tfar=" << rayhit.ray.tfar;
      //   rtcIntersect1(scene, &context, &rayhit);
      //   // BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "After rtcIntersect1()";
      //   // BOOST_LOG_SEV(rti::mRLogger, blt::trace) <<  "rayhit.ray.tnear=" << rayhit.ray.tnear
      //   //                                          << " rayhit.ray.tfar=" << rayhit.ray.tfar
      //   //   // << " rayhit.hit.Ng_x=" << rayhit.hit.Ng_x
      //   //   // << " rayhit.hit.Ng_y=" << rayhit.hit.Ng_y
      //   //   // << " rayhit.hit.Ng_z=" << rayhit.hit.Ng_z
      //   //                                          << " rayhit.hit.primID=" << rayhit.hit.primID;
      //   // unsigned int primID = rayhit.hit.primID;
      //   // BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "hit-primID=" << primID
      //   //                                          << " " << mGeo.prim_to_string(primID);
      // }

      test_result dummyResult {};
      return dummyResult;
    }
  private:
    i_geometry_from_gmsh& mGeo;
    i_ray_source& mRaySource;

    //unsigned int g_seed = 1234;
    inline int fastrand(unsigned int pSeed) {
      pSeed = (214013*pSeed+2531011);
      return (pSeed>>16)&0x7FFF; // bs
    }
  };
}
