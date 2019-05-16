#pragma once

#include<boost/timer/timer.hpp>

#include <cmath>

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

      auto geomID = rtcAttachGeometry(scene, geometry);
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

      size_t nrexp = 20;
      size_t numRays = std::pow(2,nrexp);

      BOOST_LOG_SEV(rti::mRLogger, blt::warning)
        << "\n>>>>>>>> TODO: Compare oo-source-design with template-source-design and "
        << "no-abstraction-source-design\n";

      // Start timing until end of variable scope.
      boost::timer::auto_cpu_timer timer;

      std::cout << "Running 2**" << nrexp << " ray traces on "
                << boost::core::demangle(typeid(mGeo).name()) << std::endl << std::flush;
      for (size_t idx = 0; idx < numRays; ++idx) { // Do some rays for now.
        //RTCRay ray = rti::utils::constructRay(0,0,1, 1,1,-1); // origin (0,0,1) direction (1,1,-1)
        //RTCRay ray = rti::utils::constructRay(0,0,1, 1,0,-1); // origin (0,0,1) direction (1,0,-1)
        //RTCRay ray = rti::utils::constructRay(0,0,0.5, 1,0,0);
        //RTCRay ray = rti::utils::constructRay(0,0,0.5, 1,0.1,0); // origin (0,0,1) direction (1,1,-1)
        RTCRay ray = mRaySource.get_ray();

        RTCRayHit rayhit {ray, RTCHit {}};
        //RTCRayHit rayhit = rti::utils::constructRayHit(ray, rti::utils::constructHit()); // TODO
        BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "Before rtcIntersect1()";
        BOOST_LOG_SEV(rti::mRLogger, blt::trace) <<  "rayhit.ray.tfar=" << rayhit.ray.tfar;
        rtcIntersect1(scene, &context, &rayhit);
        BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "After rtcIntersect1()";
        BOOST_LOG_SEV(rti::mRLogger, blt::trace) <<  "rayhit.ray.tnear=" << rayhit.ray.tnear
                                                 << " rayhit.ray.tfar=" << rayhit.ray.tfar
                                                 // << " rayhit.hit.Ng_x=" << rayhit.hit.Ng_x
                                                 // << " rayhit.hit.Ng_y=" << rayhit.hit.Ng_y
                                                 // << " rayhit.hit.Ng_z=" << rayhit.hit.Ng_z
                                                 << " rayhit.hit.primID=" << rayhit.hit.primID;
        unsigned int primID = rayhit.hit.primID;
        BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "hit-primID=" << primID
                                                 << " " << mGeo.prim_to_string(primID);
      }

      test_result dummyResult {};
      return dummyResult;
    }
  private:
    i_geometry_from_gmsh& mGeo;
    i_ray_source& mRaySource;
  };
}
