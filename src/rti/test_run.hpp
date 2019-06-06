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
#include "rti/triangle_geometry_from_gmsh.hpp" // only debug

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
      // Selecting higher build quality results in better rendering performance but slower3
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      RTCBuildQuality bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);

      RTCGeometry geometry = mGeo.get_rtc_geometry();
      rtcCommitGeometry(geometry);

      (void) rtcAttachGeometry(scene, geometry);
      rtcCommitScene(scene);

      //////////////////////////////
      // Ray queries
      //////////////////////////////


      // size_t nrexp = 27;
      // size_t nrexp = 24;
      size_t nrexp = 20;
      size_t numRays = std::pow(2,nrexp);
      result.numRays = numRays; // Save the number of rays also to the test result

      // std::cout << "Running 2**" << nrexp << " ray traces on "
      //           << boost::core::demangle(typeid(mGeo).name()) << std::endl << std::flush;

      // const std::vector<rti::pair<float>> ddxys
      // {{-3.f, -3.f}, {-1.f, -3.f}, { 1.f, -3.f}, {3.f, -3.f},
      //  {-3.f, -1.f}, {-1.f, -1.f}, { 1.f, -1.f}, {3.f, -1.f},
      //  {-3.f,  1.f}, {-1.f,  1.f}, { 1.f,  1.f}, {3.f,  1.f},
      //  {-3.f,  3.f}, {-1.f,  3.f}, { 1.f,  3.f}, {3.f,  3.f}};
      // std::vector<rti::pair<float>> ddxys
      //   {{-3.f, -3.f}, {3.f, -3.f},
      //    {-3.f,  3.f}, {3.f,  3.f}};

      // 16 coordinates evenly distributed along a square of radius 10
      std::vector<rti::pair<float> > ddxys
      {{10.f, 0.f},
         {9.238795325112868f, 3.826834323650898f},
         {7.0710678118654755f, 7.0710678118654755f},
         {3.826834323650898f, 9.238795325112868f},
       {0.f, 10.f},
         {-3.826834323650898f, 9.238795325112868f},
         {-7.0710678118654755f, 7.0710678118654755f},
         {-9.238795325112868f, 3.826834323650898f},
       {-10.f, 0.f},
         {-9.238795325112868f, -3.826834323650898f},
         {-7.0710678118654755f, -7.0710678118654755f},
         {-3.826834323650898f, -9.238795325112868f},
       {0.f, -10.f},
         {3.826834323650898f, -9.238795325112868f},
         {7.0710678118654755f, -7.0710678118654755f},
         {9.238795325112868f, -3.826834323650898f}};

      // Thread local variables
      tbb::enumerable_thread_specific<unsigned int> randGrp(1); // some number
      tbb::enumerable_thread_specific<rti::cosine_direction_source> sourceGrp;
      tbb::enumerable_thread_specific<RTCRayHit> rayhitGrp(RTCRayHit {RTCRay {}, RTCHit {}});
      tbb::enumerable_thread_specific<size_t> rpidxGrp(0);
      tbb::enumerable_thread_specific<size_t> hitcGrp(0);
      tbb::enumerable_thread_specific<size_t> nonhitcGrp(0);

      tbb::enumerable_thread_specific<std::vector<rti::triple<float> > > hitpointgrp;

      // Initialize rays
      // Initializing rays here does not work, cause TBB creates the members of an
      // enumerable_thread_specific object only when a thread requests it. That is,
      // at this position all enumerable_thread_specific containers are empty.

      // Start timing
      boost::timer::cpu_timer timer; // also calls start

      tbb::parallel_for(
        tbb::blocked_range<size_t>(0, numRays, 64),
        //[raypack = raypack] // C++14; use generalized lambda capture to copy raypack
        //[&] // Capture all by refernce can easily cause errors
        // capture by refernce
        [&scene, &ddxys, &randGrp, &sourceGrp, &rayhitGrp, &rpidxGrp, &hitcGrp, &nonhitcGrp, &hitpointgrp]
        (const tbb::blocked_range<size_t>& range) {
          // Here it is important to use the refernce type explicitly, otherwise we cannot
          // modify the content.
          //auto& rand = randGrp.local();
          auto& source = sourceGrp.local();
          auto& rayhit = rayhitGrp.local();
          auto& rpidx = rpidxGrp.local();
          auto& nonhitc = nonhitcGrp.local();
          auto& hitc = hitcGrp.local();

          // 0.1f == 0.100000001490116119384765625
          rayhit.ray.org_x = 0.1f;
          rayhit.ray.org_y = 0;
          rayhit.ray.org_z = 0;
          rayhit.ray.tnear = 0;
          // rayhit.ray.dir_x = 12.0f; // magic number
          // rayhit.ray.dir_y = 0;
          // rayhit.ray.dir_z = 0;

          rayhit.ray.time = 0;
          rayhit.ray.tfar = std::numeric_limits<float>::max();
          rayhit.ray.mask = 0;
          rayhit.ray.id = 0;
          rayhit.ray.flags = 0;

          // The range often contains only one single element. That is, the tbb::blocked_range() function
          // really is the loop.
          for(size_t idx = range.begin(); idx < range.end(); ++idx, ++rpidx /* we also maintain the rpidx here */) {

            // // Generate random direction:w

            // if (rpidx >= ddxys.size()) {
            //   rpidx = 0;
            //   rand = rand_r(&rand); // stdlib.h
            // }

            // rayhit.ray.dir_x = static_cast<float>(rand) * 4e-7f; // Number chosen by experimentation
            // // Read access to ddxys; no write.
            // rayhit.ray.dir_y = ddxys[rpidx].frst;
            // rayhit.ray.dir_z = ddxys[rpidx].scnd;

            source.set_ray(rayhit.ray);


            rayhit.ray.tfar = std::numeric_limits<float>::max();
            rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

            // The Embree tutorial uses one context per ray.
            // I think that is not necessary, though.
            RTCIntersectContext context;
            rtcInitIntersectContext(&context);

            // performing ray queries in a scene is thread-safe
            rtcIntersect1(scene, &context, &rayhit);

            // if (rayhit.ray.tfar * raypack[rpidx].dir_x > 1) {
            // if (rayhit.hit.primID == 0) {
            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
              nonhitc += 1;
              // std::cerr
              //   << std::setprecision(128)
              //   << ">>>>> ERROR: ray did not hit anything!"
              //   << " rayhit.ray.tfar == " << rayhit.ray.tfar
              //   << " origin: (" << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z << ")"
              //   << " direction: (" << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z << ")"
              //   << " primID == " << rayhit.hit.primID
              //   // << " idx == " << idx
              //   << std::endl;

              // std::cerr << "#";

              // Calculate where the hit should be and display it in Gmsh.

              // double xxhitp = 1; // Assuming hit at that x-value
              // double tt = (xxhitp - (double) rayhit.ray.org_x) / (double) rayhit.ray.dir_x;
              // double yyhitp = (double) rayhit.ray.org_y + (double) rayhit.ray.dir_y * tt;
              // double zzhitp = (double) rayhit.ray.org_z + (double) rayhit.ray.dir_z * tt;
              // gmsh::model::geo::addPoint(xxhitp, yyhitp, zzhitp);

            } else {
              hitc += 1;
              // /************************************************************************************
              //  * The following works (almost obviously) only with the triangl_geometry_from_gmsh !!
              //  ************************************************************************************/
              // rti::triple<rti::triple<float> > coords =
              //   dynamic_cast<rti::triangle_geometry_from_gmsh*>(&mGeo)->prim_to_coords(rayhit.hit.primID);
              // // gmsh::model::geo::addPoint(coords.frst.frst, coords.frst.scnd, coords.frst.thrd);
              // // gmsh::model::geo::addPoint(coords.scnd.frst, coords.scnd.scnd, coords.scnd.thrd);
              // // gmsh::model::geo::addPoint(coords.thrd.frst, coords.thrd.scnd, coords.thrd.thrd);

              // // Barycentric coordinates
              // // Note the ordering here is important! rayhit.ray.u corresponds to the y-ccordinate and
              // // rayhit.ray.v corresponds to the z-coordinate.
              // float ss = rayhit.hit.u;
              // float tt = rayhit.hit.v;
              // float rr = (1.0f - ss) - tt;
              // float xxbhit = rr * coords.frst.frst + ss * coords.scnd.frst + tt * coords.thrd.frst;
              // float yybhit = rr * coords.frst.scnd + ss * coords.scnd.scnd + tt * coords.thrd.scnd;
              // float zzbhit = rr * coords.frst.thrd + ss * coords.scnd.thrd + tt * coords.thrd.thrd;
              // float xxahit = rayhit.ray.org_x + rayhit.ray.dir_x * rayhit.ray.tfar;
              // float yyahit = rayhit.ray.org_y + rayhit.ray.dir_y * rayhit.ray.tfar;
              // float zzahit = rayhit.ray.org_z + rayhit.ray.dir_z * rayhit.ray.tfar;

              // double epsilon = 1e-3;
              // assert(std::abs(xxahit - xxbhit) < epsilon && "Error in ray hit computation");
              // assert(std::abs(yyahit - yybhit) < epsilon && "Error in ray hit computation");
              // assert(std::abs(zzahit - zzbhit) < epsilon && "Error in ray hit computation");

              // // Calling gmsh concurrently causes segmentation fault
              // //gmsh::model::geo::addPoint(xxbhit, yybhit, zzbhit);
              // //gmsh::model::geo::addPoint(xxahit, yyahit, zzahit);
              // std::vector<rti::triple<float> >& hitpoints = hitpointgrp.local();
              // hitpoints.push_back({xxahit, yyahit, zzahit});
              // //hitpointgrp.local().push_back({xxahit, yyahit, zzahit});



              // std::cerr
              //   //   << "<<<<<< no error in this iteration"
              //   //   << " rayhit.ray.tfar == " << rayhit.ray.tfar
              //   << " origin: (" << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z << ")"
              //   << " direction: (" << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z << ")"
              //   << " hit: (" << xxbhit << " " << yybhit << " " << zzbhit << ")"
              //   //   << " primID == " << rayhit.hit.primID
              //   //   << " " << mGeo.prim_to_string(rayhit.hit.primID)
              //   << std::endl;
              //std::cerr << "*";

              // double xxhitp = (double)rayhit.ray.org_x + (double)rayhit.ray.dir_x * (double)rayhit.ray.tfar;
              // double yyhitp = (double)rayhit.ray.org_y + (double)rayhit.ray.dir_y * (double)rayhit.ray.tfar;
              // double zzhitp = (double)rayhit.ray.org_z + (double)rayhit.ray.dir_z * (double)rayhit.ray.tfar;

              // std::cerr << ">>>>>>> writing hit point (" << xxhitp << " " << yyhitp << " " << zzhitp << ")"
              //           << std::endl;

              // gmsh::model::geo::addPoint(xxhitp, yyhitp, zzhitp);
            }
          }
        });

        result.timeNanoseconds = timer.elapsed().wall; // Using boost timer


        // Release geomtery.
        // The geometry object will be destructed when the scene object is
        // destructed, because the geometry object is attached to the scene
        // object.
        rtcReleaseGeometry(geometry);


        // for (auto& grp : hitpointgrp) {
        //   for (auto& xyz : grp) {
        //     gmsh::model::geo::addPoint(xyz.frst, xyz.scnd, xyz.thrd);
        //   }
        // }
        // gmsh::model::geo::synchronize();

        result.hitc = hitcGrp.combine([](size_t xx, size_t yy){ return xx + yy; });
        result.nonhitc = nonhitcGrp.combine([](size_t xx, size_t yy){ return xx + yy; });

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
