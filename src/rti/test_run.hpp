#pragma once

#include <boost/core/demangle.hpp>
#include <boost/math/constants/constants.hpp> // pi
#include <boost/timer/timer.hpp>

#include <cmath>
#include <chrono>
#include <iomanip> // std::setprecision()

#include "tbb/enumerable_thread_specific.h"
#include "tbb/tbb.h"

#include "rti/i_geometry.hpp"
#include "rti/i_ray_source.hpp"
#include "rti/test_result.hpp"

// Debug
#include <gmsh.h>

#include "rti/triangle_geometry_from_gmsh.hpp"

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
      //rtcReleaseGeometry(geometry);

      //////////////////////////////
      // Ray queries
      //////////////////////////////
      // The Embree tutorial uses one context per ray.
      // I think that is not necessary, though.
      RTCIntersectContext context;
      rtcInitIntersectContext(&context);

      // size_t nrexp = 26;
      // size_t nrexp = 24;
      size_t nrexp = 20;
      size_t numRays = std::pow(2,nrexp);
      result.numRays = numRays; // Save the number of rays also to the test result

      // std::cout << "Running 2**" << nrexp << " ray traces on "
      //           << boost::core::demangle(typeid(mGeo).name()) << std::endl << std::flush;

      // std::vector<rti::pair<float>> ddxys
      // {{-3.f, -3.f}, {-1.f, -3.f}, { 1.f, -3.f}, {3.f, -3.f},
      //  {-3.f, -1.f}, {-1.f, -1.f}, { 1.f, -1.f}, {3.f, -1.f},
      //  {-3.f,  1.f}, {-1.f,  1.f}, { 1.f,  1.f}, {3.f,  1.f},
      //  {-3.f,  3.f}, {-1.f,  3.f}, { 1.f,  3.f}, {3.f,  3.f}};
      std::vector<rti::pair<float>> ddxys
      {{-3.f, -3.f}, {3.f, -3.f},
       {-3.f,  3.f}, {3.f,  3.f}};

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
      // tbb::enumerable_thread_specific<unsigned int> randGrp(1234); // some number
      // tbb::enumerable_thread_specific<RTCRayHit> rayhitGrp; //(RTCRayHit {RTCRay {}, RTCHit {}});
      // tbb::enumerable_thread_specific<size_t> rpidxGrp(0);

      // std::cerr
      //   << "rayhitGrp.size() == " << rayhitGrp.size()
      //   << std::endl;
      // all enumerable_thread_specific<> objects are empty before a thread requests them.

      // // Initialize rays
      // for (auto& rayhit : rayhitGrp) {
      //   rayhit.ray.org_x = 1.0f;
      //   rayhit.ray.org_y = 0;
      //   rayhit.ray.org_z = 0;
      //   rayhit.ray.tnear = 0;
      //   rayhit.ray.dir_x = 0.5f; // magic number
      //   rayhit.ray.dir_y = 0;
      //   rayhit.ray.dir_z = 0;
      //   rayhit.ray.time = 0;
      //   rayhit.ray.tfar = std::numeric_limits<float>::max();
      //   rayhit.ray.mask = 0;
      //   rayhit.ray.id = 0;
      //   rayhit.ray.flags = 0;
      // }

      // Start timing until end of variable scope.
      //boost::timer::auto_cpu_timer timerOld;

      boost::timer::cpu_timer timer; // also calls start

      result.startTime = std::chrono::high_resolution_clock::now();

      unsigned int rand = 1;
      RTCRayHit rayhit;
      rayhit.ray.dir_x = 0;

      uint64_t hitcc = 0;
      uint64_t errcc = 0;

      size_t idx1  = 0;
      size_t rpidx = 0;
      size_t rpidx2 = 0;
      for (idx1 = 0, rpidx = 0; idx1 < numRays; ++idx1, ++rpidx) {
        if (idx1 % std::max<size_t>(1, numRays / 100) == 0) {std::cerr << "*";}

      // tbb::parallel_for(
      //   tbb::blocked_range<size_t>(0, numRays),
        //[raypack = raypack] // C++14; use generalized lambda capture to copy raypack
        //[&] // Capture all by refernce can easily cause errors
        // [&scene, &context, &ddxys, &randGrp, &rayhitGrp, &rpidxGrp] // capture by refernce
        // (const tbb::blocked_range<size_t>& range) {
          // tbb::enumerable_thread_specific<unsigned int>::reference rand = randGrp.local();
          // tbb::enumerable_thread_specific<RTCRayHit>::reference rayhit = rayhitGrp.local();
          // tbb::enumerable_thread_specific<size_t>::reference rpidx = rpidxGrp.local();
          //
          // Here it is important to use the refernce type explicitly, otherwise we cannot
          // modify the content.
          // auto& rand = randGrp.local();
          // auto& rayhit = rayhitGrp.local();
          // auto& rpidx = rpidxGrp.local();


          // 0.1f == 0.100000001490116119384765625
          rayhit.ray.org_x = 0.1f;
          rayhit.ray.org_y = 0;
          rayhit.ray.org_z = 0;
          rayhit.ray.tnear = 0;
          rayhit.ray.dir_x = 12.0f; // magic number
          rayhit.ray.dir_y = 0;
          rayhit.ray.dir_z = 0;
          rayhit.ray.time = 0;
          rayhit.ray.tfar = std::numeric_limits<float>::max();
          rayhit.ray.mask = 0;
          rayhit.ray.id = 0;
          rayhit.ray.flags = 0;


          // The range often contains only one single element. That is, the tbb::blocked_range() function
          // really is the loop.
          //for(size_t idx = range.begin(); idx < range.end(); ++idx) {

            // We use an index over the coordinates saved in ddxys, which
            // we have to maintain our selfs.
            if (rpidx >= ddxys.size()) {
              //std::cout << "[rand == " << rand;
              rpidx = 0;
              rpidx2 += 1;

              double rotfrac = 32;
              // if (rpidx2 > rotfrac - 1) {
              //   rpidx2 = 0;
              //   rand += 1;
              // }

              //rand = fastrand(rand);
              for (auto & xxyy : ddxys) {
                double cos = std::cos(boost::math::constants::pi<double>() / (rotfrac / 2));
                double sin = std::sin(boost::math::constants::pi<double>() / (rotfrac / 2));
                double frst = cos * xxyy.frst - sin * xxyy.scnd;
                double scnd = sin * xxyy.frst + cos * xxyy.scnd;
                xxyy.frst = frst;
                xxyy.scnd = scnd;
              }
              //rand = rand_r(&rand); // stdlib.h
              //std::cout << " new rand == " << rand << "]";
            }
            //BOOST_LOG_SEV(rti::mRLogger, blt::warning) << "rpidx == " << rpidx;
            rand += 1;


            // A ray from (0.1 0 0) in direction (1493484032 1 1) produces an error.

            // rayhit.ray.dir_x = std::max<float>(static_cast<float>(rand), 21);
            rayhit.ray.dir_x = static_cast<float>(rand) / 10000;
            rayhit.ray.dir_y = ddxys[rpidx].frst;
            rayhit.ray.dir_z = ddxys[rpidx].scnd;

            // rayhit.ray.dir_x = 704282496; // static_cast<float>(rand);
            // rayhit.ray.dir_y = 3;  // ddxys[rpidx].frst;
            // rayhit.ray.dir_z = -1; // ddxys[rpidx].scnd;

            // rayhit.ray.dir_x = 1;
            // rayhit.ray.dir_y = 0;
            // rayhit.ray.dir_z = 0;

            rayhit.ray.tfar = std::numeric_limits<float>::max();
            rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

            // debug
            rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;
            rayhit.ray.tnear = 0;

            // Do we need that?
            RTCIntersectContext context;
            rtcInitIntersectContext(&context);

            //std::cerr << " [idx == " << idx << "] ";

            rtcIntersect1(scene, &context, &rayhit);

            // if (rayhit.ray.tfar * raypack[rpidx].dir_x > 1) {
            // if (rayhit.hit.primID == 0) {
            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
              errcc += 1;
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
              hitcc += 1;
              rti::triple<rti::triple<float> > coords =
                dynamic_cast<rti::triangle_geometry_from_gmsh*>(&mGeo)->prim_to_coords(rayhit.hit.primID);
              // gmsh::model::geo::addPoint(coords.frst.frst, coords.frst.scnd, coords.frst.thrd);
              // gmsh::model::geo::addPoint(coords.scnd.frst, coords.scnd.scnd, coords.scnd.thrd);
              // gmsh::model::geo::addPoint(coords.thrd.frst, coords.thrd.scnd, coords.thrd.thrd);

              // Barycentric coordinates
              // Note the ordering here is important! rayhit.ray.u corresponds to the y-ccordinate and
              // rayhit.ray.v corresponds to the z-coordinate.
              float ss = rayhit.hit.u;
              float tt = rayhit.hit.v;
              float rr = (1.0f - ss) - tt;
              float xxbhit = rr * coords.frst.frst + ss * coords.scnd.frst + tt * coords.thrd.frst;
              float yybhit = rr * coords.frst.scnd + ss * coords.scnd.scnd + tt * coords.thrd.scnd;
              float zzbhit = rr * coords.frst.thrd + ss * coords.scnd.thrd + tt * coords.thrd.thrd;
              float xxahit = rayhit.ray.org_x + rayhit.ray.dir_x * rayhit.ray.tfar;
              float yyahit = rayhit.ray.org_y + rayhit.ray.dir_y * rayhit.ray.tfar;
              float zzahit = rayhit.ray.org_z + rayhit.ray.dir_z * rayhit.ray.tfar;


              double epsilon = 1e-6;
              assert(xxahit - xxbhit < epsilon && "Error in ray hit computation");
              assert(yyahit - yybhit < epsilon && "Error in ray hit computation");
              assert(zzahit - zzbhit < epsilon && "Error in ray hit computation");

              gmsh::model::geo::addPoint(xxbhit, yybhit, zzbhit);
              //gmsh::model::geo::addPoint(xxahit, yyahit, zzahit);

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
            //});


        // std::cerr
        //   << "rayhitGrp.size() == " << rayhitGrp.size()
        //   << std::endl;

        result.endTime = std::chrono::high_resolution_clock::now();

        result.timeNanoseconds = timer.elapsed().wall; // Using boost timer

        rtcReleaseGeometry(geometry);

        std::cerr << "hitcc == " << hitcc << " errcc == " << errcc << std::endl;

        gmsh::model::geo::synchronize();
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
