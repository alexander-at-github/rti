#pragma once

#include <boost/core/demangle.hpp>

#include <cmath>
#include <chrono>
#include <omp.h>

#include <embree3/rtcore.h>

#include "rti/geo/i_boundary.hpp"
#include "rti/geo/i_geometry.hpp"
#include "rti/ray/i_source.hpp"
#include "rti/reflection/diffuse.hpp"
//#include "rti/reflection/specular.hpp"
//#include "rti/trace/counter.hpp"
#include "rti/trace/context.hpp"
#include "rti/trace/hit_accumulator.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/result.hpp"
#include "rti/util/logger.hpp"
#include "rti/util/ray_logger.hpp"
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
      assert(rtcGetDeviceProperty(pGeo.get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) != 0 &&
             "Error: Embree filter functions are not supported by your Embree instance.");
      std::cerr << "RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED == " << rtcGetDeviceProperty(pGeo.get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) << std::endl;
    }

    rti::trace::result<Ty> run() {
      // Prepare a data structure for the result.
      auto result = rti::trace::result<Ty> {};
      result.inputFilePath = mGeo.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(mGeo).name());

      // Prepare Embree
      auto device = mGeo.get_rtc_device();
      auto scene = rtcNewScene(device);

      //std::cerr << "rtc get scene flags == " << rtcGetSceneFlags(scene) << std::endl;
      // scene flags
      //rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE | RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      //std::cerr << "rtc get scene flags == " << rtcGetSceneFlags(scene) << std::endl;

      // Selecting higher build quality results in better rendering performance but slower
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      auto bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);
      auto geometry = mGeo.get_rtc_geometry();
      // rtcCommitGeometry(geometry); // Removed; should be done in the implementations of i_geometry
      auto boundary = mBoundary.get_rtc_geometry();

      // rtcAttachGeometry() is thread safe
      auto boundaryID = rtcAttachGeometry(scene, boundary);
      auto geometryID = rtcAttachGeometry(scene, geometry);

      rti::trace::context<Ty>::register_intersect_filter_funs(mGeo, mBoundary);
      assert(rtcGetDeviceError(device) == RTC_ERROR_NONE);

      // Use openMP for parallelization
      #pragma omp parallel
      {
        rtcJoinCommitScene(scene);
        // TODO: move to the other parallel region at the bottom
      }

      // *Ray queries*
      //size_t nrexp = 27;
      auto nrexp = 22; // int
      auto numRays = std::pow(2.0, nrexp); // returns a double
      result.numRays = numRays; // Save the number of rays also to the test result

      //rti::reflection::specular reflectionModel;
      //rti::reflection::diffuse reflectionModel(0.015625);
      auto reflectionModel = rti::reflection::diffuse<Ty> {};
      auto boundaryReflection = rti::reflection::specular<Ty> {};

      auto geohitc = 0ull; // unsigned long long int
      auto nongeohitc = 0ull;
      //auto hitCounter = rti::trace::counter {mGeo.get_num_primitives()};
      auto hitAccumulator = rti::trace::hit_accumulator<Ty> {mGeo.get_num_primitives()};

      // TODO: in order to make this reduction independent of the concrete type one could
      // use an abstract factory.
      #pragma omp declare \
        reduction(hit_accumulator_combine : \
                  rti::trace::hit_accumulator<Ty> : \
                  omp_out = rti::trace::hit_accumulator<Ty>(omp_out, omp_in)) \
        initializer(omp_priv = rti::trace::hit_accumulator<Ty>(omp_orig))
        //initializer(omp_priv = rti::trace::counter(omp_orig))

      // omp_set_dynamic(false);

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      auto rng = std::make_unique<rti::rng::cstdlib_rng>();

      // Start timing
      auto timer = rti::util::timer {};

      #pragma omp parallel \
        reduction(+ : geohitc, nongeohitc) \
        reduction(hit_accumulator_combine : hitAccumulator)
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

        // A dummy counter for the boundary
        auto boundaryCntr = rti::trace::dummy_counter {};

        // We will attach our data to the memory immediately following the context as described
        // in https://www.embree.org/api.html#rtcinitintersectcontext .
        auto rtiContext = rti::trace::context<Ty> {geometryID, mGeo, reflectionModel, hitAccumulator,
                                                   boundaryID, mBoundary, boundaryReflection,
                                                   *rng, *rngSeed2};
        // Initialize (also takes care for the initialization of the Embree context)
        rtiContext.init();

        std::cerr << "Thread number " << omp_get_thread_num() << std::endl;


        // std::cerr << "rtc get scene flags == " << rtcGetSceneFlags(scene) << std::endl;
        // std::cerr << "RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION == "
        // << RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION << std::endl;

        // assert((rtcGetSceneFlags(scene) & RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION) != 0 &&
        //        "Assumption: context filter-function is enabled");
        // rtcContext.filter = &rti::trace::context<Ty>::filter_fun;

        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) numRays; ++idx) {

          rayhit.ray.tnear = 0;
          rayhit.ray.time = 0;
          // prepare our custom ray tracing context
          rtiContext.rayWeight = rtiContext.INITIAL_RAY_WEIGHT;

          RLOG_DEBUG << "NEW: Preparing new ray from source" << std::endl;
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed1);

          RAYSRCLOG(rayhit);

          auto reflect = false; // initialize to some value
          do {
            RLOG_DEBUG
              << "preparing ray == ("
              << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z
              << ") ("
              << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z
              << ")" << std::endl;

            // rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            // rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            // rtiContext.geoFirstHit = true;
            // rtiContext.boundFirstHit = true;

            // // performing ray queries in a scene is thread-safe
            // rtcIntersect1(scene, &rtcContext, &rayhit);
            // rtiContext.post_process_intersect(rayhit);

            rayhit.ray.tfar = std::numeric_limits<Ty>::max();

            // Runn the intersection
            rtiContext.intersect1(scene, rayhit);

            //std::cout << "AFTER INTERSECT" << std::endl;

            RAYLOG(rayhit, (rtiContext.tfar < 50 ? rtiContext.tfar : 50));
            //RAYLOG(rayhit, rtiContext.tfar);

            // // TODO: is that correct? What if we hit surface elements first and only then the invalid area?
            // // Can that happen?
            // if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            //   // No  hit
            //   nongeohitc += 1;
            //   break; // break do-while loop
            // }

            // else
            // A hit
            reflect = rtiContext.reflect;
            RLOG_DEBUG << "tracer::run().reflect == " << (reflect ? "true" : "false")  << std::endl;

            //std::cout << "reflect == " << reflect << std::endl;

            // the following data is actually only used if a reflection happens
            // if a new ray is started from the source, these data values will be overwritten
            rayhit.ray.org_x = rtiContext.rayout[0][0];
            rayhit.ray.org_y = rtiContext.rayout[0][1];
            rayhit.ray.org_z = rtiContext.rayout[0][2];
            rayhit.ray.dir_x = rtiContext.rayout[1][0];
            rayhit.ray.dir_y = rtiContext.rayout[1][1];
            rayhit.ray.dir_z = rtiContext.rayout[1][2];
            // if (rayhit.hit.geomID == boundaryID) {
            //   // Ray hit the boundary
            //   reflect = boundaryReflection.use(rayhit, *rng, *rngSeed2, this->mBoundary);
            // } else {
            //   geohitc += 1;
            //   RLOG_DEBUG << "rayhit.hit.primID == " << rayhit.hit.primID << std::endl;
            //   RLOG_DEBUG << "prim == " << mGeo.prim_to_string(rayhit.hit.primID) << std::endl;
            //   reflect = reflectionModel.use(rayhit, *rng, *rngSeed2, this->mGeo);
            // }
          } while (reflect);
        }
      }

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitAccumulator = std::make_unique<rti::trace::hit_accumulator<Ty> >(hitAccumulator);
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      // // Turn dynamic adjustment of number of threads on again
      // omp_set_dynamic(true);

      // The geometry object will be destructed when the scene object is
      // destructed, because the geometry object is attached to the scene
      // object.
      rtcReleaseGeometry(geometry);
      rtcReleaseGeometry(boundary);


      return result;
    }
  private:
    rti::geo::i_geometry<Ty>& mGeo;
    rti::geo::i_boundary<Ty>& mBoundary;
    rti::ray::i_source& mSource;
  };
}} // namespace
