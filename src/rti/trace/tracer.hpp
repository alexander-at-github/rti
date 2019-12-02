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
#include "rti/rng/cstdlib_rng.hpp"
//#include "rti/trace/counter.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/hit_accumulator.hpp"
#include "rti/trace/hit_accumulator_with_checks.hpp"
#include "rti/trace/point_cloud_context.hpp"
#include "rti/trace/result.hpp"
#include "rti/util/logger.hpp"
#include "rti/util/ray_logger.hpp"
#include "rti/util/timer.hpp"

namespace rti { namespace trace {
  template<typename Ty>
  class tracer {
  public:
    tracer(rti::geo::i_factory<Ty>& pFactory,
           rti::geo::i_boundary<Ty>& pBoundary,
           rti::ray::i_source& pSource,
           size_t pNumRays) :
      mFactory(pFactory),
      mBoundary(pBoundary),
      mSource(pSource),
      mNumRays(pNumRays) {
      assert (mFactory.get_geometry().get_rtc_device() == pBoundary.get_rtc_device() &&
              "the geometry and the boundary need to refer to the same Embree (rtc) device");
      assert(rtcGetDeviceProperty(mFactory.get_geometry().get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) != 0 &&
             "Error: Embree filter functions are not supported by your Embree instance.");
      // std::cerr << "RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED == " << rtcGetDeviceProperty(mFactory.get_geometry().get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) << std::endl;
      auto& device = mFactory.get_geometry().get_rtc_device();
      assert(rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_VERSION) >= 30601 && "Error: The minimum version of Embree is 3.6.1");
      std::cerr << "Warning: tnear set to a constant! FIX" << std::endl;
    }

    rti::trace::result<Ty> run() {
      // Prepare a data structure for the result.
      auto result = rti::trace::result<Ty> {};
      auto& geo = mFactory.get_geometry();
      result.inputFilePath = geo.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(geo).name());

      // Prepare Embree
      auto device = geo.get_rtc_device();
      auto scene = rtcNewScene(device);

      // scene flags
      // Does one need to set this flag if one uses the registered call-back functions only?
      //rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE | RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      //std::cerr << "rtc get scene flags == " << rtcGetSceneFlags(scene) << std::endl;

      // Selecting higher build quality results in better rendering performance but slower
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      auto bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);
      auto geometry = geo.get_rtc_geometry();
      // rtcCommitGeometry(geometry); // Removed; should be done in the implementations of i_geometry
      auto boundary = mBoundary.get_rtc_geometry();

      // rtcAttachGeometry() is thread safe
      auto boundaryID = rtcAttachGeometry(scene, boundary);
      auto geometryID = rtcAttachGeometry(scene, geometry);

      mFactory.register_intersect_filter_funs(mBoundary);
      assert(rtcGetDeviceError(device) == RTC_ERROR_NONE && "Error");

      // Use openMP for parallelization
      #pragma omp parallel
      {
        rtcJoinCommitScene(scene);
        // TODO: move to the other parallel region at the bottom
      }

      // Ray queries
      result.numRays = mNumRays; // Save the number of rays also to the test result

      auto reflectionModel = rti::reflection::diffuse<Ty> {};
      auto boundaryReflection = rti::reflection::specular<Ty> {};

      auto geohitc = 0ull; // unsigned long long int
      auto nongeohitc = 0ull;
      //auto hitCounter = rti::trace::counter {mGeo.get_num_primitives()};
      auto hitAccumulator = rti::trace::hit_accumulator_with_checks<Ty> {geo.get_num_primitives()};

      // in order to make this reduction independent of the concrete type one could
      // use an abstract factory.
      #pragma omp declare \
        reduction(hit_accumulator_combine : \
                  rti::trace::hit_accumulator_with_checks<Ty> : \
                  omp_out = rti::trace::hit_accumulator_with_checks<Ty>(omp_out, omp_in)) \
        initializer(omp_priv = rti::trace::hit_accumulator_with_checks<Ty>(omp_orig))
        //initializer(omp_priv = rti::trace::counter(omp_orig))

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
        alignas(128) auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

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
        // Note: the memory layout only works with plain old data, that is, C-style structs.
        // Otherwise the compiler might change the memory layout, e.g., with the vtable.
        auto rtiContext = mFactory.get_new_context(geometryID, geo, reflectionModel,
                                                   hitAccumulator, boundaryID, mBoundary,
                                                   boundaryReflection, *rng, *rngSeed2);

        // Initialize (also takes care for the initialization of the Embree context)
        rtiContext->init();

        //std::cerr << "Thread number " << omp_get_thread_num() << std::endl;

        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) mNumRays; ++idx) {

          // Note: Embrees backface culling does not solve our problem of intersections
          // when starting a new ray very close to or a tiny bit below the surface.
          // For that reason we set tnear to some value.
          // There is a risk though: when setting tnear to some strictly positive value
          // we depend on the length of the direction vector of the ray (rayhit.ray.dir_X).

          // prepare our custom ray tracing context
          rtiContext->init_ray_weight();

          RLOG_DEBUG << "NEW: Preparing new ray from source" << std::endl;
          // TODO: FIX: tnear set to a constant
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed1); // fills also tnear!

          RAYSRCLOG(rayhit);

          auto reflect = false; // initialize to some value
          do {
            RLOG_DEBUG
              << "preparing ray == ("
              << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z
              << ") ("
              << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z
              << ")" << std::endl;

            rayhit.ray.tfar = std::numeric_limits<Ty>::max();

            // Runn the intersection
            rtiContext->intersect1(scene, rayhit);

            RAYLOG(rayhit, rtiContext->tfar);

            // // TODO: is that correct? What if we hit surface elements first and only then the invalid area?
            // // Can that happen?
            // if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            //   // No  hit
            //   nongeohitc += 1;
            //   break; // break do-while loop
            // }

            // else
            // A hit
            reflect = rtiContext->reflect;
            auto hitpoint = rti::util::triple<float> {rayhit.ray.org_x + rayhit.ray.dir_x * rtiContext->tfar,
                                                      rayhit.ray.org_y + rayhit.ray.dir_y * rtiContext->tfar,
                                                      rayhit.ray.org_z + rayhit.ray.dir_z * rtiContext->tfar};
            RLOG_DEBUG
              << "tracer::run(): hit-point: " << hitpoint[0] << " " << hitpoint[1] << " " << hitpoint[2]
              << " reflect == " << (reflect ? "true" : "false")  << std::endl;

            //std::cout << "reflect == " << reflect << std::endl;

            // the following data is actually only used if a reflection happens
            // if a new ray is started from the source, these data values will be overwritten

            // ATTENTION tnear is set in another function, too! When the ray starts from the source, then
            // the source class also sets tnear!
            auto tnear = 1e-4f; // float
            // Same holds for time
            auto time = 0.0f; // float
            // rayhit.ray.org_x = rtiContext->rayout[0][0];
            // rayhit.ray.org_y = rtiContext->rayout[0][1];
            // rayhit.ray.org_z = rtiContext->rayout[0][2];
            // rayhit.ray.tnear = tnear;
            // rayhit.ray.dir_x = rtiContext->rayout[1][0];
            // rayhit.ray.dir_y = rtiContext->rayout[1][1];
            // rayhit.ray.dir_z = rtiContext->rayout[1][2];
            // rayhit.time = time;


            // float vara[4] = {(float) rtiContext->rayout[0][0],
            //                  (float) rtiContext->rayout[0][1],
            //                  (float) rtiContext->rayout[0][2],
            //                  tnear};
            // float varb[4] = {(float) rtiContext->rayout[1][0],
            //                  (float) rtiContext->rayout[1][1],
            //                  (float) rtiContext->rayout[1][2],
            //                  time};
            // reinterpret_cast<__m128&>(rayhit.ray) = _mm_load_ps(vara);
            // reinterpret_cast<__m128&>(rayhit.ray.dir_x) = _mm_load_ps(varb);

            reinterpret_cast<__m128&>(rayhit.ray) = _mm_set_ps(tnear, (float)
                                                         (float) rtiContext->rayout[0][2],
                                                         (float) rtiContext->rayout[0][1],
                                                         (float) rtiContext->rayout[0][0]
                                                         );
            reinterpret_cast<__m128&>(rayhit.ray.dir_x) = _mm_set_ps(time,
                                                               (float) rtiContext->rayout[1][2],
                                                               (float) rtiContext->rayout[1][1],
                                                               (float) rtiContext->rayout[1][0]
                                                               );





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
      result.hitAccumulator = std::make_unique<rti::trace::hit_accumulator_with_checks<Ty> >(hitAccumulator);
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      // The geometry object will be destructed when the scene object is
      // destructed, because the geometry object is attached to the scene
      // object.
      rtcReleaseGeometry(geometry);
      rtcReleaseGeometry(boundary);


      // { // Count buckets
      //   auto numOfBuckets = 1024;
      //   auto zvalToValueMap = std::vector<rti::util::pair<Ty> > {};
      //   for (unsigned int primid = 0 ; primid < mGeo.get_num_primitives(); ++primid) {
      //     auto value = hitAccumulator.get_values()[primid];
      //     auto coords = mGeo.get_prim(primid);
      //     auto zcoord = coords[2];
      //     zvalToValueMap.push_back({zcoord, value});
      //   }
      //   // sort vector by first element of pair
      //   std::sort(zvalToValueMap.begin(), zvalToValueMap.end(),
      //             [](rti::util::pair<Ty> p1, rti::util::pair<Ty> p2) { return p1[0] < p2[0]; });
      //   //
      //   auto deltaZ = std::fabs(zvalToValueMap.front()[0] - zvalToValueMap.back()[0]) / numOfBuckets;
      //   auto cBucketMinZ = zvalToValueMap.front()[0];
      //   //
      //   auto lastzcoord = std::numeric_limits<Ty>::lowest();
      //   auto valuesum = (Ty) 0;
      //   for (auto const& pair : zvalToValueMap) {
      //     auto zcoord = pair[0];
      //     auto value = pair[1];
      //
      //     assert(cBucketMinZ <= zcoord && "Error");
      //     if (zcoord < cBucketMinZ + deltaZ) {
      //       valuesum += value;
      //       continue;
      //     }
      //     // else
      //     assert(cBucketMinZ + deltaZ <= zcoord && "Error");
      //     // write bucket out
      //     std::cerr << cBucketMinZ + (deltaZ/2) << " " << valuesum << std::endl;
      //     // move to next bucket
      //     cBucketMinZ += deltaZ;
      //     valuesum = value;
      //   }
      //   std::cerr << cBucketMinZ + (deltaZ/2) << " " << valuesum << std::endl;
      // }

      return result;
    }
  private:
    rti::geo::i_factory<Ty>& mFactory;
    rti::geo::i_boundary<Ty>& mBoundary;
    rti::ray::i_source& mSource;
    size_t mNumRays;
  };
}} // namespace
