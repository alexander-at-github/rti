#pragma once

// Debug
#include <ctime>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>


#include <chrono>
#include <cmath>
#include <iostream>
#include <omp.h>

#include <embree3/rtcore.h>

#include "../geo/absc_point_cloud_geometry.hpp"
#include "../geo/disc_bounding_box_intersector.hpp"
#include "../geo/i_boundary.hpp"
#include "../geo/i_geometry.hpp"
#include "../geo/point_cloud_disc_factory.hpp"
#include "../io/vtp_writer.hpp"
#include "../particle/i_particle.hpp"
#include "../particle/i_particle_factory.hpp"
#include "../ray/constant_direction.hpp"
#include "../ray/disc_origin.hpp"
#include "../ray/i_source.hpp"
#include "../reflection/diffuse.hpp"
//#include "../reflection/specular.hpp"
#include "../rng/cstdlib_rng.hpp"
#include "../rng/mt64_rng.hpp"
#include "dummy_counter.hpp"
#include "hit_accumulator.hpp"
#include "point_cloud_context.hpp"
#include "result.hpp"
#include "../util/logger.hpp"
#include "../util/ray_logger.hpp"
#include "../util/timer.hpp"

namespace rti { namespace trace {
  template<typename numeric_type>
  class tracer {

  public:
    tracer(rti::geo::point_cloud_disc_factory<numeric_type,  rti::trace::point_cloud_context<numeric_type> >& pFactory,
           rti::geo::i_boundary<numeric_type>& pBoundary,
           rti::ray::i_source& pSource,
           size_t pNumRays,
           rti::particle::i_particle_factory<numeric_type>& particlefactory) :
      mFactory(pFactory),
      mBoundary(pBoundary),
      mSource(pSource),
      mNumRays(pNumRays),
      particlefactory(particlefactory) {
      assert (mFactory.get_geometry().get_rtc_device() == pBoundary.get_rtc_device() &&
              "the geometry and the boundary need to refer to the same Embree (rtc) device");
      assert(rtcGetDeviceProperty(mFactory.get_geometry().get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) != 0 &&
             "Error: Embree filter functions are not supported by your Embree instance.");
      // std::cerr
      //   << "RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED == "
      //   << rtcGetDeviceProperty(mFactory.get_geometry().get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED)
      //   << std::endl;
      auto& device = mFactory.get_geometry().get_rtc_device();
      assert(rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_VERSION) >= 30601 && "Error: The minimum version of Embree is 3.6.1");
      if (rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED) != 0) {
        std::cerr << "=== Warning: Embree backface culling enabled. This may result in incorrect results "
                  << "(for triangles)" << std::endl;
        assert( false && "Error: backface culling is enabled; as a consequence for triangles \
                the tracer depends on the order of the vertices of the triangles");
      }
      RLOG_WARNING << "Warning: tnear set to a constant! FIX" << std::endl;
    }

  private:
    void if_RLOG_PROGRESS_is_set_print_progress(size_t& raycnt, size_t totalnumrays) {
      if (omp_get_thread_num() != 0) {
        return;
      }
      auto barlength = 30u;
      auto barstartsymbol = '[';
      auto fillsymbol = '#';
      auto emptysymbol = '-';
      auto barendsymbol = ']';
      auto percentagestringformatlength = 3; // 3 digits

      if (raycnt % (int) std::ceil((float) totalnumrays / omp_get_num_threads() / barlength) == 0) {
        auto filllength = (int) std::ceil(raycnt / ((float)totalnumrays / omp_get_num_threads() / barlength));
        auto percentagestring = std::to_string((filllength * 100) / barlength);
        percentagestring =
          std::string(percentagestringformatlength - percentagestring.length(), ' ') +
          percentagestring + "%";
        auto bar =
          "" + std::string(1, barstartsymbol) +
          std::string(filllength, fillsymbol) + std::string(std::max(0, (int) barlength - (int) filllength), emptysymbol) +
          std::string(1, barendsymbol) + " " + percentagestring ;
        RLOG_PROGRESS << "\r" << bar;
        if (filllength >= barlength) {
          RLOG_PROGRESS << std::endl;
        }
      }
      raycnt += 1;
    }

    void
    compute_exposed_areas_by_sampling
    (rti::geo::i_geometry<numeric_type>& geo,
     RTCScene& scene,
     rti::trace::i_hit_accumulator<numeric_type>& hitacc,
     unsigned int geometryID,
     rti::rng::i_rng& rng,
     rti::rng::i_rng::i_state& rngstate
     )
    {
      // Precondition: Embree has been set up properly
      //   (e.g., the geometry has been built)
      // Precondition: We are in an OpenMP parallel block
      // Precondition: the following cast characterizes a precondition:
      auto pcgeo =
        dynamic_cast<rti::geo::absc_point_cloud_geometry<numeric_type>*> (&geo);

      std::cerr << "### FIX: This function does not consider that areas of "
                << "discs may be located outside of the boundary" << std::endl;

      // We use a new RTC context object.
      auto context = RTCIntersectContext {};
      rtcInitIntersectContext(&context);
      auto numOfSamples = 1024u;

      auto numofprimitives = pcgeo->get_num_primitives();
      auto areas = std::vector<double> (numofprimitives, 0);
      #pragma omp for
      for (size_t primidx = 0; primidx < numofprimitives; ++primidx) {
        auto pointradius = pcgeo->get_prim(primidx);
        auto normal = pcgeo->get_normal(primidx);
        auto point = rti::util::triple<numeric_type>
          {pointradius[0], pointradius[1], pointradius[2]};
        auto radius = pointradius[3];
        auto origincenter = rti::util::sum(point, rti::util::scale(2*radius, normal));
        auto invnormal = rti::util::inv(normal);
        auto origin = rti::ray::disc_origin<numeric_type>
          {origincenter, invnormal, radius};
        auto direction = rti::ray::constant_direction<numeric_type> {invnormal};
        auto source = rti::ray::source<numeric_type> {origin, direction};

        alignas(128) auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        auto hits = 0u;
        for (size_t sidx = 0; sidx < numOfSamples; ++sidx) {
          rayhit.ray.tnear = 0;
          rayhit.ray.time = 0;
          rayhit.ray.tfar = std::numeric_limits<float>::max();
          rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
          rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
          source.fill_ray(rayhit.ray, rng, rngstate);

          RAYSRCLOG(rayhit);

          rtcIntersect1(scene, &context, &rayhit);

          RAYLOG(rayhit, rayhit.ray.tfar);

          if (rayhit.hit.geomID != geometryID)
            continue;
          if (rayhit.hit.primID != primidx)
            continue;
          hits += 1;
        }
        assert (0 <= hits && hits <= numOfSamples && "Correctness assertion");
        //std::cerr << " " << (double) hits / numOfSamples;
        areas[primidx] = pcgeo->get_area(primidx) * (double) hits / numOfSamples;
      }
      hitacc.set_exposed_areas(areas);
    }

    void
    use_entire_areas_of_primitives_as_exposed
    (rti::geo::i_geometry<numeric_type>& geo,
     rti::trace::i_hit_accumulator<numeric_type>& hitacc)
    {
      // Precondition: We are in an OpenMP parallel block
      auto numofprimitives = geo.get_num_primitives();
      auto areas = std::vector<double> (numofprimitives, 0);
      #pragma omp for
      for (size_t idx = 0; idx < numofprimitives; ++idx) {
        areas[idx] = geo.get_area(idx);
      }
      // Note: effectively in each thread only some of the area-values have
      // been set. Also each thread holds its own hit-accumulator.
      hitacc.set_exposed_areas(areas);
    }

  public:
    rti::trace::result<numeric_type> run()
    {
      // Prepare a data structure for the result.
      auto result = rti::trace::result<numeric_type> {};
      auto& geo = *dynamic_cast<rti::geo::absc_point_cloud_geometry<numeric_type>* > (&mFactory.get_geometry());
      result.inputFilePath = geo.get_input_file_path();
      result.geometryClassName = typeid(geo).name();

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
      auto boundary = mBoundary.get_rtc_geometry();

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

      result.numRays = mNumRays;

      auto reflectionModel = rti::reflection::diffuse<numeric_type> {};
      auto boundaryReflection = rti::reflection::specular<numeric_type> {};

      auto geohitc = 0ull;
      auto nongeohitc = 0ull;
      auto hitAccumulator = rti::trace::hit_accumulator<numeric_type> {geo.get_num_primitives()};

      #pragma omp declare \
        reduction(hit_accumulator_combine : \
                  rti::trace::hit_accumulator<numeric_type> : \
                  omp_out = rti::trace::hit_accumulator<numeric_type>(omp_out, omp_in)) \
        initializer(omp_priv = rti::trace::hit_accumulator<numeric_type>(omp_orig))

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      // auto rng = std::make_unique<rti::rng::cstdlib_rng>();
      auto rng = std::make_unique<rti::rng::mt64_rng>();

      // Start timing
      auto timer = rti::util::timer {};

      #pragma omp parallel \
        reduction(+ : geohitc, nongeohitc) \
        reduction(hit_accumulator_combine : hitAccumulator)
      {
        // Thread local data goes here, if it is not needed anymore after the execution of the parallel region.
        alignas(128) auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        // REMARK: All data which is modified in the parallel loop should be
        // handled here explicitely.
        auto seed = (unsigned int) ((omp_get_thread_num() + 1) * 29); // multiply by magic number (prime)
        // It seems really important to use two separate seeds / states for
        // sampling the source and sampling reflections. When we use only one
        // state for both, then the variance is very high.
        auto rngSeed1 = std::make_unique<rti::rng::mt64_rng::state>(seed);
        auto rngSeed2 = std::make_unique<rti::rng::mt64_rng::state>(seed+2);

        // A dummy counter for the boundary
        auto boundaryCntr = rti::trace::dummy_counter {};

        // thread-local particle
        auto particle = particlefactory.create();

        // We will attach our data to the memory immediately following the context as described
        // in https://www.embree.org/api.html#rtcinitintersectcontext .
        // Note: the memory layout only works with plain old data, that is, C-style structs.
        // Otherwise the compiler might change the memory layout, e.g., with the vtable.
        auto rtiContext = mFactory.get_new_context(geometryID, geo, reflectionModel,
                                                   hitAccumulator, boundaryID, mBoundary,
                                                   boundaryReflection, *rng, *rngSeed2,
                                                   *particle);

        // Initialize (also takes care for the initialization of the Embree context)
        rtiContext->init();

        size_t raycnt = 0;

        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) mNumRays; ++idx) {

          // Note: Embrees backface culling does not solve our problem of intersections
          // when starting a new ray very close to or a tiny bit below the surface.
          // For that reason we set tnear to some value.
          // There is a risk though: when setting tnear to some strictly positive value
          // we depend on the length of the direction vector of the ray (rayhit.ray.dir_X).

          particle->init_new();
          // prepare our custom ray tracing context
          rtiContext->init_ray_weight();

          RLOG_DEBUG << "NEW: Preparing new ray from source" << std::endl;
          // TODO: FIX: tnear set to a constant
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed1); // fills also tnear!

          RAYSRCLOG(rayhit);

          if_RLOG_PROGRESS_is_set_print_progress(raycnt, mNumRays);

          auto reflect = false;
          do {
            RLOG_DEBUG
              << "preparing ray == ("
              << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z
              << ") ("
              << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z
              << ")" << std::endl;

            rayhit.ray.tfar = std::numeric_limits<numeric_type>::max();
            
            // Runn the intersection
            rtiContext->intersect1(scene, rayhit);

            RAYLOG(rayhit, rtiContext->tfar);

            reflect = rtiContext->reflect;
            auto hitpoint = rti::util::triple<numeric_type> {rayhit.ray.org_x + rayhit.ray.dir_x * rtiContext->tfar,
                                                   rayhit.ray.org_y + rayhit.ray.dir_y * rtiContext->tfar,
                                                   rayhit.ray.org_z + rayhit.ray.dir_z * rtiContext->tfar};
            RLOG_DEBUG
              << "tracer::run(): hit-point: " << hitpoint[0] << " " << hitpoint[1] << " " << hitpoint[2]
              << " reflect == " << (reflect ? "true" : "false")  << std::endl;

            // ATTENTION tnear is set in another function, too! When the ray starts from the source, then
            // the source class also sets tnear!
            auto tnear = 1e-4f; // float
            // Same holds for time
            auto time = 0.0f; // float
            // reinterpret_cast<__m128&>(rayhit.ray) = _mm_load_ps(vara);
            // reinterpret_cast<__m128&>(rayhit.ray.dir_x) = _mm_load_ps(varb);

            // For the store-forwarding stall also see the following link:
            // https://stackoverflow.com/questions/49265634/what-is-the-difference-between-loadu-ps-and-set-ps-when-using-unformatted-data
            // OPTIMIZATION: It might be better to change the memory layout of the rayout variable
            // such that tnear (a constant) is saved right in front of rayout. Then we can use two
            // __m128d _mm_store_sd (__m128d a, __m128d b) to write 128 bits in two chunks of 64 bits
            // into the destination.
            reinterpret_cast<__m128&>(rayhit.ray) =
              _mm_set_ps(tnear,
                         (float) rtiContext->rayout[0][2],
                         (float) rtiContext->rayout[0][1],
                         (float) rtiContext->rayout[0][0]);
            reinterpret_cast<__m128&>(rayhit.ray.dir_x) =
              _mm_set_ps(time,
                         (float) rtiContext->rayout[1][2],
                         (float) rtiContext->rayout[1][1],
                         (float) rtiContext->rayout[1][0]);

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
        // if (rtiContext->compute_exposed_areas_by_sampling()) {
        //   std::cerr
        //     << "###############"
        //     << "### WARNING ### Computing exposed area by sampling does not work"
        //     << "###############" << std::endl;
        //   // // Embree Documentation: "Passing NULL as function pointer disables the registered callback function."
        //   // rtcSetGeometryIntersectFilterFunction(geometry, nullptr);
        //   // rtcSetGeometryIntersectFilterFunction(boundary, nullptr);
        //   // rtcCommitGeometry(geometry);
        //   // rtcCommitGeometry(boundary);
        //   // compute_exposed_areas_by_sampling(geo, scene, hitAccumulator, geometryID, *rng, *rngSeed1);
        // } else {
        //   use_entire_areas_of_primitives_as_exposed(geo, hitAccumulator);
        // }
        auto discareas = compute_disc_areas(geo, mBoundary);
        hitAccumulator.set_exposed_areas(discareas);
      }
      // Assertion: hitAccumulator is reduced to one instance by openmp reduction

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitAccumulator = std::make_unique<rti::trace::hit_accumulator<numeric_type> >(hitAccumulator);
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      rtcReleaseGeometry(geometry);
      rtcReleaseGeometry(boundary);

      auto raylog = RAYLOG_GET_PTR();
      if (raylog != nullptr) {
        auto raylogfilename = "raylog.vtp";
        std::cout << "Writing ray log to " << raylogfilename << std::endl;
        rti::io::vtp_writer<float>::write(raylog, raylogfilename);
      }
      auto raysrclog = RAYSRCLOG_GET_PTR();
      if (raysrclog != nullptr) {
        auto raysrclogfilename = "raysrclog.vtp";
        std::cout << "Writing ray src log to " << raysrclogfilename << std::endl;
        rti::io::vtp_writer<float>::write(raysrclog, raysrclogfilename);
      }

      { // Debug
        std::cout << "[Alex] V 4" << std::endl;
        auto dirpath = std::string {"/home/alexanders/vtk/outputs/current/"};
        auto fileidx = 0u;
        auto file1path = std::string {};
        while (true) {
          file1path = dirpath + "result-file-" + std::to_string(fileidx) + ".vtp";
          if (rti::util::file_exists(file1path)) {
            fileidx += 1;
            continue;
          }
          break;
        }
        auto file2path = dirpath + "bounding-box-" + std::to_string(fileidx) + ".vtp";
        assert( ! rti::util::file_exists(file1path) && ! rti::util::file_exists(file2path) && "Assumption");
        std::cout << "Writing file " << file1path << std::endl;
        auto metadata =  std::vector<rti::util::pair<std::string> > {};
        rti::io::vtp_writer<float>::write(geo, hitAccumulator, file1path, metadata);
        std::cout << "Writing bounding box to " << file2path << std::endl;
        rti::io::vtp_writer<numeric_type>::write(mBoundary, file2path);

        
        // auto time = std::time(nullptr);
        // auto loctime = *std::localtime(&time);
        // // auto strstr = std::istringstream {};
        // // strstr << std::put_time(&loctime, "%F_%H-%M-%S");
        // char timestr [64];
        // // timestr.reserve(32); // magic number
        // // std::cout << timestr.capacity() << std::endl;
        // auto success = strftime(&timestr[0], 64, "%F_%H-%M-%S", &loctime);
        // if ( ! success) {
        //   std::cout << "Error creating file name" << std::endl;
        // }
        
        // std::cout << "[Alex] timestr == " << timestr << std::endl;
        // auto dirpath = std::string {"/home/alexanders/vtk/outputs/"};
        // // if (filesystem::exists(dirpath)) {
        // //   std::cerr << "Error unique output directory " << dirpath
        // //             << " for Alex's VTK files already exists." << std::endl;
        // // }
        // // fs::create_directories(dirpath);
        
        // // auto path = "/home/alexanders/vtk/outputs/result-file.vtp";
        // auto path = dirpath + "result-file-" + timestr + ".vtp";
        // if (rti::util::file_exists(path)) {
        //   std::cerr << "Error: file " << path << " already exists" << std::endl;
        // } else {
        //   std::cout << "Writing file " << path << std::endl;
        //   auto metadata =  std::vector<rti::util::pair<std::string> > {};
        //   rti::io::vtp_writer<float>::write(geo, hitAccumulator, path, metadata);
        // }
        
        // path = dirpath + "bounding-box-" + timestr + ".vtp";
        // if (rti::util::file_exists(path)) {
        //   std::cerr << "Error: file " << path << " already exists" << std::endl;
        // } else {
        //   // auto path = "/home/alexanders/vtk/outputs/bounding-box.vtp";
        //   std::cout << "Writing bounding box to " << path << std::endl;
        //   rti::io::vtp_writer<numeric_type>::write(mBoundary, path);
        // }
      }

      return result;
    }

  private:

    std::vector<numeric_type>
    compute_disc_areas
    (rti::geo::absc_point_cloud_geometry<numeric_type>& geo,
     rti::geo::i_boundary<numeric_type>& boundary)
    {
      auto xmin = std::numeric_limits<numeric_type>::max();
      auto ymin = std::numeric_limits<numeric_type>::max();
      auto xmax = std::numeric_limits<numeric_type>::lowest();
      auto ymax = std::numeric_limits<numeric_type>::lowest();
      for (auto const& vv : boundary.get_vertices()) {
        if (vv[0] < xmin) xmin = vv[0];
        if (vv[0] > xmax) xmax = vv[0];
        if (vv[1] < ymin) ymin = vv[1];
        if (vv[1] > ymax) ymax = vv[1];
      }
      for (auto const& vv : boundary.get_vertices()) {
        assert( (vv[0] == xmin || vv[0] == xmax) && "Bounding box assumption");
        assert( (vv[1] == ymin || vv[1] == ymax) && "Bounding box assumption");
      }
      auto dbbi = rti::geo::disc_bounding_box_intersector(xmin, ymin, xmax, ymax);
      auto numofprimitives = geo.get_num_primitives();
      auto areas = std::vector<numeric_type> (numofprimitives, 0);
      #pragma omp for
      for (size_t idx = 0; idx < numofprimitives; ++idx) {
        areas[idx] = dbbi.area_inside(geo.get_prim_ref(idx), geo.get_normal_ref(idx));
      }
      return areas;
    }

    
  private:
    
    rti::geo::i_factory<numeric_type>& mFactory;
    rti::geo::i_boundary<numeric_type>& mBoundary;
    rti::ray::i_source& mSource;
    size_t mNumRays;
    rti::particle::i_particle_factory<numeric_type>& particlefactory;
  };
}}
