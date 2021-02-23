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

#include "dummy_counter.hpp"
#include "hit_accumulator.hpp"
#include "local_intersector.hpp"
//#include "point_cloud_context.hpp"
#include "result.hpp"
//#include "../geo/absc_point_cloud_geometry.hpp"
#include "../geo/disc_bounding_box_intersector.hpp"
#include "../geo/boundary_x_y.hpp"
#include "../geo/absc_geometry.hpp"
#include "../io/vtp_writer.hpp"
#include "../mc/rejection_control.hpp"
#include "../particle/i_particle.hpp"
//#include "../ray/constant_direction.hpp"
#include "../ray/disc_origin.hpp"
#include "../ray/i_source.hpp"
#include "../reflection/i_reflection.hpp"
#include "../rng/mt64_rng.hpp"
#include "../util/logger.hpp"
#include "../util/ray_logger.hpp"
#include "../util/timer.hpp"

namespace rti { namespace trace {
    
    template<typename numeric_type, typename particle_type, typename reflection_type>
  class tracer {
    
    static_assert(std::is_base_of<particle::i_particle<numeric_type>, particle_type>::value, "Precondition");
    static_assert(std::is_base_of<reflection::i_reflection<numeric_type>, reflection_type>::value, "Precondition");
    
  public:
    
    tracer
    (geo::point_cloud_disc_geometry<numeric_type>& pGeometry,
     geo::boundary_x_y<numeric_type>& pBoundary,
     ray::i_source& pSource,
     size_t pNumRays) :
      mGeometry(pGeometry),
      mBoundary(pBoundary),
      mSource(pSource),
      mNumRays(pNumRays)
    {
      assert(mGeometry.get_rtc_device() == pBoundary.get_rtc_device() &&
              "the geometry and the boundary need to refer to the same Embree (rtc) device");
      auto& rtcdevice = mGeometry.get_rtc_device();
      assert(rtcGetDeviceProperty(rtcdevice, RTC_DEVICE_PROPERTY_VERSION) >= 30601 &&
             "Error: The minimum version of Embree is 3.6.1");
      RLOG_WARNING << "Warning: tnear set to a constant! FIX" << std::endl;
    }

    trace::result<numeric_type> run()
    {
      // Prepare a data structure for the result.
      auto result = trace::result<numeric_type> {};
      result.inputFilePath = mGeometry.get_input_file_path();
      result.geometryClassName = typeid(mGeometry).name();

      // Prepare Embree
      auto rtcdevice = mGeometry.get_rtc_device();
      auto rtcscene = rtcNewScene(rtcdevice);

      // scene flags
      // Does one need to set this flag if one uses the registered call-back functions only?
      //rtcSetSceneFlags(rtcscene, RTC_SCENE_FLAG_NONE | RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
      rtcSetSceneFlags(rtcscene, RTC_SCENE_FLAG_NONE);
      //std::cerr << "rtc get scene flags == " << rtcGetSceneFlags(rtcscene) << std::endl;

      // Selecting higher build quality results in better rendering performance but slower
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      auto bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(rtcscene, bbquality);
      auto rtcgeometry = mGeometry.get_rtc_geometry();
      auto rtcboundary = mBoundary.get_rtc_geometry();

      auto boundaryID = rtcAttachGeometry(rtcscene, rtcboundary);
      auto geometryID = rtcAttachGeometry(rtcscene, rtcgeometry);

      assert(rtcGetDeviceError(rtcdevice) == RTC_ERROR_NONE && "Error");

      // Use openMP for parallelization
      #pragma omp parallel
      {
        rtcJoinCommitScene(rtcscene);
        // TODO: move to the other parallel region at the bottom
      }

      result.numRays = mNumRays;

      auto boundaryReflection = reflection::specular<numeric_type> {};

      auto geohitc = 0ull;
      auto nongeohitc = 0ull;
      auto hitAccumulator = trace::hit_accumulator<numeric_type> {mGeometry.get_num_primitives()};

      #pragma omp declare \
        reduction(hit_accumulator_combine : \
                  trace::hit_accumulator<numeric_type> : \
                  omp_out = trace::hit_accumulator<numeric_type>(omp_out, omp_in)) \
        initializer(omp_priv = trace::hit_accumulator<numeric_type>(omp_orig))

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      // auto rng = std::make_unique<rng::cstdlib_rng>();
      auto rng = rng::mt64_rng {};

      // Start timing
      auto timer = util::timer {};

      #pragma omp parallel \
        reduction(+ : geohitc, nongeohitc) \
        reduction(hit_accumulator_combine : hitAccumulator)
      {
        // Thread local data goes here, if it is not needed anymore after the execution
        // of the parallel region.
        alignas(128) auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        // 5393
        // 115249
        // 2147483647
        // 1442968193
        auto seed = (unsigned int) ((omp_get_thread_num() + 1) *  31); // multiply by magic number (prime)
        auto rngstate1 = rng::mt64_rng::state { seed + 0 };
        auto rngstate2 = rng::mt64_rng::state { seed + 1 };
        auto rngstate3 = rng::mt64_rng::state { seed + 2 };
        auto rngstate4 = rng::mt64_rng::state { seed + 3 };
        auto rngstate5 = rng::mt64_rng::state { seed + 4 };
        auto rngstate6 = rng::mt64_rng::state { seed + 5 };
        auto rngstate7 = rng::mt64_rng::state { seed + 6 };

        // A dummy counter for the boundary
        auto boundaryCntr = trace::dummy_counter {};

        // thread-local particle and reflection object
        auto particle = particle_type {};
        auto surfreflect = reflection_type {};
          
        // probabilistic weight
        auto rayweight = (numeric_type) 1;

        auto rtccontext = RTCIntersectContext {};
        rtcInitIntersectContext(&rtccontext);

        size_t progresscnt = 0;

        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) mNumRays; ++idx) {      
          particle.init_new();
          rayweight = get_init_ray_weight();
          auto lastinitRW = rayweight;
          mSource.fill_ray(rayhit.ray, rng, rngstate1, rngstate2, rngstate3, rngstate4); // fills also tnear
          RAYSRCLOG(rayhit);
          if_RLOG_PROGRESS_is_set_print_progress(progresscnt, mNumRays);
          auto reflect = false;
          do {
            RLOG_DEBUG
              << "preparing ray == ("
              << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z
              << ") ("
              << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z
              << ")" << std::endl;
            rayhit.ray.tfar = std::numeric_limits<float>::max(); // Embree uses float
            rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rayhit.ray.tnear = 1e-4; // tnear is also set in the particle source
            // Run the intersection
            rtcIntersect1(rtcscene, &rtccontext, &rayhit);

            RAYLOG(rayhit, rayhit.ray.tfar);
            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
              // No  hit
              nongeohitc += 1;
              reflect = false;
              RLOG_TRACE << "i";
              break; // break do-while loop
            }
            // else
            // A hit
            if (rayhit.hit.geomID == boundaryID) {
              // // Ray hit the boundary
              // reflect = true;
              // auto orgdir = boundaryReflection.use (rayhit.ray, rayhit.hit, mBoundary, rng, rngstate2);
              
              reflect = true;
              auto orgdir = mBoundary.process_hit(rayhit.ray, rayhit.hit);
              // TODO: optimize
              rayhit.ray.org_x = orgdir[0][0];
              rayhit.ray.org_y = orgdir[0][1];
              rayhit.ray.org_z = orgdir[0][2];
              rayhit.ray.dir_x = orgdir[1][0];
              rayhit.ray.dir_y = orgdir[1][1];
              rayhit.ray.dir_z = orgdir[1][2];
              RLOG_TRACE << "b";
              continue;
            } 
            assert (rayhit.hit.geomID == geometryID && "Correctness Assumption");

            // If the dot product of the ray direction and the surface normal is greater than zero, then
            // we hit the back face of the disc.
            auto const& ray = rayhit.ray;
            auto const& hit = rayhit.hit;
            if (rti::util::dot_product(rti::util::triple<numeric_type> {ray.dir_x, ray.dir_y, ray.dir_z},
                                       mGeometry.get_normal(hit.primID)) > 0) {
              // Hit from the back
              RLOG_TRACE << "a";
              // Let ray through, i.e., continue.
              reflect = true; // reflect means continue
              rayhit.ray.org_x = ray.org_x + ray.dir_x * ray.tfar;
              rayhit.ray.org_y = ray.org_y + ray.dir_y * ray.tfar;
              rayhit.ray.org_z = ray.org_z + ray.dir_z * ray.tfar;
              // keep ray direction as it is
              continue;
            }
            RLOG_TRACE << "h";
            geohitc += 1;
            RLOG_DEBUG << "rayhit.hit.primID == " << rayhit.hit.primID << std::endl;
            RLOG_DEBUG << "prim == " << mGeometry.prim_to_string(rayhit.hit.primID) << std::endl;
            // auto sticking = particle.process_hit(hit.primID, {ray.dir_x, ray.dir_y, ray.dir_z});
            auto sticking = particle.get_sticking_probability(rayhit.ray, rayhit.hit, mGeometry, rng, rngstate5);
            auto valuetodrop = rayweight * sticking;
            hitAccumulator.use(rayhit.hit.primID, valuetodrop);
            check_for_additional_intersections(rayhit.ray, rayhit.hit.primID, hitAccumulator, valuetodrop);
            rayweight -= valuetodrop;
            if (rayweight == 0) {
              break;
            }
            reflect = mc::rejection_control<numeric_type>::check_weight_reweight_or_kill
              (rayweight, lastinitRW, rng, rngstate6);
            if ( ! reflect ) {
              break;
            }
            auto orgdir = surfreflect.use (rayhit.ray, rayhit.hit, mGeometry, rng, rngstate7);
            // TODO: optimize
            rayhit.ray.org_x = orgdir[0][0];
            rayhit.ray.org_y = orgdir[0][1];
            rayhit.ray.org_z = orgdir[0][2];
            rayhit.ray.dir_x = orgdir[1][0];
            rayhit.ray.dir_y = orgdir[1][1];
            rayhit.ray.dir_z = orgdir[1][2];

            // // ATTENTION tnear is set in another function, too! When the ray starts from the source, then
            // // the source class also sets tnear!
            // auto tnear = 1e-4f; // float
            // // Same holds for time
            // auto time = 0.0f; // float
            // // reinterpret_cast<__m128&>(rayhit.ray) = _mm_load_ps(vara);
            // // reinterpret_cast<__m128&>(rayhit.ray.dir_x) = _mm_load_ps(varb);

            // // For the store-forwarding stall also see the following link:
            // // https://stackoverflow.com/questions/49265634/what-is-the-difference-between-loadu-ps-and-set-ps-when-using-unformatted-data
            // // OPTIMIZATION: It might be better to change the memory layout of the rayout variable
            // // such that tnear (a constant) is saved right in front of rayout. Then we can use two
            // // __m128d _mm_store_sd (__m128d a, __m128d b) to write 128 bits in two chunks of 64 bits
            // // into the destination.
            // reinterpret_cast<__m128&>(rayhit.ray) =
            //   _mm_set_ps(tnear,
            //              (float) rtiContext->rayout[0][2],
            //              (float) rtiContext->rayout[0][1],
            //              (float) rtiContext->rayout[0][0]);
            // reinterpret_cast<__m128&>(rayhit.ray.dir_x) =
            //   _mm_set_ps(time,
            //              (float) rtiContext->rayout[1][2],
            //              (float) rtiContext->rayout[1][1],
            //              (float) rtiContext->rayout[1][0]);
          } while (reflect);
        }
        auto discareas = compute_disc_areas(mGeometry, mBoundary);
        hitAccumulator.set_exposed_areas(discareas);
      }
      // Assertion: hitAccumulator is reduced to one instance by openmp reduction

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitAccumulator = std::make_unique<trace::hit_accumulator<numeric_type> >(hitAccumulator);
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      rtcReleaseGeometry(rtcgeometry);
      rtcReleaseGeometry(rtcboundary);

      auto raylog = RAYLOG_GET_PTR();
      if (raylog != nullptr) {
        auto raylogfilename = "raylog.vtp";
        std::cout << "Writing ray log to " << raylogfilename << std::endl;
        io::vtp_writer<float>::write(raylog, raylogfilename);
      }
      auto raysrclog = RAYSRCLOG_GET_PTR();
      if (raysrclog != nullptr) {
        auto raysrclogfilename = "raysrclog.vtp";
        std::cout << "Writing ray src log to " << raysrclogfilename << std::endl;
        io::vtp_writer<float>::write(raysrclog, raysrclogfilename);
      }

      // { // Debug
      //   std::cout << "[Alex] V 5" << std::endl;
      //   auto dirpath = std::string {"/home/alexanders/vtk/outputs/current/"};
      //   auto fileidx = 0u;
      //   auto file1path = std::string {};
      //   while (true) {
      //     file1path = dirpath + "result-file-" + std::to_string(fileidx) + ".vtp";
      //     if (util::file_exists(file1path)) {
      //       fileidx += 1;
      //       continue;
      //     }
      //     break;
      //   }
      //   auto file2path = dirpath + "bounding-box-" + std::to_string(fileidx) + ".vtp";
      //   assert( ! util::file_exists(file1path) && ! util::file_exists(file2path) && "Assumption");
      //   std::cout << "Writing file " << file1path << std::endl;
      //   auto metadata =  std::vector<util::pair<std::string> > {};
      //   io::vtp_writer<float>::write(mGeometry, hitAccumulator, file1path, metadata);
      //   std::cout << "Writing bounding box to " << file2path << std::endl;
      //   io::vtp_writer<numeric_type>::write(mBoundary, file2path);
      // }

      return result;
    }

  private:

    constexpr numeric_type get_init_ray_weight()
    {
      return 1;
    }

    void check_for_additional_intersections
    (RTCRay& ray,
     unsigned int hit1id,
     trace::hit_accumulator<numeric_type>& hitAcc,
     numeric_type valuetodrop)
    {
      // std::cout << "neighborhoodsize == " << mGeometry.get_neighbors(hit1id).size() << std::endl;

      // std::cout << "disc " << hit1id << " == " << disc[0] << " " << disc[1] << " " << disc[2] << " " << disc[3] << std::endl;

      // { // Debug
      //   std::cerr << "## Debug" << std::endl;
      //   local_intersector::intersect(ray, disc, mGeometry.get_normal_ref(hit1id));
      //   std::cerr << "## Debug" << std::endl;
      // }

      // { // Debug
      //   std::cout << "check_for_additional_intersections(): " << hit1id << " ";
      // }
      
      auto cnt = 0u;
      for (auto const& id : mGeometry.get_neighbors(hit1id)) {
        // std::cout << "using prim id " << id << std::endl;
        // auto& printdisc = mGeometry.get_prim_ref(id);
        // std::cout
        //   << "printdisc == " << printdisc[0] << " " << printdisc[1]
        //   << " " << printdisc[2] << " " << printdisc[3] << std::endl;

        auto const& disc = mGeometry.get_prim_ref(id);
        auto const& dnormal = mGeometry.get_normal_ref(id);
        auto intersects = local_intersector::intersect(ray, disc, dnormal);
        if ( intersects ) {
          hitAcc.use(id, valuetodrop);
          cnt += 1;
          // { // Debug
          //   std::cout << id << " ";
          // }

        }
      }
      // std::cout << " hits == " << cnt << std::endl;
      // { // Debug
      //   std::cout << std::endl;
      // }
    }
      
    std::vector<numeric_type>
    compute_disc_areas
    (geo::point_cloud_disc_geometry<numeric_type>& geometry,
     geo::boundary_x_y<numeric_type>& boundary)
    {
      // TODO: Instead of comparing introduce a method boundary.get_bounding_box()
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
      auto dbbi = geo::disc_bounding_box_intersector(xmin, ymin, xmax, ymax);
      auto numofprimitives = geometry.get_num_primitives();
      auto areas = std::vector<numeric_type> (numofprimitives, 0);
      #pragma omp for
      for (size_t idx = 0; idx < numofprimitives; ++idx) {
        areas[idx] = dbbi.area_inside(geometry.get_prim_ref(idx), geometry.get_normal_ref(idx));
      }
      return areas;
    }

    void if_RLOG_PROGRESS_is_set_print_progress(size_t& raycnt, size_t const& totalnumrays)
    {
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
          std::string(filllength, fillsymbol) +
          std::string(std::max(0, (int) barlength - (int) filllength), emptysymbol) +
          std::string(1, barendsymbol) + " " + percentagestring ;
        RLOG_PROGRESS << "\r" << bar;
        if (filllength >= barlength) {
          RLOG_PROGRESS << std::endl;
        }
      }
      raycnt += 1;
    }
    
  private:
    
    geo::point_cloud_disc_geometry<numeric_type>& mGeometry;
    geo::boundary_x_y<numeric_type>& mBoundary;
    ray::i_source& mSource;
    size_t mNumRays;
  };
}}
