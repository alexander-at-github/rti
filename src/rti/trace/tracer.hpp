#pragma once

//#define STATS_ENABLE_MATRIX_FEATURES
//#define STATS_ENABLE_STDVEC_WRAPPERS
#define STATS_ENABLE_ARMA_WRAPPERS
#define STATS_GO_INLINE
//#define STATS_USE_OPENMP

#include <boost/core/demangle.hpp>
#include <boost/current_function.hpp>

#include <cmath>
#include <chrono>
#include <omp.h>

#include <embree3/rtcore.h>
#include <RcppArmadillo.h>
#include <RInside.h>
#include <stats.hpp>

#include "rti/geo/absc_point_cloud_geometry.hpp"
#include "rti/geo/i_boundary.hpp"
#include "rti/geo/i_geometry.hpp"
#include "rti/io/vtp_writer.hpp"
#include "rti/particle/i_particle.hpp"
#include "rti/particle/i_particle_factory.hpp"
#include "rti/ray/constant_direction.hpp"
#include "rti/ray/disc_origin.hpp"
#include "rti/ray/i_source.hpp"
#include "rti/reflection/diffuse.hpp"
//#include "rti/reflection/specular.hpp"
#include "rti/rng/cstdlib_rng.hpp"
#include "rti/rng/mt64_rng.hpp"
//#include "rti/trace/counter.hpp"
#include "rti/trace/dummy_counter.hpp"
#include "rti/trace/hit_accumulator.hpp"
#include "rti/trace/point_cloud_context.hpp"
#include "rti/trace/result.hpp"
#include "rti/trace/triangle_context_simplified.hpp"
#include "rti/util/logger.hpp"
#include "rti/util/ray_logger.hpp"
#include "rti/util/timer.hpp"

namespace rti { namespace trace {
  template<typename numeric_type>
  class tracer {

  public:
    tracer(rti::geo::i_factory<numeric_type>& pFactory,
           rti::geo::i_boundary<numeric_type>& pBoundary,
           rti::ray::i_source& pSource,
           rti::particle::i_particle_factory<numeric_type>& particlefactory) :
      mFactory(pFactory),
      mBoundary(pBoundary),
      mSource(pSource),
      particlefactory(particlefactory),
      rtigeometry(mFactory.get_geometry()),
      rtcgeometry(rtigeometry.get_rtc_geometry()),
      rtcboundary(mBoundary.get_rtc_geometry()),
      device(rtigeometry.get_rtc_device()),
      scene(rtcNewScene(device))
    {
      assert (mFactory.get_geometry().get_rtc_device() == pBoundary.get_rtc_device() &&
              "the geometry and the boundary need to refer to the same Embree (rtc) device");
      assert(rtcGetDeviceProperty(mFactory.get_geometry().get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) != 0 &&
             "Error: Embree filter functions are not supported by your Embree instance.");
      // std::cerr
      //   << "RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED == "
      //   << rtcGetDeviceProperty(mFactory.get_geometry().get_rtc_device(), RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED)
      //   << std::endl;
      assert(rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_VERSION) >= 30601 && "Error: The minimum version of Embree is 3.6.1");
      if (rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED) != 0) {
        std::cerr << "=== Warning: Embree backface culling enabled. This may result in incorrect results "
                  << "(for triangles)" << std::endl;
        assert( false && "Error: backface culling is enabled; as a consequence for triangles \
                the tracer depends on the order of the vertices of the triangles");
      }
      RLOG_WARNING << "Warning: tnear set to a constant! FIX" << std::endl;

      prepare_embree();
    }

  private:
    void if_RLOG_PROGRESS_is_set_print_progress(size_t& raycnt, size_t totalnumrays) {
      if (omp_get_thread_num() != 0) {
        return;
      }
      auto barlength = 60u;
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
          "Monte Carlo Progress: " + std::string(1, barstartsymbol) +
          std::string(filllength, fillsymbol) + std::string(std::max(0, (int) barlength - (int) filllength), emptysymbol) +
          std::string(1, barendsymbol) + " " + percentagestring ;
        RLOG_PROGRESS << "\r" << bar;
        if (filllength >= barlength) {
          RLOG_PROGRESS << std::endl;
        }
      }
      raycnt += 1;
    }

    // void
    // compute_exposed_areas_by_sampling
    // (rti::geo::i_geometry<numeric_type>& geo,
    //  RTCScene& scene,
    //  rti::trace::i_hit_accumulator<numeric_type>& hitacc,
    //  unsigned int geometryID,
    //  rti::rng::i_rng& rng,
    //  rti::rng::i_rng::i_state& rngstate
    //  )
    // {
    //   assert(false && "Deprecated");

    //   // Precondition: Embree has been set up properly
    //   //   (e.g., the geometry has been built)
    //   // Precondition: We are in an OpenMP parallel block
    //   // Precondition: the following cast characterizes a precondition:
    //   auto pcgeo =
    //     dynamic_cast<rti::geo::absc_point_cloud_geometry<numeric_type>*> (&geo);

    //   std::cerr << "### FIX: This function does not consider that areas of "
    //             << "discs may be located outside of the boundary" << std::endl;

    //   // We use a new RTC context object.
    //   auto context = RTCIntersectContext {};
    //   rtcInitIntersectContext(&context);
    //   auto numOfSamples = 1024u;

    //   auto numofprimitives = pcgeo->get_num_primitives();
    //   auto areas = std::vector<double> (numofprimitives, 0);
    //   #pragma omp for
    //   for (size_t primidx = 0; primidx < numofprimitives; ++primidx) {
    //     auto pointradius = pcgeo->get_prim(primidx);
    //     auto normal = pcgeo->get_normal(primidx);
    //     auto point = rti::util::triple<numeric_type>
    //       {pointradius[0], pointradius[1], pointradius[2]};
    //     auto radius = pointradius[3];
    //     auto origincenter = rti::util::sum(point, rti::util::scale(2*radius, normal));
    //     auto invnormal = rti::util::inv(normal);
    //     auto origin = rti::ray::disc_origin<numeric_type>
    //       {origincenter, invnormal, radius};
    //     auto direction = rti::ray::constant_direction<numeric_type> {invnormal};
    //     auto source = rti::ray::source<numeric_type> {origin, direction};

    //     alignas(128) auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //     auto hits = 0u;
    //     for (size_t sidx = 0; sidx < numOfSamples; ++sidx) {
    //       rayhit.ray.tnear = 0;
    //       rayhit.ray.time = 0;
    //       rayhit.ray.tfar = std::numeric_limits<float>::max();
    //       rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    //       rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    //       source.fill_ray(rayhit.ray, rng, rngstate);

    //       //RAYSRCLOG(rayhit);

    //       rtcIntersect1(scene, &context, &rayhit);

    //       //RAYLOG(rayhit, rayhit.ray.tfar);

    //       if (rayhit.hit.geomID != geometryID)
    //         continue;
    //       if (rayhit.hit.primID != primidx)
    //         continue;
    //       hits += 1;
    //     }
    //     assert (0 <= hits && hits <= numOfSamples && "Correctness assertion");
    //     //std::cerr << " " << (double) hits / numOfSamples;
    //     areas[primidx] = pcgeo->get_area(primidx) * (double) hits / numOfSamples;
    //   }
    //   hitacc.set_exposed_areas(areas);
    // }

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

    rti::util::pair<double>
    compute_mean_coord
    (std::vector<std::pair<rti::util::triple<float>, double> > coordvaluepairs)
    {
      { // hack
        static auto firstexecution = true;
        if (firstexecution)
          std::cout
            << ">>> Warning: This code assumes a source perpendicular to the z-achsis!" << std::endl;
        firstexecution = false;
      }
      auto acc = rti::util::pair<double> {0, 0}; // x and y
      for (auto const& cv : coordvaluepairs) {
        acc[0] += cv.first[0];
        acc[1] += cv.first[1];
      }
      auto num = coordvaluepairs.size();
      return {acc[0]/num, acc[1]/num};
    }

    rti::util::pair<double>
    compute_variance_of_coords
    (std::vector<std::pair<rti::util::triple<float>, double> > coordvaluepairs,
     rti::util::triple<double> meancoord)
    {
      { // hack
        static auto firstexecution = true;
        if (firstexecution)
          std::cout
            << ">>> Warning: This code assumes a source perpendicular to the z-achsis!" << std::endl;
        firstexecution = false;
      }
      auto acc = rti::util::pair<double> {0, 0};
      for (auto const& cv : coordvaluepairs) {
        auto tmp0 = (cv.first[0] - meancoord[0]);
        acc[0] += tmp0 * tmp0;
        auto tmp1 = (cv.first[1] - meancoord[1]);
        acc[1] += tmp1 * tmp1;
      }
      auto num = coordvaluepairs.size();
      return {acc[0]/num, acc[1]/num};
    }

    std::vector<rti::util::tripleN<double, rti::util::pair<double>, rti::util::pair<double> > >
    //std::vector<std::tuple<double, rti::util::pair<double>, rti::util::pair<double> > >
    //rti::util::pair<std::vector<rti::util::pair<double> > >
    compute_GMM_using_R
    (std::vector<std::pair<rti::util::triple<float>, double> > sourcesamples)
    {
      auto argc = 0;
      auto argv = (char**) nullptr;
      auto rplotfilename1 = "r-plot_1.eps";
      auto rplotfilename2 = "r-plot_2.eps";
      auto ri = RInside (argc, argv);
      assert (sourcesamples.size() <= std::numeric_limits<int>::max() && "Correctness Assumption");
      auto nrows = (int) sourcesamples.size();
      auto ncols = 2;
      auto data = Rcpp::NumericMatrix {nrows, ncols};
      for (size_t idx = 0; idx < nrows; ++idx) {
        data(idx, 0) = sourcesamples[idx].first[0]; // 0 === x
        data(idx, 1) = sourcesamples[idx].first[1]; // 1 === y
        //std::cout << data(idx, 0) << " " << data(idx, 1) << std::endl;
      }
      ri["data"] = data;
      //ri.parseEvalQ("print('FOO FROM R'); print(data)");
      //std::string cmd0 = "library(mclust); density <- densityMclust(data); summary(density, parameters=TRUE);";
      //ri.parseEvalQ(cmd0);
      // std::string mclustfitcmd = "library(mclust); density <- densityMclust(data, prior=priorControl(), G=1, modelNames='VVI');";

      auto modelNames = (std::string) "VVI";

      //ri.parseEvalQ("print(\"just printing something\")");
      ri.parseEvalQ("library(mclust)");
      ri.parseEvalQ(std::string("mclusticl <- mclustICL(data, modelNames='") + modelNames + "')");
      //ri.parseEvalQ("print(mclusticl[1])"); // R uses 1-based index
      std::cout << "HERE" << std::endl;
      auto rrr = Rcpp::as<std::vector<double> > (ri["mclusticl"]);
      { // Debug
        std::cout << "rrr: ";
        for (auto const& rrrelem : rrr) {
          std::cout << rrrelem << " ";
        }
        std::cout << std::endl;
      }
      auto idxofmax = std::max_element(rrr.begin(), rrr.end()) - rrr.begin();
      std::cout << "index of max element == " << idxofmax << std::endl;

      auto numModelMixComponents = idxofmax + 1; // +1 because of 0-based indexing
      std::cout << "Using " << numModelMixComponents << " bivariate Gaussians" << std::endl;
      auto densitycmdstr = std::string("density <- densityMclust(data, modelNames='") +
        modelNames + "', G=" + std::to_string(numModelMixComponents) + ", prior=priorControl())";
      std::cout << densitycmdstr << std::endl;
      ri.parseEvalQ(densitycmdstr);

      std::cout << std::endl;
      ri.parseEvalQ("print(summary(density$parameters))");
      std::cout << std::endl;
      ri.parseEvalQ("print(density$parameters$pro)");
      std::cout << std::endl;
      ri.parseEvalQ("print(density$parameters$mean)");
      std::cout << std::endl;
      ri.parseEvalQ("print(density$parameters$variance$sigma)");

      std::cout << "printing GMM to file " << rplotfilename1 << " and " <<rplotfilename2 << std::endl;
      std::cout << "for the plot I am using hard-coded limits of the finstack geoemtry" << std::endl;
      ri.parseEvalQ(
         std::string("postscript('") + rplotfilename1 + "', horizontal=F, width=8, height=8, paper='special', onefile=F); \
         plot(density, what='density', type='persp', xlim=c(-0.004501, 0.328006), ylim=c(0.0005, 0.2875));  \
         dev.off();");
      ri.parseEvalQ(
         std::string("postscript('") + rplotfilename2 + "', horizontal=F, width=8, height=8, paper='special', onefile=F); \
         plot(density, what='density', data=data, xlim=c(-0.004501, 0.328006), ylim=c(0.0005, 0.2875));  \
         dev.off();");

      auto result =
        std::vector<rti::util::tripleN<double, rti::util::pair<double>, rti::util::pair<double> > > (numModelMixComponents);

      std::cout << "HERE 0" << std::endl;
      auto componentprobabilities = Rcpp::as<std::vector<double> > (ri.parseEval("density$parameters$pro"));
      assert(componentprobabilities.size() == numModelMixComponents && "Correctness Assertion");

      std::cout << "HERE 1" << std::endl;
      auto means = Rcpp::as<arma::mat> (ri.parseEval("density$parameters$mean"));
      assert(means.n_cols == numModelMixComponents && "Correctness Assertion");
      assert(means.n_rows == 2 && "Assumption: Model has two dimensions");
      // auto meansvec = std::vector<rti::util::pair<double> > (means.n_cols);
      // for (size_t idxc = 0; idxc < means.n_cols; ++idxc) {
      //   meansvec[idxc] = {means(0, idxc), means(1, idxc)};
      // }

      std::cout << "HERE 2" << std::endl;
      auto variances = Rcpp::as<arma::cube> (ri.parseEval("density$parameters$variance$sigma"));
      //auto variances = Rcpp::as<arma::mat> (ri.parseEval("density$parameters$variance$sigma[,,1]")); // also works
      // Note that the ri[] syntax (as shown in the commented-out code below) does not work.
      //auto variances = Rcpp::as<arma::mat> (ri["density$parameters$variance$sigma[,,1]"]); // does not work!
      std::cout << "HERE 3" << std::endl;
      assert(variances.n_slices == numModelMixComponents && "Correctness Assertion");
      assert(variances.n_rows == 2 && variances.n_cols == 2 && "Assumption: Model has two dimensions");
      assert(variances.n_rows == variances.n_cols && "Correctness Assertion");
      // auto variancesvec = std::vector<rti::util::pair<double> > (variances.n_slices);
      // // for (size_t idxr = 0; idxr < variances.n_rows; ++idxr) {
      // //   for (size_t idxc = 0; idxc < variances.n_cols; ++idxc) {
      // //     if (idxr == idxc) {
      // //       variancesvec[idxr] = variances(idxr, idxc, 0);
      // //       continue;
      // //     }
      // //     assert(variances(idxr, idxc, 0) == 0 && "Correctness Assertion");
      // //   }
      // // }
      // for (size_t idxs = 0; idxs < variances.n_slices; ++idxs) {
      //   assert(variances(0, 1, idxs) == 0 && variances(1, 0, idxs) == 0 && "Correctness Assertion");
      //   variancesvec[idxs] = {variances(0, 0, idxs), variances(1, 1, idxs)};
      // }
      for(size_t idx = 0; idx < numModelMixComponents; ++idx) {
        assert(variances(0, 1, idx) == 0 && variances(1, 0, idx) == 0 && "Correctness Assertion; VVI model");
        result[idx] = {componentprobabilities[idx],
                       {means(0, idx), means(1, idx)},
                       // Note: We are saving the variance here. Not the standatd deviation.
                       {variances(0, 0, idx), variances(1, 1, idx)}};
      }

      std::cerr << "DONE WITH THE R STUFF" << std::endl;
      return result;
      // return {meansvec, variancesvec};
      //return {std::vector<double> {}, std::vector<double> {}};

      // std::string mclustfitcmd = "library(mclust); density <- densityMclust(data, prior=priorControl(), modelNames='VVI');";
      // std::string getmeancmd = "density$parameters$mean;";
      // std::string getvariance1cmd = "density$parameters$variance$sigma[,,1]";
      // ri.parseEvalQ(mclustfitcmd);
      // auto mean = Rcpp::as<arma::mat>(ri.parseEval(getmeancmd));
      // auto variance1 = Rcpp::as<arma::mat>(ri.parseEval(getvariance1cmd));
      // std::cerr
      //   << "mean:" << std::endl << mean << std::endl
      //   << "variance1:" << std::endl << variance1 << std::endl;
      // // plot(density, what='density', data=data);
      // // plot(density, what='denstiy', type='persp')
      // ri.parseEvalQ(
      //    std::string("postscript('") + rplotfilename + "', horizontal=F, width=4, height=4,paper='special', onefile=F); \
      //    plot(density, what='density', type='persp'); \
      //    dev.off();");

      // std::string cmd = "set.seed(42); matrix(rnorm(9),3,3)"; 	// create a random Matrix in r
      // arma::mat m = Rcpp::as<arma::mat>(ri.parseEval(cmd)); // parse, eval + return result
      // std::cout << "callR() arma::mat: " << m(0, 0) << std::endl;
    }

  public:

    rti::trace::result<numeric_type> run_plain(size_t numrays)
    {
      try {
        return run_plain_aux(numrays);
      } catch (std::exception& ex) {
        std::cerr << ">>> exception: " << ex.what() << std::endl;
      }
      return {};
    }

    template<typename parameter_type>
    rti::trace::result<numeric_type>
    run_adaptive
    (std::vector<parameter_type> relativeerrors, size_t numrays)
    {
      auto rethreshold = (parameter_type) 0.2;
      auto numprerays = (size_t) (128 * 1024);
      std::cout << "Adaptive source sampling: using " << rethreshold << " as threshold on the relative error." << std::endl;
      try {
        return run_adaptive_aux(relativeerrors, rethreshold, numprerays, numrays);
      } catch (std::exception& ex) {
        std::cerr << ">>> exception: " << ex.what() << std::endl;
      }
      return {};
    }

  private:
    void prepare_embree()
    {
      // scene flags
      // Does one need to set this flag if one uses the registered call-back functions only?
      //rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE | RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
      // Selecting higher build quality results in better rendering performance but slower
      // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
      auto bbquality = RTC_BUILD_QUALITY_HIGH;
      rtcSetSceneBuildQuality(scene, bbquality);
      boundaryID = rtcAttachGeometry(scene, rtcboundary);
      geometryID = rtcAttachGeometry(scene, rtcgeometry);

      mFactory.register_intersect_filter_funs(mBoundary);
      assert(rtcGetDeviceError(device) == RTC_ERROR_NONE && "Error");
      // Use openMP for parallelization
      #pragma omp parallel
      {
        rtcJoinCommitScene(scene);
        // TODO: move to the other parallel region at the bottom
      }
    }

  public:
    void destroy_data()
    {
      rtcReleaseGeometry(rtcgeometry);
      rtcReleaseGeometry(rtcboundary);
    }

  private:
    rti::trace::result<numeric_type>
    run_plain_aux
    (size_t numrays)
    {
      // Prepare a data structure for the result.
      auto result = rti::trace::result<numeric_type> {};
      result.inputFilePath = rtigeometry.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(rtigeometry).name());
      result.tracerFunctionName = BOOST_CURRENT_FUNCTION;

      result.numRays = numrays;

      auto reflectionModel = rti::reflection::diffuse<numeric_type> {};
      auto boundaryReflection = rti::reflection::specular<numeric_type> {};

      auto geohitc = 0ull;
      auto nongeohitc = 0ull;
      auto hitAccumulator = rti::trace::hit_accumulator<numeric_type> {rtigeometry.get_num_primitives()};

      #pragma omp declare \
        reduction(hit_accumulator_combine : \
                  rti::trace::hit_accumulator<numeric_type> : \
                  omp_out = rti::trace::hit_accumulator<numeric_type>(omp_out, omp_in)) \
        initializer(omp_priv = rti::trace::hit_accumulator<numeric_type>(omp_orig))

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      // auto rng = std::make_unique<rti::rng::cstdlib_rng>();
      auto rng = std::make_unique<rti::rng::mt64_rng>();

      std::cout << "Starting stop watch" << std::endl;
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
        auto rtiContext = mFactory.get_new_context(geometryID, rtigeometry, reflectionModel,
                                                   hitAccumulator, boundaryID, mBoundary,
                                                   boundaryReflection, *rng, *rngSeed2,
                                                   *particle);

        // Initialize (also takes care for the initialization of the Embree context)
        rtiContext->init();
        rtiContext->set_initial_ray_weight(1); // standard weight for a ray

        std::cerr << "Starting " << numrays << " rays in " << BOOST_CURRENT_FUNCTION << std::endl;

        auto raycnt = (size_t) 0;
        auto rejectedsamples = (size_t) 0;
        auto cnt = 0u;
        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) numrays; ++idx) {

          // Note: Embrees backface culling does not solve our problem of intersections
          // when starting a new ray very close to or a tiny bit below the surface.
          // For that reason we set tnear to some value.
          // There is a risk though: when setting tnear to some strictly positive value
          // we depend on the length of the direction vector of the ray (rayhit.ray.dir_X).

          particle->init_new();

          RLOG_DEBUG << "NEW: Preparing new ray from source" << std::endl;
          // TODO: FIX: tnear set to a constant
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed1); // fills also tnear!
          rtiContext->init_ray_weight();

          //RAYSRCLOG(rayhit);

          if_RLOG_PROGRESS_is_set_print_progress(raycnt, numrays);

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

            //RAYLOG(rayhit, rtiContext->tfar);

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
        RLOG_TRACE << "before area calc" << std::endl << std::flush;
        if (rtiContext->compute_exposed_areas_by_sampling()) {
          std::cerr
            << "###############"
            << "### WARNING ### Computing exposed area by sampling does not work"
            << "###############" << std::endl;
          // // Embree Documentation: "Passing NULL as function pointer disables the registered callback function."
          // rtcSetGeometryIntersectFilterFunction(rtcgeometry, nullptr);
          // rtcSetGeometryIntersectFilterFunction(rtcboundary, nullptr);
          // rtcCommitGeometry(rtcgeometry);
          // rtcCommitGeometry(rtcboundary);
          // compute_exposed_areas_by_sampling(geo, scene, hitAccumulator, geometryID, *rng, *rngSeed1);
        } else {
          use_entire_areas_of_primitives_as_exposed(rtigeometry, hitAccumulator);
        }
        RLOG_TRACE << "after area calc" << std::endl << std::flush;
      }
      // Assertion: hitAccumulator is reduced to one instance by openmp reduction

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitAccumulator = std::make_unique<rti::trace::hit_accumulator<numeric_type> >(hitAccumulator);

      return result;
    }

    template<typename parameter_type>
    rti::trace::result<numeric_type>
    run_adaptive_aux
    (std::vector<parameter_type> relativeerrors, parameter_type rethreshold, size_t numprerays, size_t numrays)
    {
      // Prepare a data structure for the result.
      auto result = rti::trace::result<numeric_type> {};
      result.inputFilePath = rtigeometry.get_input_file_path();
      result.geometryClassName = boost::core::demangle(typeid(rtigeometry).name());
      result.tracerFunctionName = BOOST_CURRENT_FUNCTION;

      result.numRays = numrays;

      auto reflectionModel = rti::reflection::diffuse<numeric_type> {};
      auto boundaryReflection = rti::reflection::specular<numeric_type> {};

      auto geohitc = 0ull;
      auto nongeohitc = 0ull;
      auto hitAccumulator = rti::trace::hit_accumulator<numeric_type> {rtigeometry.get_num_primitives()};

      // vector of pairs of source sample points and "energy" delivered to the surface.
      auto relevantSourceSamples = std::vector<std::pair<rti::util::triple<float>, double> > {};
      auto maxrelevantsourcesamplessizeperthread = (size_t) (4 * 1024 / omp_get_max_threads()); // per thread

      #pragma omp declare \
        reduction(hit_accumulator_combine : \
                  rti::trace::hit_accumulator<numeric_type> : \
                  omp_out = rti::trace::hit_accumulator<numeric_type>(omp_out, omp_in)) \
        initializer(omp_priv = rti::trace::hit_accumulator<numeric_type>(omp_orig))
      #pragma omp declare \
        reduction(merge_vectors : \
                  decltype(relevantSourceSamples) : \
                  omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) \
        initializer(omp_priv = omp_orig)

      // The random number generator itself is stateless (has no members which
      // are modified). Hence, it may be shared by threads.
      // auto rng = std::make_unique<rti::rng::cstdlib_rng>();
      auto rng = std::make_unique<rti::rng::mt64_rng>();

      std::cout << "Starting stop watch" << std::endl;
      auto timer = rti::util::timer {};

      /*** Sample to find relevant area of the source ***/
      #pragma omp parallel \
        reduction(+ : geohitc, nongeohitc)                  \
        reduction(hit_accumulator_combine : hitAccumulator) \
        reduction(merge_vectors : relevantSourceSamples)
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
        auto rtiContext = mFactory.get_new_context(geometryID, rtigeometry, reflectionModel,
                                                   hitAccumulator, boundaryID, mBoundary,
                                                   boundaryReflection, *rng, *rngSeed2,
                                                   *particle);

        // Initialize (also takes care for the initialization of the Embree context)
        rtiContext->init();

        // Estimate the distribution
        rtiContext->set_initial_ray_weight(1); // standard weight for a ray
        size_t raycnt = 0;
        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) numprerays; ++idx) {

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
          auto rayOriginCopy = rti::util::triple<float>
            {rayhit.ray.org_x, rayhit.ray.org_y, rayhit.ray.org_z};

          //RAYSRCLOG(rayhit);

          // std::cout
          //   << "new ray == ("
          //   << rayhit.ray.org_x << " " << rayhit.ray.org_y << " " << rayhit.ray.org_z
          //   << ") ("
          //   << rayhit.ray.dir_x << " " << rayhit.ray.dir_y << " " << rayhit.ray.dir_z
          //   << ")" << std::endl << std::flush;


          //if_RLOG_PROGRESS_is_set_print_progress(raycnt, numprerays);

          auto reflect = false;
          auto sumOfRelevantValuesDroped = 0.0;
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


            if (rayhit.hit.geomID == geometryID) {
              auto hitprimID = (dynamic_cast<rti::trace::point_cloud_context<numeric_type>*>
                                (rtiContext.get())->mGeoHitPrimIDs)[0];
              if (relativeerrors[hitprimID] > rethreshold) {
              // if (rtigeometry.get_relevance(hitprimID)) {
                sumOfRelevantValuesDroped +=
                  rtiContext->get_value_of_last_intersect_call();
              }
            }

            // { // Debug
            //   //RAYLOG(rayhit, rtiContext->tfar);
            //   auto p1 = rti::util::triple<float> {rayhit.ray.org_x, rayhit.ray.org_y, rayhit.ray.org_z};
            //   auto tfar_ = (float) (rtiContext->tfar > 10 ? 10.0f : rtiContext->tfar);
            //   auto p2 = rti::util::triple<float> {p1[0] + tfar_ * rayhit.ray.dir_x,
            //                                       p1[1] + tfar_ * rayhit.ray.dir_y,
            //                                       p1[2] + tfar_ * rayhit.ray.dir_z};
            //   rayLogBuff.push_back({p1, p2});
            // }

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
          } while (reflect);
          assert(sumOfRelevantValuesDroped >= 0 && "Correctness Assumption");
          if (sumOfRelevantValuesDroped > 0) {
            relevantSourceSamples.push_back({rayOriginCopy, sumOfRelevantValuesDroped});
            // { // Debug
            //   RTCRayHit rayhit2;
            //   rayhit2.ray.org_x = pointLogBuff[0];
            //   rayhit2.ray.org_y = pointLogBuff[1];
            //   rayhit2.ray.org_z = pointLogBuff[2];
            //   //RAYSRCLOG(rayhit2);
            //   for (auto& seg : rayLogBuff) {
            //     rti::util::sRayLogVec.push_back(seg);
            //   }
            //   rayLogBuff.clear();
            // }
            // std::cerr << "relevantSourceSamples.size() == " <<  relevantSourceSamples.size() << std::endl << std::flush;

            // if (relevantSourceSamples.size() >= (maxrelevantsourcesamplessizeperthread)) {
            //   break;
            // }
          }
        }
      }
      assert(omp_get_num_threads() == 1 && "Assertion to support code readability");
      std::cerr << "relevantSourceSamples.size() == " <<  relevantSourceSamples.size() << std::endl << std::flush;
      // Gaussian Mixture Model
      //auto gmm = rti::util::pair<std::vector<double> > {};
      auto gmm = compute_GMM_using_R(relevantSourceSamples);

      // Prepare distributions from GMM
      auto numMixtureComponents = gmm.size();
      auto mixtureprobabilities = std::vector<numeric_type> (numMixtureComponents);
      auto gaussdists = std::vector<rti::util::pair<std::normal_distribution<numeric_type> > > (numMixtureComponents);
      auto gaussdistsArma = std::vector<rti::util::pair<arma::mat> > (numMixtureComponents);
      for (size_t idx = 0; idx < gmm.size(); ++idx) {
        auto probability = gmm[idx].get0();
        auto mean = gmm[idx].get1();
        auto variance = gmm[idx].get2();
        mixtureprobabilities[idx] = probability;
        /*****************************************************************************/
        /*  Note: std::normal_distribution takes a standard deviation! Not variance! */
        /*****************************************************************************/
        gaussdists[idx] =
          {std::normal_distribution<numeric_type> {(numeric_type) mean[0], (numeric_type) std::sqrt(variance[0])},
           std::normal_distribution<numeric_type> {(numeric_type) mean[1], (numeric_type) std::sqrt(variance[1])}};
        gaussdistsArma[idx][0] = arma::mat {2, 1}; // 2x1 matrix
        gaussdistsArma[idx][0] << mean[0] << arma::endr << mean[1] << arma::endr;
        gaussdistsArma[idx][1] = arma::mat {2, 2};
        gaussdistsArma[idx][1] << variance[0] << 0.0 << arma::endr << 0.0 << variance[1] << arma::endr;
      }
      auto mixdist = std::discrete_distribution<int>
        (mixtureprobabilities.begin(), mixtureprobabilities.end()); // int is the type of the result
      std::cerr << "Finished computation of sampling distribution" << std::endl << std::flush;

      /*** Do the adaptive (primary) sampling ***/

      // Using a new hit accumulator
      hitAccumulator = rti::trace::hit_accumulator<numeric_type> {rtigeometry.get_num_primitives()};

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
        auto rngSeed3 = std::make_unique<rti::rng::mt64_rng::state>(seed+4);
        auto rng3 = std::make_unique<std::mt19937_64>(); // default seed

        // A dummy counter for the boundary
        auto boundaryCntr = rti::trace::dummy_counter {};

        // thread-local particle
        auto particle = particlefactory.create();

        // We will attach our data to the memory immediately following the context as described
        // in https://www.embree.org/api.html#rtcinitintersectcontext .
        // Note: the memory layout only works with plain old data, that is, C-style structs.
        // Otherwise the compiler might change the memory layout, e.g., with the vtable.
        auto rtiContext = mFactory.get_new_context(geometryID, rtigeometry, reflectionModel,
                                                   hitAccumulator, boundaryID, mBoundary,
                                                   boundaryReflection, *rng, *rngSeed2,
                                                   *particle);

        // Initialize (also takes care for the initialization of the Embree context)
        rtiContext->init();

        // The Cross Entropy Method for Optimization, Botev, Krose, Rubinstein, L'Ecuyer
        // In Handbook of statistics (Vol. 31, pp. 35-59). Elsevier.
        // page 17
        // "The generation of a random vector X=(X1, . . . , Xn) ∈ R^n in Step 2 of
        // Algorithm 2.2 is most easily performed by drawing then coordinates independently
        // from some 2-parameter distribution. In most applications a normal (Gaussian)
        // distribution is employed for each component. Thus, the sampling density f(·;v)
        // of X is characterized by a vector of means μ and a vector of variances σ^2 (and
        // we may write v=(μ,σ2)). The choice of the normal distribution is motivated by
        // the availability of fast normal random number generators on modern statistical
        // software and the fact that the maximum likelihood maximization (or cross-entropy
        // minimization) in (9) yields a very simple solution — at each iteration of the CE
        // algorithm the parameter vectors μ and σ^2 are the vectors of sample means and
        // sample variance of the elements of the set of N^e best performing vectors
        // (that is, the elite set); see, for example, Kroese et al. (2006). In summary,
        // the CE method for continuous optimization with a Gaussian sampling density is as
        // follows."

        // auto coordmean = compute_mean_coord(relevantSourceSamples);
        // auto coordvariance = compute_variance_of_coords(relevantSourceSamples, {coordmean[0], coordmean[1], 0.0});

        // The Cross Entropy Method for Optimization, Botev, Kroese, Rubinstein, L'Ecuyer
        // In Handbook of statistics (Vol. 31, pp. 35-59). Elsevier.
        // page 18
        // "For constrained continuous optimization problems, where the samples are
        // restricted to a subset X⊂R^n, it is often possible to replace the normal
        // sampling with sampling from a truncated normal distribution while retaining
        // the updating formulas (19)–(20). An alternative is to use a beta
        // distribution. Smoothing, as in Step 4, is often crucial to prevent premature
        // shrinking of the sampling distribution. Another approach is toinjectextra
        // variance into thesampling distribution, for example by increasing the components
        // ofσ2, oncethe distribution has degenerated; see the examples below and Botev and
        // Kroese(2004)"

        // Also the following paper contains very relevant information for the application of the CE method.
        // The Cross-Entropy Method for Continiuous Multi-extremal Optimization, Kroese, Porotsky, Rubinstein

        // // vector draw
        // template<typename T>
        // statslib_inline
        // T rmvnorm(const T& mu_par, const T& Sigma_par, const bool pre_chol = false);
        //auto onesample = stats::rmvnorm<std::array<double, 2> >(coordmean, coordvariance);

        auto csrcp = dynamic_cast<rti::ray::source<numeric_type>*>(&mSource);
        auto recorgp = dynamic_cast<rti::ray::rectangle_origin_z<numeric_type>*>(&csrcp->get_origin());
        auto zval = recorgp->get_z_val();
        auto orgdel = recorgp->get_x_y_delimiters();
        auto c1 = orgdel[0];
        auto c2 = orgdel[1];
        assert(c1[0] <= c2[0] && c1[1] <= c2[1] && "Assumption");

        // { // Debug
        //   int tid = omp_get_thread_num();
        //   #pragma omp barrier
        //   if (tid == 0) {
        //     std::cerr << "relevantSourceSamples.size() == " << relevantSourceSamples.size() << std::endl;
        //     std::cerr << "means == {" << means(0) << ", " << means(1) << "}" << std::endl;
        //     std::cerr << "variances == {" << covmat(0,0) << ", " << covmat(1,1) << "}" << std::endl;
        //     std::cerr << "zval == " << zval
        //               << " c1 == {" << c1[0] << ", " << c1[1] << "}"
        //               << " c2 == {" << c2[0] << ", " << c2[1] << "}"
        //               << std::endl << std::flush;
        //   }
        //   #pragma omp barrier
        // }

        std::cerr << "Starting " << numrays << " rays" << std::endl;

        auto raycnt = (size_t) 0;
        auto rejectedsamples = (size_t) 0;
        auto debugcnt = 0u;
        #pragma omp for
        for (size_t idx = 0; idx < (unsigned long long int) numrays; ++idx) {

          // Note: Embrees backface culling does not solve our problem of intersections
          // when starting a new ray very close to or a tiny bit below the surface.
          // For that reason we set tnear to some value.
          // There is a risk though: when setting tnear to some strictly positive value
          // we depend on the length of the direction vector of the ray (rayhit.ray.dir_X).

          particle->init_new();

          RLOG_DEBUG << "NEW: Preparing new ray from source" << std::endl;
          // TODO: FIX: tnear set to a constant
          mSource.fill_ray(rayhit.ray, *rng, *rngSeed1); // fills also tnear!

          { // Adaptive sampling
            //auto rng3 = dynamic_cast<rti::rng::mt64_rng::state*>(rngSeed3.get())->get_mt19937_64_ptr();
            // TODO: There does not seem to be a way to use your own RNG for
            // MVN distributions.
            do {
              // auto sample = stats::rmvnorm (means, covmat);
              // // uses arma::randn()
              // // To change the RNG seed, use arma_rng::set_seed(value)

              // // std::cout << "sample.n_elem == " << sample.n_elem << std::endl;
              // // std::cout << "sample.n_cols == " << sample.n_cols << std::endl;
              // // std::cout << "sample.n_rows == " << sample.n_rows << std::endl;
              // auto sx_ = sample(0);
              // auto sy_ = sample(1);

              auto mixidx = mixdist(*rng3);
              assert(0 <= mixidx && mixidx <= numMixtureComponents && "Correctness Assertion");

              // auto xdist = std::normal_distribution<numeric_type> {(numeric_type) meansvec[0], (numeric_type) variancesvec[0]};
              // auto ydist = std::normal_distribution<numeric_type> {(numeric_type) meansvec[1], (numeric_type) variancesvec[1]};
              auto xdist = gaussdists[mixidx][0];
              auto ydist = gaussdists[mixidx][1];
              auto sx = xdist(*rng3);
              auto sy = ydist(*rng3);

              auto sample = arma::mat(2,1); // number of rows and number of columns
              sample(0,0) = sx;
              sample(1,0) = sy;

              { // Debug
                int tid = omp_get_thread_num();
                if (tid == 0) {
                  if (debugcnt < 1024) {
                    debugcnt += 1;
                    std::cerr << "debugcnt == " << debugcnt << " sample == {" << sx << ", " << sy << "}";
                    if ( !(c1[0] <= sx && sx <= c2[0] && c1[1] <= sy && sy <= c2[1])) {
                      std::cerr << " rejected";
                    }
                    std::cerr << std::endl << std::flush;
                  }
                }
              }
              if ( ! (c1[0] <= sx && sx <= c2[0] &&
                      c1[1] <= sy && sy <= c2[1])) {
                rejectedsamples += 1;
                continue;
              }
              // else; sample accepted
              rayhit.ray.org_x = sx;
              rayhit.ray.org_y = sy;
              rayhit.ray.org_z = zval;
              //std::cerr << "just before break" << std::endl;

              // Compute weight for unbiasing
              auto probgmm = 0.0;
              for(size_t idx = 0; idx < numMixtureComponents; ++idx) {
                auto pnorm = rti::util::normal_pdf_2d
                  (rti::util::pair<double> {sample(0), sample(1)},
                   rti::util::pair<double> {gaussdistsArma[idx][0](0), gaussdistsArma[idx][0](1)},
                   rti::util::pair<double> {gaussdistsArma[idx][1](0,0), gaussdistsArma[idx][1](1,1)});
                probgmm += mixtureprobabilities[mixidx] * pnorm;
              }
              //auto dmvUniv = 1.0 / (((double) c2[0] - c1[0]) * ((double) c2[1] - c1[1]));
              auto rayweight = 1.0 / probgmm;
              rtiContext->set_initial_ray_weight(rayweight);
              rtiContext->init_ray_weight();

              { // Debug
                int tid = omp_get_thread_num();
                if (tid == 0) {
                  if (debugcnt < 1024) {
                    std::cerr
                      << "probgmm == " << probgmm << " "
                      << "rayweight == " << rtiContext->get_initial_ray_weight() << std::endl;
                  }
                }
              }
              break;
            } while (true);
          }
          //std::cerr << "HERE" << std::endl << std::flush;

          RAYSRCLOG(rayhit);

          if_RLOG_PROGRESS_is_set_print_progress(raycnt, numrays);

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
        RLOG_TRACE << "before area calc" << std::endl << std::flush;
        if (rtiContext->compute_exposed_areas_by_sampling()) {
          std::cerr
            << "###############"
            << "### WARNING ### Computing exposed area by sampling does not work"
            << "###############" << std::endl;
          assert(false && "Correctness Assertion");
          // // Embree Documentation: "Passing NULL as function pointer disables the registered callback function."
          // rtcSetGeometryIntersectFilterFunction(rtcgeometry, nullptr);
          // rtcSetGeometryIntersectFilterFunction(rtcboundary, nullptr);
          // rtcCommitGeometry(rtcgeometry);
          // rtcCommitGeometry(rtcboundary);
          // compute_exposed_areas_by_sampling(geo, scene, hitAccumulator, geometryID, *rng, *rngSeed1);
        } else {
          use_entire_areas_of_primitives_as_exposed(rtigeometry, hitAccumulator);
        }
        RLOG_TRACE << "after area calc" << std::endl << std::flush;
      }
      // Assertion: hitAccumulator is reduced to one instance by openmp reduction

      result.timeNanoseconds = timer.elapsed_nanoseconds();
      result.hitAccumulator = std::make_unique<rti::trace::hit_accumulator<numeric_type> >(hitAccumulator);
      result.hitc = geohitc;
      result.nonhitc = nongeohitc;

      // auto raylog = RAYLOG_GET_PTR();
      // if (raylog != nullptr) {
      //   auto raylogfilename = "raylog.vtp";
      //   std::cout << "Writing ray log to " << raylogfilename << std::endl;
      //   rti::io::vtp_writer<float>::write(raylog, raylogfilename);
      // }
      // auto raysrclog = RAYSRCLOG_GET_PTR();
      // if (raysrclog != nullptr) {
      //   auto raysrclogfilename = "raysrclog.vtp";
      //   std::cout << "Writing ray src log to " << raysrclogfilename << std::endl;
      //   rti::io::vtp_writer<float>::write(raysrclog, raysrclogfilename);
      // }

      return result;
    }

  private:
    rti::geo::i_factory<numeric_type>& mFactory;
    rti::geo::i_boundary<numeric_type>& mBoundary;
    rti::ray::i_source& mSource;
    rti::particle::i_particle_factory<numeric_type>& particlefactory;

    rti::geo::i_geometry<numeric_type>& rtigeometry;
    RTCGeometry& rtcgeometry;
    RTCGeometry& rtcboundary;
    unsigned int geometryID = std::numeric_limits<unsigned int>::max(); // some initial value
    unsigned int boundaryID = std::numeric_limits<unsigned int>::max(); // some initial value
    RTCDevice& device;
    RTCScene scene; // no pointer or reference. The instance is saved here.

  };
}}
