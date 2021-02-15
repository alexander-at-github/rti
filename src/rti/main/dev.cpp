#include <algorithm>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
//#include <xmmintrin.h> // SSE
#include <immintrin.h> // AVX
#include <pmmintrin.h> // SSE3
#include <tmmintrin.h> // SSSE3
#include <smmintrin.h> // SSE4
#include <nmmintrin.h> // SSE4
// Can we replace all the includes for the instrinsics with one single include?
//#include <x86intrin.h>

#include "../geo/boundary_x_y.hpp"
#include "../geo/point_cloud_disc_factory.hpp"
#include "../geo/point_cloud_disc_geometry.hpp"
#include "../geo/point_cloud_sphere_geometry.hpp"
#include "../geo/triangle_factory.hpp"
#include "../geo/triangle_geometry.hpp"
#include "../io/vtp_point_cloud_reader.hpp"
#include "../io/christoph/vtu_point_cloud_reader.hpp"
#include "../io/christoph/vtu_triangle_reader.hpp"
#include "../io/vtp_triangle_reader.hpp"
#include "../io/vtp_writer.hpp"
#include "../io/xaver/vtu_point_cloud_reader.hpp"
#include "../particle/i_particle.hpp"
#include "../ray/constant_origin.hpp"
#include "../ray/cosine_direction_z.hpp"
#include "../ray/disc_origin_x.hpp"
#include "../ray/disc_origin_z.hpp"
#include "../ray/nonrng_loc_rng_dir_source.hpp"
#include "../ray/non_mc_source.hpp"
#include "../ray/source.hpp"
#include "../ray/rectangle_origin_z.hpp"
#include "../ray/rectangle_origin_z_v2.hpp" // testing
#include "../reflection/diffuse.hpp"
#include "../trace/point_cloud_context_simplified.hpp"
#include "../trace/tracer.hpp"
#include "../trace/triangle_context_simplified.hpp"
#include "../trace/result.hpp"
#include "../util/clo.hpp"
#include "../util/logger.hpp"
#include "../util/ray_logger.hpp"
#include "../util/utils.hpp"

namespace rti {
  namespace main {

    std::unique_ptr<rti::util::clo::manager> init(int argc, char* argv[]) {
      // Setup command line arguments parser
      auto optMan = std::make_unique<rti::util::clo::manager>();
      optMan->addCmlParam(rti::util::clo::string_option
        {"MAX_THREADS", {"--max-threads", "-m"}, "specifies the maximum number of threads used", false});
      optMan->addCmlParam(rti::util::clo::string_option
        {"INPUT_FILE", {"--infile", "-i"}, "specifies the path of the input file", true});
      // We might want the output file to be mandatory in the future
      optMan->addCmlParam(rti::util::clo::string_option
        {"OUTPUT_FILE", {"--outfile", "-o"}, "specifies the path of the output file", false});
      optMan->addCmlParam(rti::util::clo::string_option
        {"NUM_RAYS", {"--number-of-rays", "--n-rays", "-r"}, "specifies the number of rays to use", false});
      optMan->addCmlParam(rti::util::clo::bool_option
        {"TRIANGLES", {"--triangles"}, "sets triangles as surface primitives"});
      optMan->addCmlParam(rti::util::clo::bool_option
        {"DISCS", {"--discs"}, "sets discs as surface primitives"});
      optMan->addCmlParam(rti::util::clo::string_option
        {"STICKING_COEFFICIENT", {"--sticking-coefficient", "--sticking-c", "--sticking", "-s"},
         "specifies the sticking coefficient of the surface", false});
      optMan->addCmlParam(rti::util::clo::bool_option
        {"SINGLE_HIT", {"--single-hit", "--single"}, "sets single-hit intersections for the ray tracer"});
      bool succ = optMan->parse_args(argc, argv);
      if (!succ) {
        std::cout << optMan->get_usage_msg();
        exit(EXIT_FAILURE);
      }

      std::string maxThreadsStr = optMan->get_string_option_value("MAX_THREADS");
      if ( ! maxThreadsStr.empty() ) {
        int maxThreads = std::stoi(maxThreadsStr);
        if (maxThreads < omp_get_max_threads()) {
          omp_set_num_threads(maxThreads);
        }
      }
      std::cout << "Maximum number of threads used == " << omp_get_max_threads() << std::endl;

      // Enable Flush-to-Zero and Denormals-are-Zero for the MXSCR status and
      // control registers for performance reasons.
      //
      // "If using a different tasking system, make sure each rendering thread has the proper mode set."
      // See: https://embree.github.io/api.html#mxcsr-control-and-status-register
      #pragma omp parallel
      {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
      }

      return optMan;
    }

    void print_rtc_device_info(RTCDevice pDevice) {
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED) << std::endl;
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED) << std::endl;
      // RLOG_INFO
      //   << "RTC_DEVICE_PROPERTY_TASKING_SYSTEM == "
      //   << rtcGetDeviceProperty(pDevice,RTC_DEVICE_PROPERTY_TASKING_SYSTEM) << std::endl
      //   << "0 indicates internal tasking system" << std::endl
      //   << "1 indicates Intel Threading Building Blocks (TBB)" << std::endl
      //   << "2 indicates Parallel Patterns Library (PPL)" << std::endl;
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED) << std::endl;
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_RAY_STREAM_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_RAY_STREAM_SUPPORTED) << std::endl;
    }

    std::string get_git_hash() {
      auto cmd = "git rev-parse HEAD";
      auto result = std::string {};
      auto buffer = std::array<char, 128> {};
      auto pipe = std::unique_ptr<FILE, decltype(&pclose)> (popen(cmd, "r"), pclose);
      if (!pipe)
        return "";
      while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
      return result;
    }
  }
}

int main(int argc, char* argv[]) {

  // We are using floats. There would actually not be any benefit of using
  // double because Embree can work only with floats internally.
  using numeric_type = float;
  using namespace rti;

  auto cmlopts = main::init(argc, argv);
  auto infilename = cmlopts->get_string_option_value("INPUT_FILE");
  auto outfilename = cmlopts->get_string_option_value("OUTPUT_FILE");

  // Enable huge page support.
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  main::print_rtc_device_info(device);

  auto stickingCDefault = 0.8f;
  auto stickingC = stickingCDefault;
  try {
    stickingC = static_cast<float> (std::stod(cmlopts->get_string_option_value("STICKING_COEFFICIENT")));
  } catch (...) {}
  if ( ! (0 < stickingC && stickingC <= 1)) {
    std::cout << "Warning: sticking coefficient has been reset to sane value." << std::endl;
    stickingC = 1.0f; // reset to sane value
  }
  std::cout << "sticking coefficient == " << stickingC << std::endl;
  
  // Handle input file
  if (vtksys::SystemTools::GetFilenameLastExtension(infilename) == ".vtk") {
    std::cout
      << "Recognized \".vtk\" input file file extension. "
      << "Please convert the input file to vtp (e.g., using Paraview) and rerun "
      << argv[0] << " with the new file." << std::endl;
    std::cout << "terminating" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (cmlopts->get_bool_option_value("TRIANGLES")) {
    std::cerr << "Triangles in this version of " << argv[0] << " not supported " << std::endl;
    exit(EXIT_FAILURE);
  }
  auto reader = io::vtp_point_cloud_reader<numeric_type> {infilename};
  auto geometry = geo::point_cloud_disc_geometry<numeric_type> {device, reader};
  
  // Compute bounding box
  auto bdbox = geometry.get_bounding_box();

  auto boundary = geo::boundary_x_y<numeric_type>
    {device, bdbox, geo::bound_condition::PERIODIC, geo::bound_condition::PERIODIC};
  // auto boundary = geo::boundary_x_y<numeric_type>
  //   {device, bdbox, geo::bound_condition::REFLECTIVE, geo::bound_condition::REFLECTIVE};
  // Prepare source
  auto zmax = std::max(bdbox[0][2], bdbox[1][2]);
  auto originC1 = util::pair<numeric_type> {bdbox[0][0], bdbox[0][1]};
  auto originC2 = util::pair<numeric_type> {bdbox[1][0], bdbox[1][1]};

  // ASSUMPTION: the source is on a plain above (positive values) the structure
  // such that z == c for some constant c. (That is in accordance with the silvaco
  // verification instances.)
  //

  //auto origin = ray::rectangle_origin_z<numeric_type> {zmax, originC1, originC2};
  // auto origin = ray::rectangle_origin_z_v2<numeric_type> {zmax, originC1, originC2};

  // // Cosine direction in the opposite direction of the z-axis
  // // auto direction = ray::cosine_direction<numeric_type> {
  // //   {util::triple<numeric_type> {0.f, 0.f, -1.f},
  // //    util::triple<numeric_type> {0.f, 1.f,  0.f},
  // //    util::triple<numeric_type> {1.f, 0.f,  0.f}}};
  // auto direction = ray::cosine_direction_z<numeric_type> {};
  // auto source = ray::source<numeric_type> {origin, direction};

  auto numrays = 128 * 1024ull; // default value // magic number
  auto numraysstr = cmlopts->get_string_option_value("NUM_RAYS");
  try {
    numrays = std::stoull(numraysstr);
  } catch (...) {}

  auto direction = ray::cosine_direction_z<numeric_type> {};
  // auto source = ray::non_mc_source<numeric_type> {zmax, originC1, originC2, numrays, direction};
  auto source = ray::nonrng_loc_rng_dir_source<numeric_type> {zmax, originC1, originC2, numrays, direction};

  //// Define particle
  // In order to pass a value into a local class we need to declare a static
  // variable. This is like capturing variables in a lambda expression.
  static numeric_type stickingStatic = stickingC;
  class particle_t : public particle::i_particle<numeric_type> {
  public:
    numeric_type
    get_sticking_probability
    (RTCRay& rayin,
     RTCHit& hitin,
     geo::meta_geometry<numeric_type>& geometry,
     rng::i_rng& rng,
     rng::i_rng::i_state& rngstate) override final {
      return stickingStatic;
    }

    void init_new() override final {}
  };

  using reflection = reflection::diffuse<numeric_type>;
  auto tracer = trace::tracer<numeric_type, particle_t, reflection>
    {geometry, boundary, source, numrays};
  auto result = tracer.run();
  std::cout << result << std::endl;
  //std::cout << *result.hitAccumulator << std::endl;

  std::cout << "after tracer" << std::endl << std::flush;

  if ( ! outfilename.empty()) {
    // Write output to file
    if (vtksys::SystemTools::GetFilenameLastExtension(outfilename) != ".vtp") {
      std::cout << "Appending .vtp to the given file name" << std::endl;
      outfilename.append(".vtp");
    }
    auto bbpath = vtksys::SystemTools::GetFilenamePath(outfilename);
    if( ! bbpath.empty()) {
      bbpath.append("/");
    }
    auto bbfilename = bbpath + vtksys::SystemTools::
      GetFilenameWithoutExtension(outfilename).append(".bounding-box.vtp");

    std::cout << "Writing output to " << outfilename << std::endl << std::flush;
    auto cmdstr = util::foldl<std::string, std::string>
      ([](auto p1, auto const p2){return (p1+=" ")+=p2;}, "", std::vector<std::string> (argv, argv+argc));
    std::cerr << "cmdstr == " << cmdstr << std::endl;
    auto geoname = typeid(&geometry).name();    
    io::vtp_writer<numeric_type>::write
      (geometry, 
       *result.hitAccumulator,
       outfilename,
       // std::vector<rti::util::pair<std::string> >
       {{"running-time[ns]", std::to_string(result.timeNanoseconds)},
        {"git-hash", main::get_git_hash()},
        {"cmd", cmdstr},
        {"geo-name", geoname}});
    std::cout << "Writing bounding box to " << bbfilename << std::endl;
    io::vtp_writer<numeric_type>::write(boundary, bbfilename);

    auto raylog = RAYLOG_GET_PTR();
    if (raylog != nullptr) {
      auto raylogfilename = vtksys::SystemTools::
        GetFilenameWithoutExtension(outfilename).append(".ray-log.vtp");
      std::cout << "Writing ray log to " << raylogfilename << std::endl;
      io::vtp_writer<numeric_type>::write(raylog, raylogfilename);
    }
    auto raysrclog = RAYSRCLOG_GET_PTR();
    if (raysrclog != nullptr) {
      auto raysrclogfilename = vtksys::SystemTools::
        GetFilenameWithoutExtension(outfilename).append(".ray-src-log.vtp");
      std::cout << "Writing ray src log to " << raysrclogfilename << std::endl;
      io::vtp_writer<numeric_type>::write(raysrclog, raysrclogfilename);
    }
  }

  rtcReleaseDevice(device);

  return EXIT_SUCCESS;
}

