#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
//#include <xmmintrin.h> // SSE
#include <immintrin.h> // AVX
#include <pmmintrin.h> // SSE3
#include <tmmintrin.h> // SSSE3
#include <smmintrin.h> // SSE4
#include <nmmintrin.h> // SSE4
// Can we replace all the includes for the instrinsics with one single include?
//#include <x86intrin.h>

#include "rti/geo/boundary_x_y.hpp"
#include "rti/geo/disc_geometry_from_gmsh.hpp"
#include "rti/geo/oriented_disc_geometry_from_gmsh.hpp"
#include "rti/geo/point_cloud_disc_factory.hpp"
#include "rti/geo/point_cloud_disc_geometry.hpp"
#include "rti/geo/point_cloud_sphere_geometry.hpp"
#include "rti/geo/sphere_geometry_from_gmsh.hpp"
#include "rti/geo/triangle_factory.hpp"
#include "rti/geo/triangle_geometry_from_gmsh.hpp"
#include "rti/geo/triangle_geometry.hpp"
#include "rti/io/vtp_point_cloud_reader.hpp"
#include "rti/io/christoph/vtu_point_cloud_reader.hpp"
#include "rti/io/christoph/vtu_triangle_reader.hpp"
#include "rti/io/vtp_triangle_reader.hpp"
#include "rti/io/vtp_writer.hpp"
#include "rti/ray/constant_origin.hpp"
#include "rti/ray/cosine_direction.hpp"
#include "rti/ray/disc_origin_x.hpp"
#include "rti/ray/disc_origin_z.hpp"
#include "rti/ray/dummy_direction.hpp"
#include "rti/ray/source.hpp"
#include "rti/ray/rectangle_origin_z.hpp"
#include "rti/test_and_benchmark/test_result.hpp"
#include "rti/trace/tracer.hpp"
#include "rti/trace/result.hpp"
#include "rti/util/clo.hpp"
#include "rti/util/logger.hpp"
#include "rti/util/ray_logger.hpp"
#include "rti/util/utils.hpp"

namespace rti {
  namespace main_rt {

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
        {"STICKING_COEFFICIENT", {"--sticking-coefficient", "--sticking-c", "--sticking", "-s"}, "specifies the sticking coefficient of the surface", false});
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

  auto optMan = rti::main_rt::init(argc, argv);
  auto infilename = optMan->get_string_option_value("INPUT_FILE");
  auto outfilename = optMan->get_string_option_value("OUTPUT_FILE");

  // Enable huge page support.
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  rti::main_rt::print_rtc_device_info(device);
  assert(rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED) != 0 &&
         "Error: backface culling is not enabled");

  auto stickingCDefault = 0.8f; // float // default value
  auto stickingC = stickingCDefault;
  try {
    stickingC = static_cast<float> (std::stod(optMan->get_string_option_value("STICKING_COEFFICIENT")));
  } catch (...) {}
  if ( ! (0 < stickingC && stickingC <= 1)) {
    std::cout << "Warning: sticking coefficient has been reset to sane value." << std::endl;
    stickingC = 1.0f; // reset to sane value
  }
  std::cout << "sticking coefficient == " << stickingC << std::endl;
  // create variable
  auto geoFactory = std::unique_ptr<rti::geo::i_factory<numeric_type> > (nullptr);
  if (optMan->get_bool_option_value("TRIANGLES")) {
    RLOG_DEBUG << "recognized cmd line option TRIANGLES" << std::endl;
    if (vtksys::SystemTools::GetFilenameLastExtension(infilename) == ".vtp") {
      RLOG_DEBUG << "recognized .vtp file" << std::endl;
      auto reader = rti::io::vtp_triangle_reader<numeric_type> {infilename};
      geoFactory = std::make_unique<rti::geo::triangle_factory<numeric_type, rti::trace::triangle_context_simplified<numeric_type> > > (device, reader, stickingC);
    } else if (vtksys::SystemTools::GetFilenameLastExtension(infilename) == ".vtu") {
      RLOG_DEBUG << "recognized .vtu file" << std::endl;
      auto reader = rti::io::christoph::vtu_triangle_reader<numeric_type> {infilename};
      geoFactory = std::make_unique<rti::geo::triangle_factory<numeric_type, rti::trace::triangle_context_simplified<numeric_type> > > (device, reader, stickingC);
    }
  } else { // Default
  //} else if (optMan->get_bool_option_value("DISCS")) {
    RLOG_DEBUG << "using DISCS (default option)" << std::endl;
    if (vtksys::SystemTools::GetFilenameLastExtension(infilename) == ".vtp") {
      RLOG_DEBUG << "recognized .vtp file" << std::endl;
      auto reader = rti::io::vtp_point_cloud_reader<numeric_type> {infilename};
      geoFactory = std::make_unique<rti::geo::point_cloud_disc_factory<numeric_type> > (device, reader, stickingC);
    //} else {
    } else if (vtksys::SystemTools::GetFilenameLastExtension(infilename) == ".vtu") {
      RLOG_DEBUG << "recognized .vtu file" << std::endl;
      auto reader = rti::io::christoph::vtu_point_cloud_reader<numeric_type> {infilename};
      geoFactory = std::make_unique<rti::geo::point_cloud_disc_factory<numeric_type> > (device, reader, stickingC);
    }
  }

  // Compute bounding box
  auto bdBox = geoFactory->get_geometry().get_bounding_box();
  // Increase the size of the bounding box by an epsilon on the z achsis.
  auto epsilon = 0.1; //0.1; // -0.1;
  if (bdBox[0][2] > bdBox[1][2]) {
    bdBox[0][2] += epsilon;
  } else {
    bdBox[1][2] += epsilon;
  }
  // std::cerr << "[main] bdBox: ";
  // for (auto const& bb : bdBox)
  //   for (auto const& cc : bb)
  //     std::cerr << cc << " ";
  // std::cerr << std::endl;

  // Prepare boundary
  auto bdBoxAltered = bdBox;
  // bdBoxAltered[0][0] -= 0.1; // x
  // bdBoxAltered[0][1] -= 0.1; // y
  // bdBoxAltered[1][0] += 0.1; // x
  // bdBoxAltered[1][1] += 0.1; // y
  auto boundary = rti::geo::boundary_x_y<numeric_type> {device, bdBoxAltered};

  // Prepare source
  auto zmax = std::max(bdBox[0][2], bdBox[1][2]);
  auto originC1 = rti::util::pair<numeric_type> {bdBox[0][0], bdBox[0][1]};
  auto originC2 = rti::util::pair<numeric_type> {bdBox[1][0], bdBox[1][1]};

  // ASSUMPTION: the source is on a plain above (positive values) the structure
  // such that z == c for some constant c. (That is in accordance with the silvaco
  // verification instances.)
  //
  auto origin = rti::ray::rectangle_origin_z<numeric_type> {zmax, originC1, originC2};
  // auto origin = rti::ray::disc_origin_z<numeric_type> {(originC1[0] + originC2[0])/2,
  //                                                      (originC1[1] + originC2[1])/2,
  //                                                      zmax,
  //                                                      (originC2[0] - originC1[0])/2};

  // Cosine direction in the opposite direction of the z-axis
  auto direction = rti::ray::cosine_direction<numeric_type> {
    {rti::util::triple<numeric_type> {0.f, 0.f, -1.f},
     rti::util::triple<numeric_type> {0.f, 1.f,  0.f},
     rti::util::triple<numeric_type> {1.f, 0.f,  0.f}}};
  auto source = rti::ray::source<numeric_type> {origin, direction};

  auto numraysstr = optMan->get_string_option_value("NUM_RAYS");
  auto numrays = 128 * 1024ull; // default value // magic number
  try {
    numrays = std::stoull(numraysstr);
  } catch (...) {}


  auto tracer = rti::trace::tracer<numeric_type> {*geoFactory, boundary, source, numrays};
  auto result = tracer.run();
  std::cout << result << std::endl;
  //std::cout << *result.hitAccumulator << std::endl;

  if ( ! outfilename.empty()) {
    // Write output to file
    if (vtksys::SystemTools::GetFilenameLastExtension(outfilename) != ".vtp") {
      std::cout << "Appending .vtp to the given file name" << std::endl;
      outfilename.append(".vtp");
    }
    auto bbpath = vtksys::SystemTools::GetFilenamePath(outfilename);
    auto bbfilename = vtksys::SystemTools::
      GetFilenameWithoutExtension(outfilename).append(".bounding-box.vtp");
    { // boost uses portable path separator
      namespace bfs = boost::filesystem;
      bbfilename = (bfs::path{bbpath} / bfs::path{bbfilename}).string();
    }

    std::cout << "Writing output to " << outfilename << std::endl;
    auto cmdstr = rti::util::foldl<std::string, std::string>
      ([](auto p1, auto const p2){return (p1+=" ")+=p2;}, "", std::vector<std::string> (argv, argv+argc));
    std::cerr << "cmdstr == " << cmdstr << std::endl;
    geoFactory->write_to_file(*result.hitAccumulator, outfilename,
                              {{"running-time[ns]", std::to_string(result.timeNanoseconds)},
                               {"git-hash", rti::main_rt::get_git_hash()},
                               {"cmd", cmdstr},
                               // the following line does not really give you the most derived type. FIX
                               {"geo-factory-name", boost::core::demangle(typeid(geoFactory.get()).name())},
                               {"geo-name", boost::core::demangle(typeid(geoFactory->get_geometry()).name())}});
    std::cout << "Writing bounding box to " << bbfilename << std::endl;
    rti::io::vtp_writer<numeric_type>::write(boundary, bbfilename);

    auto raylog = RAYLOG_GET_PTR();
    if (raylog != nullptr) {
      auto raylogfilename = vtksys::SystemTools::
        GetFilenameWithoutExtension(outfilename).append(".ray-log.vtp");
      std::cout << "Writing ray log to " << raylogfilename << std::endl;
      rti::io::vtp_writer<numeric_type>::write(raylog, raylogfilename);
    }
    auto raysrclog = RAYSRCLOG_GET_PTR();
    if (raysrclog != nullptr) {
      auto raysrclogfilename = vtksys::SystemTools::
        GetFilenameWithoutExtension(outfilename).append(".ray-src-log.vtp");
      std::cout << "Writing ray src log to " << raysrclogfilename << std::endl;
      rti::io::vtp_writer<numeric_type>::write(raysrclog, raysrclogfilename);
    }
  }

  rtcReleaseDevice(device);

  return EXIT_SUCCESS;
}
