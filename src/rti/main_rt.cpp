#include <algorithm>
#include <iostream>
#include <omp.h>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <immintrin.h> // Need? Added later
#include <pmmintrin.h>
#include <xmmintrin.h>

#include "rti/boundary.hpp"
#include "rti/clo.hpp"
#include "rti/disc_geometry_from_gmsh.hpp"
#include "rti/constant_origin.hpp"
#include "rti/cosine_direction.hpp"
#include "rti/cosine_direction.hpp"
#include "rti/disc_origin_x.hpp"
#include "rti/dummy_direction.hpp"
#include "rti/logger.hpp"
#include "rti/oriented_disc_geometry_from_gmsh.hpp"
#include "rti/point_cloud_geometry.hpp"
#include "rti/ray_source.hpp"
#include "rti/rectangle_origin_z.hpp"
#include "rti/sphere_geometry_from_gmsh.hpp"
#include "rti/test_result.hpp"
#include "rti/tracer.hpp"
#include "rti/trace_result.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp"
#include "rti/vtp_point_cloud_reader.hpp"

namespace rti {
  namespace main_rt {

    std::unique_ptr<rti::clo::manager> init(int argc, char* argv[]) {
      // Setup command line arguments parser
      auto optMan = std::make_unique<rti::clo::manager>();
      optMan->addCmlParam(rti::clo::string_option
        {"MAX_THREADS", {"--max-threads", "-m"}, "specifies the maximum number of threads used", false});
      optMan->addCmlParam(rti::clo::string_option
        {"INPUT_FILE", {"--infile", "-i"}, "specifies the path of the input file", true});
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
      std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl;

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
    }
  }
}

int main(int argc, char* argv[]) {

  // We are using floats. There would actually not be any benefit of using
  // double because Embree can work only with floats internally.
  using numeric_type = float;

  auto optMan = rti::main_rt::init(argc, argv);
  //rti::gmsh_reader& gmshReader = rti::gmsh_reader::getInstance();
  auto infilename = optMan->get_string_option_value("INPUT_FILE");

  // Enable huge page support.
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  rti::main_rt::print_rtc_device_info(device);

  // rti::ray_source<numeric_type> source(
  //   std::make_unique<rti::disc_origin_x<numeric_type> >(0, 0, 0, 0.5),
  //   std::make_unique<rti::cosine_direction<numeric_type> >(
  //     rti::triple<rti::triple<numeric_type> > {
  //       {1.f, 0.f, 0.f},
  //       {0.f, 1.f, 0.f},
  //       {0.f, 0.f, 1.f}}));

  //rti::triangle_geometry_from_gmsh geometry(device, gmshReader);
  //rti::oriented_disc_geometry_from_gmsh orntdDiscGeo(device, gmshReader);
  //rti::disc_geometry_from_gmsh discGeo(device, gmshReader);
  //rti::sphere_geometry_from_gmsh geometry(device, gmshReader);
  //rti::disc_origin_x<numeric_type> origin(0, 0, 0, 0.5);

  auto pntCldReader = rti::vtp_point_cloud_reader<numeric_type> {infilename};
  auto geometry = rti::point_cloud_geometry<numeric_type> {device, pntCldReader};

  // Compute bounding box
  auto bdBox = geometry.get_bounding_box();
  // for (auto const& bb : bdBox)
  //   for (auto const& cc : bb)
  //     std::cout << cc << " ";
  // std::cout << std::endl;

  // Prepare boundary
  auto boundary = rti::boundary_x_y<numeric_type> {device, bdBox};

  // Prepare source
  auto zmax = std::max(bdBox[0][2], bdBox[1][2]);
  auto originC1 = rti::pair<numeric_type> {bdBox[0][0], bdBox[0][1]};
  auto originC2 = rti::pair<numeric_type> {bdBox[1][0], bdBox[1][1]};

  // ASSUMPTION: the source is on a plain above (positive values) the structure
  // such that z == c for some constant c. (That is in accordance with the silvaco
  // verification instances.)
  auto origin = rti::rectangle_origin_z<numeric_type> {zmax, originC1, originC2};
  // Cosine direction in the opposite direction of the z-axis
  auto direction = rti::cosine_direction<numeric_type> {
    {rti::triple<numeric_type> {0.f, 0.f, -1.f},
     rti::triple<numeric_type> {0.f, 1.f,  0.f},
     rti::triple<numeric_type> {1.f, 0.f,  0.f}}};
  auto source = rti::ray_source<numeric_type> {origin, direction};
  auto tracer = rti::tracer<numeric_type> {geometry, boundary, source};
  auto result = tracer.run();
  //std::cout << result << std::endl;
  return EXIT_SUCCESS;
}
