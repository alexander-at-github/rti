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

#include "rti/command_line_options.hpp"
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
#include "rti/sphere_geometry_from_gmsh.hpp"
#include "rti/test_result.hpp"
#include "rti/tracer.hpp"
#include "rti/trace_result.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp"
#include "rti/vtp_point_cloud_reader.hpp"

namespace rti {
  namespace main_rt {

    void init(int argc, char* argv[]) {
      rti::command_line_options::init(argc, argv);
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

      std::string maxThreadsStr = rti::command_line_options::get_instance().
        get_option_value(rti::command_line_options::option_type::MAX_THREADS);
      if ( ! maxThreadsStr.empty() ) {
        // unsigned long maxThreadsTmp = std::stoul(maxThreadsStr);
        // //assert(0 <= maxThreadsTmp); // Always true
        // assert(maxThreadsTmp <= std::numeric_limits<int>::max() &&
        //        "Number of threads violate assumption");
        // int maxThreads = (int) maxThreadsTmp;
        int maxThreads = std::stoi(maxThreadsStr);

        if (maxThreads < omp_get_max_threads()) {
          omp_set_num_threads(maxThreads);
        }
      }
      std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl;
    }

    void print_rtc_device_info(RTCDevice pDevice) {
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED) << std::endl;
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED) << std::endl;
      RLOG_INFO
        << "RTC_DEVICE_PROPERTY_TASKING_SYSTEM == "
        << rtcGetDeviceProperty(pDevice,RTC_DEVICE_PROPERTY_TASKING_SYSTEM) << std::endl
        << "0 indicates internal tasking system" << std::endl
        << "1 indicates Intel Threading Building Blocks (TBB)" << std::endl
        << "2 indicates Parallel Patterns Library (PPL)" << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {

  rti::main_rt::init(argc, argv);
  //rti::gmsh_reader& gmshReader = rti::gmsh_reader::getInstance();
  auto infilename = rti::command_line_options::get_instance().
    get_option_value(rti::command_line_options::option_type::INFILE_NAME);

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  rti::main_rt::print_rtc_device_info(device);

  rti::disc_origin_x<float> origin(0, 0, 0, 0.5);
  rti::cosine_direction<float> direction(
    {rti::triple<float> {1.f, 0.f, 0.f},
     rti::triple<float> {0.f, 1.f, 0.f},
     rti::triple<float> {0.f, 0.f, 1.f}});

  rti::ray_source<float> source(origin, direction);

  // rti::ray_source<float> source(
  //   std::make_unique<rti::disc_origin_x<float> >(0, 0, 0, 0.5),
  //   std::make_unique<rti::cosine_direction<float> >(
  //     rti::triple<rti::triple<float> > {
  //       {1.f, 0.f, 0.f},
  //       {0.f, 1.f, 0.f},
  //       {0.f, 0.f, 1.f}}));

  //rti::triangle_geometry_from_gmsh geometry(device, gmshReader);
  //rti::oriented_disc_geometry_from_gmsh orntdDiscGeo(device, gmshReader);
  //rti::disc_geometry_from_gmsh discGeo(device, gmshReader);
  //rti::sphere_geometry_from_gmsh geometry(device, gmshReader);

  rti::vtp_point_cloud_reader<float> pntCldReader (infilename);
  rti::point_cloud_geometry<float> geometry (device, pntCldReader);

  rti::tracer tracer (geometry, source);

  rti::trace_result result = tracer.run();

  std::cout << result << std::endl;

  return EXIT_SUCCESS;
}
