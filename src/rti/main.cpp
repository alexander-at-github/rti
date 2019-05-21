#include <algorithm>
#include <iostream>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <pmmintrin.h>
#include <tbb/tbb.h>
#include <xmmintrin.h>

#include "rti/command_line_options.hpp"
#include "rti/disc_geometry_from_gmsh.hpp"
#include "rti/dummy_ray_source.hpp"
#include "rti/logger.hpp"
#include "rti/sphere_geometry_from_gmsh.hpp"
#include "rti/test_pool.hpp"
#include "rti/test_result.hpp"
#include "rti/test_run.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp"

namespace rti {
  namespace main {

    void init(int argc, char* argv[]) {
      // Initialize global logger
      rti::logger::init();
      rti::command_line_options::init(argc, argv);
      // Enable Flush-to-Zero and Denormals-are-Zero for the MXSCR status and
      // control registers for performance reasons.
      // See: https://embree.github.io/api.html#mxcsr-control-and-status-register
      _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
      _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

      // tbb
      // do we need that? Not necessary. But one can set the number of threads.
      //static tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
      std::string maxThreadsStr = rti::command_line_options::get_instance().
        get_option_value(rti::command_line_options::option_type::MAX_THREADS);
      unsigned int maxThreads = tbb::task_scheduler_init::default_num_threads();
      if ( ! maxThreadsStr.empty()) {
        maxThreads = std::stoul(maxThreadsStr);
      }
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Using " << maxThreads << " threads";
      static tbb::task_scheduler_init init(maxThreads);
      BOOST_LOG_SEV(rti::mRLogger, blt::warning)
        << "tbb::task_scheduler_init::default_num_threads() == "
        << tbb::task_scheduler_init::default_num_threads();
    }

    void print_rtc_device_info(RTCDevice pDevice) {
      BOOST_LOG_SEV(rti::mRLogger, blt::info)
        << "RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED);
      BOOST_LOG_SEV(rti::mRLogger, blt::info)
        << "RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED);
      BOOST_LOG_SEV(rti::mRLogger, blt::info)
        << "RTC_DEVICE_PROPERTY_TASKING_SYSTEM == "
        << rtcGetDeviceProperty(pDevice,RTC_DEVICE_PROPERTY_TASKING_SYSTEM) << std::endl
        << "0 indicates internal tasking system" << std::endl
        << "1 indicates Intel Threading Building Blocks (TBB)" << std::endl
        << "2 indicates Parallel Patterns Library (PPL)";
    }
  }
}

int main(int argc, char* argv[]) {
  rti::main::init(argc, argv);
  rti::gmsh_reader& gmshReader = rti::gmsh_reader::getInstance();

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  rti::main::print_rtc_device_info(device);

  //i_geometry_from_gmsh* geo = std::make_unique<triangle_geometry_from_gmsh>(device, gmshReader);
  //i_ray_source* raySource = std::make_unique<dummy_ray_source>();
  rti::dummy_ray_source raySource;

  rti::triangle_geometry_from_gmsh triangleGeo(device, gmshReader);
  rti::disc_geometry_from_gmsh discGeo(device, gmshReader);
  rti::sphere_geometry_from_gmsh sphereGeo(device, gmshReader);

  rti::test_run testRunTriangle(triangleGeo, raySource);
  rti::test_run testRunDisc(discGeo, raySource);
  rti::test_run testRunSphere(sphereGeo, raySource);

  rti::test_pool poolTrngl;
  rti::test_pool poolDsc;
  rti::test_pool poolSphr;

  // Number of test repetitions (samples).
  size_t reps = 1;
  for(size_t nn = 0; nn < reps; ++nn) {
    poolTrngl.add_test_run(testRunTriangle);
    poolDsc.add_test_run(testRunDisc);
    poolSphr.add_test_run(testRunSphere);
  }

  for (auto & pp : std::vector<rti::test_pool> {poolTrngl, poolDsc, poolSphr}) {
    auto results = pp.run();
    for (auto & result : results) {
      std::cout << result.to_string() << std::endl;
    }
  }

  // rtiTPool.add_test_run(std::make_unique<disc_geometry_from_gmsh>(device));
  // rtiTPool.add_test_run(std::make_unique<sphere_geometry_from_gmsh>(device));
  //auto results = rtiTPool.run();
  //for (auto& result : results) {
  //  std::cout << result.to_string() << std::endl;
  //}

  // gmsh::fltk::run();
}

