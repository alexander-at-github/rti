#include <algorithm>
#include <iostream>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <immintrin.h> // Need? Added later
#include <pmmintrin.h>
#include <tbb/tbb.h>
#include <xmmintrin.h>

#include "rti/command_line_options.hpp"
#include "rti/disc_geometry_from_gmsh.hpp"
#include "rti/constant_origin.hpp"
#include "rti/cosine_direction.hpp"
#include "rti/dummy_direction.hpp"
#include "rti/logger.hpp"
#include "rti/oriented_disc_geometry_from_gmsh.hpp"
#include "rti/ray_source.hpp"
#include "rti/sphere_geometry_from_gmsh.hpp"
#include "rti/test_pool.hpp"
#include "rti/test_result.hpp"
#include "rti/test_run.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp"

namespace rti {
  namespace main_testing_rt {

    void init(int argc, char* argv[]) {
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
      //RLOG_DEBUG << "Using " << maxThreads << " threads" << std::endl;
      std::cout << "Using " << maxThreads << " threads" << std::endl;
      static tbb::task_scheduler_init init(maxThreads);
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
  rti::main_testing_rt::init(argc, argv);
  rti::gmsh_reader& gmshReader = rti::gmsh_reader::getInstance();

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  rti::main_testing_rt::print_rtc_device_info(device);

  rti::ray_source<float> source(
    std::make_unique<rti::constant_origin<float> >(0.1, 0, 0),
    //std::make_unique<rti::cosine_direction<float> >()
    std::make_unique<rti::dummy_direction>()
                                );

  rti::triangle_geometry_from_gmsh triangleGeo(device, gmshReader);
  rti::sphere_geometry_from_gmsh sphereGeo(device, gmshReader);
  rti::oriented_disc_geometry_from_gmsh orntdDiscGeo(device, gmshReader);
  rti::disc_geometry_from_gmsh discGeo(device, gmshReader);

  rti::test_run testRunTriangle(triangleGeo, source);
  rti::test_run testRunSphere(sphereGeo, source);
  rti::test_run testRunOrntdDisc(orntdDiscGeo, source);
  rti::test_run testRunDisc(discGeo, source);

  rti::test_pool poolTrngl;
  rti::test_pool poolSphr;
  rti::test_pool poolOrntdDisc;
  rti::test_pool poolDsc;

  // Number of test repetitions (samples).
  size_t reps = 35;
  //size_t reps = 1;
  for(size_t nn = 0; nn < reps; ++nn) {
    poolTrngl.add_test_run(testRunTriangle);
    poolSphr.add_test_run(testRunSphere);
    poolOrntdDisc.add_test_run(testRunOrntdDisc);
    poolDsc.add_test_run(testRunDisc);
  }

  for (auto& pp : std::vector<rti::test_pool> {poolTrngl, poolSphr, poolOrntdDisc, poolDsc}) {
    auto results = pp.run();
    for (auto& result : results) {
      std::cout << result << std::endl;
    }
  }

  // gmsh::fltk::run();
}
