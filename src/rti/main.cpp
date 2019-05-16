#include <algorithm>
#include <iostream>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

#include "rti/disc_geometry_from_gmsh.hpp"
#include "rti/dummy_ray_source.hpp"
#include "rti/logger.hpp"
#include "rti/sphere_geometry_from_gmsh.hpp"
#include "rti/test_pool.hpp"
#include "rti/test_run.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp"

namespace rti {
  namespace main {

    void init() {
      // Initialize global logger
      rti::logger::init();
      // Enable Flush-to-Zero and Denormals-are-Zero for the MXSCR status and
      // control registers for performance reasons.
      // See: https://embree.github.io/api.html#mxcsr-control-and-status-register
      _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
      _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
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
        << "and 2 indicates Parallel Patterns Library (PPL)";
    }
  }
}

int main(int argc, char* argv[]) {
  rti::main::init();
  rti::gmsh_reader& gmshReader = rti::gmsh_reader::getInstance(argc, argv);

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  rti::main::print_rtc_device_info(device);

  rti::test_pool rtiTPool;
  //i_geometry_from_gmsh* geo = std::make_unique<triangle_geometry_from_gmsh>(device, gmshReader);
  //i_ray_source* raySource = std::make_unique<dummy_ray_source>();
  rti::dummy_ray_source raySource;
  rti::triangle_geometry_from_gmsh triangleGeo(device, gmshReader);
  rti::disc_geometry_from_gmsh discGeo(device, gmshReader);
  rti::sphere_geometry_from_gmsh sphereGeo(device, gmshReader);
  rti::test_run testRunTriangle(triangleGeo, raySource);
  rti::test_run testRunDisc(discGeo, raySource);
  rti::test_run testRunSphere(sphereGeo, raySource);
  rtiTPool.add_test_run(testRunTriangle);
  rtiTPool.add_test_run(testRunDisc);
  rtiTPool.add_test_run(testRunSphere);
  //rtiTPool.add_test_run(testRunDisc);
  // rtiTPool.add_test_run(std::make_unique<disc_geometry_from_gmsh>(device));
  // rtiTPool.add_test_run(std::make_unique<sphere_geometry_from_gmsh>(device));
  auto results = rtiTPool.run();
  for (auto& result : results) {
    std::cout << result.to_string() << std::endl;
  }

  // gmsh::fltk::run();
}

