#include <algorithm>
#include <iostream>
#include <string>
#include <typeinfo> // need?

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

//#include "rti/disc_geometry_from_gmsh.hpp"
#include "rti/dummy_ray_source.hpp"
#include "rti/logger.hpp"
//#include "rti/sphere_geometry_from_gmsh.hpp"
#include "rti/test_pool.hpp"
#include "rti/test_run.hpp"
#include "rti/triangle_geometry_from_gmsh.hpp"
#include "rti/utils.hpp"

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
      BOOST_LOG_SEV(rti::mRLogger, blt::debug)
        << "RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED);
      BOOST_LOG_SEV(rti::mRLogger, blt::debug)
        << "RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED == "
        << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_POINT_GEOMETRY_SUPPORTED);
    }
  }
}

int main(int argc, char* argv[]) {
  rti::main::init();
  rti::gmsh_reader& gmshReader = rti::gmsh_reader::getInstance(argc, argv);

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  //rti::main::printI_rtc_device_info(device);

  rti::test_pool rtiTPool;
  //i_geometry_from_gmsh* geo = std::make_unique<triangle_geometry_from_gmsh>(device, gmshReader);
  //i_ray_source* raySource = std::make_unique<dummy_ray_source>();
  rti::triangle_geometry_from_gmsh geo(device, gmshReader);
  rti::dummy_ray_source raySource;
  rti::test_run testRun(geo, raySource);
  rtiTPool.add_test_run(testRun);
  // rtiTPool.add_test_run(std::make_unique<disc_geometry_from_gmsh>(device));
  // rtiTPool.add_test_run(std::make_unique<sphere_geometry_from_gmsh>(device));
  auto results = rtiTPool.run();
  for (auto& result : results) {
    std::cout << result.to_string() << std::endl;
  }

  gmsh::fltk::run();
}

