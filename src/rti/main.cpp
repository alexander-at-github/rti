#include <algorithm>
#include <iostream>
#include <string>

#include <assert.h>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

#include "rti/logger.hpp"
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
  }
}

int main(int argc, char* argv[]) {
  rti::main::init();
  
  // Set-up gmsh
  gmsh::initialize();
  // Set Gmsh to print information to the terminal
  gmsh::option::setNumber("General.Terminal", 1);
  std::string inpath;
  if (argc > 1) {
    inpath = argv[1];
  } else {
    // Try default file name
    //inpath = "../resources/test.cyl.msh";
    // Gmsh says 14 nodes and 44 elements for box.msh
    // Gmsh says 314 nodes, 704 elements, 8 points, 72 lines, and 624 triangles for box.fine.msh
    inpath = "/home/alexander/cdr/rti/resources/box.fine.msh";
  }
  std::cout << "Reading input file " << inpath << std::endl;
  // Read gmsh fille from arguments
  // Equivalent to the File->Open menu in the gmsh gui application.
  // This file should be a mesh.
  gmsh::open(inpath);

  // TODO: Verify the input.

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  RTCScene scene = rtcNewScene(device);

  // No scene flags for now.
  rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
  // Selecting higher build quality results in better rendering performance but slower
  // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
  RTCBuildQuality bbquality = RTC_BUILD_QUALITY_HIGH;
  rtcSetSceneBuildQuality(scene, bbquality);

  std::unique_ptr<rti::i_geometry_from_gmsh> geometry_from_gmsh = nullptr;
  if (true) {
    geometry_from_gmsh = std::make_unique<rti::triangle_geometry_from_gmsh>(device);
  } else if (true) { // TODO: FIX
    //geometry_from_gmsh = std::make_unique<disc_geometry_from_gmsh>(device);
  } else if (true) { // TODO: FIX
    //geometry_from_gmsh = std::make_unique<sphere_geometry_from_gmsh>(device);
  }
  //geometry_from_gmsh->read_from_gmsh();
  // Invert surface normals.
  geometry_from_gmsh->invert_surface_normals();
  RTCGeometry geometry = geometry_from_gmsh->get_rtc_geometry();
    
  rtcCommitGeometry(geometry);

  auto geomID = rtcAttachGeometry(scene, geometry);
  rtcCommitScene(scene);
  // Release geomtery.
  // The geometry object will be destructed when the scene object is
  // destructed, because the geometry object is attached to the scene
  // objec = 0t
  rtcReleaseGeometry(geometry);

  // Ray queries
  // The Embree tutorial uses one context per ray.
  RTCIntersectContext context;
  rtcInitIntersectContext(&context);

  //RTCRay ray = rti::utils::constructRay(0,0,1, 1,1,-1); // origin (0,0,1) direction (1,1,-1)
  //RTCRay ray = rti::utils::constructRay(0,0,1, 1,0,-1); // origin (0,0,1) direction (1,0,-1)
  //RTCRay ray = rti::utils::constructRay(0,0,0.5, 1,0,0);
  RTCRay ray = rti::utils::constructRay(0,0,0.5, 1,0.1,0); // origin (0,0,1) direction (1,1,-1)
  RTCRayHit rayhit = rti::utils::constructRayHit(ray, rti::utils::constructHit()); // TODO
  BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "Before rtcIntersect1()";
  BOOST_LOG_SEV(rti::mRLogger, blt::trace) <<  "rayhit.ray.tfar=" << rayhit.ray.tfar;
  rtcIntersect1(scene, &context, &rayhit);
  BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "After rtcIntersect1()";
  BOOST_LOG_SEV(rti::mRLogger, blt::trace) <<  "rayhit.ray.tnear=" << rayhit.ray.tnear
                                    << " rayhit.ray.tfar=" << rayhit.ray.tfar
                                    << " rayhit.hit.Ng_x=" << rayhit.hit.Ng_x
                                    << " rayhit.hit.Ng_y=" << rayhit.hit.Ng_y
                                    << " rayhit.hit.Ng_z=" << rayhit.hit.Ng_z
                                    << " rayhit.hit.primID=" << rayhit.hit.primID;
  unsigned int primID = rayhit.hit.primID;
  BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "hit-primID=" << primID
                                           << " " << geometry_from_gmsh->prim_to_string(primID);

  gmsh::fltk::run();
  gmsh::finalize();
}

