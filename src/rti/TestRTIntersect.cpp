#include "rti/Tools.h"

#include <embree3/rtcore.h>
#include <iostream>
#include <pmmintrin.h>
#include <string>
#include <xmmintrin.h>

int main(int argc, char* argv[]) {
  // Enable Flush-to-Zero and Denormals-are-Zero for the MXSCR status and
  // control registers for performance reasons.
  // See: https://embree.github.io/api.html#mxcsr-control-and-status-register
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  RTCScene scene = rtcNewScene(device);
  // No scene flags for now.
  //rtcSetSceneFlags(scene, NULL); // RTC_SCENE_FLAG_NONE ?
  std::cout << rtcGetSceneFlags(scene) << std::endl;
  // Set high scene build quality for better rendering (ray tracing?)
  // performance but slower scene commit.
  rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);
  // TODO: which type of geometry to use?
  //
  // "Points with per vertex radii are supported with sphere, ray-oriented
  // discs, and normal-oriented discs geometric represetntations. Such point
  // geometries are created by passing RTC_GEOMETRY_TYPE_SPHERE_POINT,
  // RTC_GEOMETRY_TYPE_DISC_POINT, or RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT to
  // the rtcNewGeometry function."
  // Source: https://embree.github.io/api.html#rtc_geometry_type_point
  RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
  // TODO: set geometry
  rtcCommitGeometry(geometry);
  auto geomID = rtcAttachGeometry(scene, geometry);
  // Release geomtery. Because the geometry object is attached to the scene
  // object the geometry object will be destructed when the scene object is
  // destructed.
  rtcReleaseGeometry(geometry);

  // Ray queries
  RTCIntersectContext context;
  rtcInitIntersectContext(&context);

  RTCRay ray = rti::Tools::getRay(0,0,0,1,1,1); //TODO
  RTCRayHit rayhit = rti::Tools::getRayHit(ray, RTCHit{}); // TODO
  rtcIntersect1(scene, &context, &rayhit);
}

