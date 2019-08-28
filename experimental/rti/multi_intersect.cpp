#include <cstdlib>
#include <iostream>
#include <limits>

#include <embree3/rtcore.h>

void occluded_filter(const struct RTCFilterFunctionNArguments* args) {
  std::cerr << "occluded function" << std::endl;
  std::cerr << "args->N == " << args->N << std::endl;
  auto primID = RTCHitN_primID(args->hit, args->N, 0);
  std::cerr << "primID == " << primID << std::endl;
  *(args->valid) = 0;
}

void intersect_filter(const struct RTCFilterFunctionNArguments* args) {
  std::cerr << "intersect function" << std::endl;
  std::cerr << "args->N == " << args->N << std::endl;
  auto primID = RTCHitN_primID(args->hit, args->N, 0);
  std::cerr << "primID == " << primID << std::endl;

  auto& tfar = RTCRayN_tfar(args->ray, args->N, 0);
  std::cerr << "tfar == " << tfar << std::endl;
  //tfar = std::numeric_limits<float>::max();
  auto valid = args->valid;
  valid[0] = 0;
}

void print_info(RTCDevice pDevice) {
  std::cerr
    << "RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED == "
    << rtcGetDeviceProperty(pDevice, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED) << std::endl;
}

int main(int argc, char* argv[]) {
  auto numpoints = 3u; // unsigned int
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  print_info(device);
  auto geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);
  struct point_4f_t {
    float xx, yy, zz, radius;
  };
  struct normal_vec_3f_t {
    float xx, yy, zz;
  };
  auto discs = (point_4f_t*)
    rtcSetNewGeometryBuffer(geometry,
                            RTC_BUFFER_TYPE_VERTEX,
                            0, // slot
                            RTC_FORMAT_FLOAT4,
                            sizeof(point_4f_t),
                            numpoints);
  auto normals = (normal_vec_3f_t*)
    rtcSetNewGeometryBuffer(geometry,
                            RTC_BUFFER_TYPE_NORMAL,
                            0, // slot
                            RTC_FORMAT_FLOAT3,
                            sizeof(normal_vec_3f_t),
                            numpoints);
  for (size_t idx = 0; idx < numpoints; ++idx) {
    discs[idx].xx = 0; discs[idx].yy = 0; discs[idx].zz = 100 - idx*10;
    discs[idx].radius = 1;
    normals[idx].xx = 0; normals[idx].yy = 0; normals[idx].zz = -1;
  }

  //rtcSetGeometryOccludedFilterFunction(geometry, &occluded_filter);
  //rtcSetGeometryIntersectFilterFunction(geometry, &intersect_filter);

  rtcCommitGeometry(geometry);

  auto scene = rtcNewScene(device);
  auto bbquality = RTC_BUILD_QUALITY_HIGH;
  rtcSetSceneBuildQuality(scene, bbquality);
  rtcSetSceneFlags(scene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION); // does not seem to be necessary
  auto geometryID = rtcAttachGeometry(scene, geometry);
  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);

  context.filter = intersect_filter; // register filter
  // Maybe it is better to use the filter on the geometry because in our application
  // the filter is not needed for the boundary conditions.

  // "Any API call that sets a property of the scene or geometries contained in
  // the scene count as scene modification, e.g. including setting of intersection
  // filter functions."
  rtcCommitScene(scene);

  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0;
  rayhit.ray.org_y = 0;
  rayhit.ray.org_z = 0;
  rayhit.ray.dir_x = 0;
  rayhit.ray.dir_y = 0;
  rayhit.ray.dir_z = 1;
  //rayhit.ray.tfar = std::numeric_limits<float>::max();
  rayhit.ray.tfar = 95;
  rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

  rtcOccluded1(scene, &context, &rayhit.ray);
  //rtcIntersect1(scene, &context, &rayhit);
  std::cout << "rayhit.ray.tfar = " << rayhit.ray.tfar << std::endl;

  rtcReleaseGeometry(geometry);
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);

  return EXIT_SUCCESS;
}
