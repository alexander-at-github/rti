#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <embree3/rtcore.h>
#include <immintrin.h> // AVX
#include <pmmintrin.h>
#include <xmmintrin.h>

int main(int argc, char* argv[])
{
  {  // if multi-threaded, then this block needs to be parallelized.
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  }

  auto points = std::vector<std::array<float, 3> > {{0,0,0}, {1,1,0}};
  auto normals = std::vector<std::array<float, 3> > {{0,0,1}, {0,0,1}};
  auto spacing = std::vector<float> {std::sqrt(2), std::sqrt(2)};
  
  auto numpoints = points.size();
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  if (rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) == 0) {
    std::cerr << "Error" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  auto geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

  struct point_4f_t {
    float xx, yy, zz, radius;
  };
  struct normal_vec_3f_t {
    float xx, yy, zz;
  };
  auto discbuff = (point_4f_t*)
    rtcSetNewGeometryBuffer(geometry,
                            RTC_BUFFER_TYPE_VERTEX,
                            0, // slot
                            RTC_FORMAT_FLOAT4,
                            sizeof(point_4f_t),
                            numpoints);
  auto normalsbuff = (normal_vec_3f_t*)
    rtcSetNewGeometryBuffer(geometry,
                            RTC_BUFFER_TYPE_NORMAL,
                            0, // slot
                            RTC_FORMAT_FLOAT3,
                            sizeof(normal_vec_3f_t),
                            numpoints);
  for (size_t idx = 0; idx < numpoints; ++idx) {
    discbuff[idx].xx = points[idx][0];
    discbuff[idx].yy = points[idx][1];
    discbuff[idx].zz = points[idx][2];
    discbuff[idx].radius = spacing[idx];
    normalsbuff[idx].xx = normals[idx][0];
    normalsbuff[idx].yy = normals[idx][1];
    normalsbuff[idx].zz = normals[idx][2];
  }

  rtcCommitGeometry(geometry);
  auto scene = rtcNewScene(device);
  auto bbquality = RTC_BUILD_QUALITY_HIGH;
  rtcSetSceneBuildQuality(scene, bbquality);
  //rtcSetSceneFlags(scene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION); // does not seem to be necessary
  auto geometryID = rtcAttachGeometry(scene, geometry);
  rtcCommitScene(scene);

  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);
  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0.5; rayhit.ray.org_y = 0.5; rayhit.ray.org_z = 2;
  rayhit.ray.dir_x = 0; rayhit.ray.dir_y = 0; rayhit.ray.dir_z = -1;
  for (size_t idx = 0; idx < 8; ++idx) {
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    //rtcOccluded1(scene, &context, &rayhit.ray);
    rtcIntersect1(scene, &context, &rayhit);
    std::cout << "rayhit.ray.tfar == " << rayhit.ray.tfar << std::endl;
  }
  rtcReleaseGeometry(geometry);
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);
}
