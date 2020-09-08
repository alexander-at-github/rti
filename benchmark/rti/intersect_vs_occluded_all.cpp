#include <benchmark/benchmark.h>

#include <embree3/rtcore.h>
#include <immintrin.h> // AVX
#include <pmmintrin.h>
#include <xmmintrin.h>

#include <random>

#include "rti/trace/local_intersector.hpp"

struct point_4f_t {
  float xx, yy, zz, radius;
};
struct normal_vec_3f_t {
  float xx, yy, zz;
};

using locint_type = rti::trace::local_intersector<float, point_4f_t, normal_vec_3f_t>;

using embree_record_t = struct {
  RTCDevice dev;
  RTCScene scn;
  RTCGeometry geo;
  locint_type locint;
};

static
void occluded_filter(const struct RTCFilterFunctionNArguments* args)
{
  auto valid = args->valid; // a pointer
  valid[0] = 0; // set invalid; causes continuation of ray
  //std::cerr << "Occluded filter. primID == " << (reinterpret_cast<RTCHit*> (args->hit))->primID << std::endl;
}

static
void intersect_filter(const struct RTCFilterFunctionNArguments* args)
{
  auto valid = args->valid; // a pointer
  valid[0] = 0; // set invalid; causes continuation of ray
}

static
float get_random_pertubation()
{
  auto static generator = std::mt19937 (std::random_device{}());
  auto static distribution = std::uniform_real_distribution<float> (-0.1, 0.1);
  return distribution(generator);
}

static
embree_record_t prepare_embree() {
  {  // if multi-threaded, then this block needs to be parallelized.
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  }
  auto numpoints = 8u; // unsigned int
  auto device_config = "hugepages=1";
  auto device = rtcNewDevice(device_config);
  assert(rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_FILTER_FUNCTION_SUPPORTED) != 0 && "Correctness Assumption");
  auto geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

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
    discs[idx].xx = get_random_pertubation(); discs[idx].yy = get_random_pertubation(); discs[idx].zz = 100 - idx*10 + get_random_pertubation();
    discs[idx].radius = 1;
    normals[idx].xx = get_random_pertubation(); normals[idx].yy = get_random_pertubation(); normals[idx].zz = -1;
  }
  rtcSetGeometryOccludedFilterFunction(geometry, &occluded_filter);
  rtcSetGeometryIntersectFilterFunction(geometry, &intersect_filter);
  rtcCommitGeometry(geometry);
  auto scene = rtcNewScene(device);
  auto bbquality = RTC_BUILD_QUALITY_HIGH;
  rtcSetSceneBuildQuality(scene, bbquality);
  //rtcSetSceneFlags(scene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION); // does not seem to be necessary
  auto geometryID = rtcAttachGeometry(scene, geometry);
  rtcCommitScene(scene);

  return {device, scene, geometry, locint_type (numpoints, discs, normals)};
}

static
void release_embree(embree_record_t pEm) {
  rtcReleaseGeometry(pEm.geo);
  rtcReleaseScene(pEm.scn);
  rtcReleaseDevice(pEm.dev);
}

static
void intersect_all(benchmark::State& pState) {
  // prepare
  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);
  auto embreerec = prepare_embree();
  auto scene = embreerec.scn;
  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0; rayhit.ray.org_y = 0; rayhit.ray.org_z = 0;
  rayhit.ray.dir_x = 0; rayhit.ray.dir_y = 0; rayhit.ray.dir_z = 1;
  // benchmark loop
  for (auto _ : pState) {
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    //rtcOccluded1(scene, &context, &rayhit.ray);
    rtcIntersect1(scene, &context, &rayhit);

    // Do we need DoNotOptimize here?
    benchmark::DoNotOptimize(rayhit);
  }
  // clean up
  release_embree(embreerec);
}

static
void occluded_all(benchmark::State& pState) {
  // prepare
  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);
  auto embreerec = prepare_embree();
  auto scene = embreerec.scn;
  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0; rayhit.ray.org_y = 0; rayhit.ray.org_z = 0;
  rayhit.ray.dir_x = 0; rayhit.ray.dir_y = 0; rayhit.ray.dir_z = 1;
  // benchmark loop
  for (auto _ : pState) {
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    rtcOccluded1(scene, &context, &rayhit.ray);
    //rtcIntersect1(scene, &context, &rayhit);

    // Do we need DoNotOptimize here?
    benchmark::DoNotOptimize(rayhit);
  }
  // clean up
  release_embree(embreerec);
}

static
void intersect_first(benchmark::State& pState) {
  // prepare
  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);
  auto embreerec = prepare_embree();
  auto scene = embreerec.scn;
  // remove the filters
  rtcSetGeometryOccludedFilterFunction(embreerec.geo, NULL);
  rtcSetGeometryIntersectFilterFunction(embreerec.geo, NULL);
  //
  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0; rayhit.ray.org_y = 0; rayhit.ray.org_z = 0;
  rayhit.ray.dir_x = 0; rayhit.ray.dir_y = 0; rayhit.ray.dir_z = 1;
  // benchmark loop
  for (auto _ : pState) {
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    //rtcOccluded1(scene, &context, &rayhit.ray);
    rtcIntersect1(scene, &context, &rayhit);

    // Do we need DoNotOptimize here?
    benchmark::DoNotOptimize(rayhit);
  }
  // clean up
  release_embree(embreerec);
}

static
void occluded_any_one(benchmark::State& pState) {
  // prepare
  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);
  auto embreerec = prepare_embree();
  auto scene = embreerec.scn;
  // remove the filters
  rtcSetGeometryOccludedFilterFunction(embreerec.geo, NULL);
  rtcSetGeometryIntersectFilterFunction(embreerec.geo, NULL);
  //
  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0; rayhit.ray.org_y = 0; rayhit.ray.org_z = 0;
  rayhit.ray.dir_x = 0; rayhit.ray.dir_y = 0; rayhit.ray.dir_z = 1;
  // benchmark loop
  for (auto _ : pState) {
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    rtcOccluded1(scene, &context, &rayhit.ray);
    //rtcIntersect1(scene, &context, &rayhit);

    // Do we need DoNotOptimize here?
    benchmark::DoNotOptimize(rayhit);
  }
  // clean up
  release_embree(embreerec);
}

static
void intersect_all_with_local_intersector(benchmark::State& pState) {
  // prepare
  auto context = RTCIntersectContext {};
  rtcInitIntersectContext(&context);
  auto embreerec = prepare_embree();
  auto scene = embreerec.scn;
  auto localintersector = embreerec.locint;
  // remove the filters
  rtcSetGeometryOccludedFilterFunction(embreerec.geo, NULL);
  rtcSetGeometryIntersectFilterFunction(embreerec.geo, NULL);
  //
  auto rayhit = RTCRayHit {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  rayhit.ray.org_x = 0; rayhit.ray.org_y = 0; rayhit.ray.org_z = 0;
  rayhit.ray.dir_x = 0; rayhit.ray.dir_y = 0; rayhit.ray.dir_z = 1;
  // benchmark loop
  for (auto _ : pState) {
    rayhit.ray.tfar = std::numeric_limits<float>::max();
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    //rtcOccluded1(scene, &context, &rayhit.ray);
    rtcIntersect1(scene, &context, &rayhit);
    auto intersects = localintersector.intersect_neighbors(rayhit);
    //std::cout << intersects.size() << std::endl;

    // Do we need DoNotOptimize here?
    benchmark::DoNotOptimize(rayhit);
    benchmark::DoNotOptimize(intersects);
  }
  // clean up
  release_embree(embreerec);
}

BENCHMARK(intersect_all)->UseRealTime();
BENCHMARK(occluded_all)->UseRealTime();
BENCHMARK(intersect_first)->UseRealTime();
BENCHMARK(occluded_any_one)->UseRealTime();
BENCHMARK(intersect_all_with_local_intersector)->UseRealTime();

