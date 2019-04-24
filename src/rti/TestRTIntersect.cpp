#include <embree3/rtcore.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

int main(int argc, char* argv[]) {
  // Enable Flush-to-Zero and Denormals-are-Zero for the MXSCR status and
  // control registers for performance reasons.
  // See: https://embree.github.io/api.html#mxcsr-control-and-status-register
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  // Enable huge page support.
  constexpr std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config);
  RTCScene scene = rtcNewScene(device);

}
