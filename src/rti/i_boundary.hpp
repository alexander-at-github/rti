#pragma once

namespace rti {
  template<typename Ty>
  class i_boundary {
  public:
    virtual ~i_boundary() {}
    // Should the following function be implemented here?
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
  };
} // namespace rti
