#pragma once

namespace rti {
  // Interface
  class i_geometry {
  public:
    // virtual destructor
    virtual ~i_geometry() {};
    virtual std::string to_string() = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual std::string get_input_file_path() = 0;
    virtual std::string prim_to_string(unsigned int) = 0;
  };
} // namespace rti
