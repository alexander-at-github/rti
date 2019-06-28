#pragma once

namespace rti {
  // Interface
  class i_geometry {
  public:
    // virtual destructor
    virtual ~i_geometry() {}
    virtual void print(std::ostream& pOs) const = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual std::string get_input_file_path() = 0;
    virtual std::string prim_to_string(unsigned int) const = 0;
    virtual rti::triple<float> get_normal(unsigned int primID) const = 0;
  };
} // namespace rti
