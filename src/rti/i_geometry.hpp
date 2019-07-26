#pragma once

namespace rti {
  // Interface
  template<typename Ty>
  class i_geometry {
  public:
    // virtual destructor
    virtual ~i_geometry() {}
    virtual void print(std::ostream& pOs) const = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual std::string get_input_file_path() = 0;
    virtual std::string prim_to_string(unsigned int primID) const = 0;
    //virtual rti::triple<Ty> get_coords(unsigned int primID) const = 0;
    virtual rti::triple<Ty> get_normal(unsigned int primID) const = 0;
    virtual rti::triple<Ty> get_new_origin(unsigned int primID) const = 0;
    virtual rti::pair<rti::triple<Ty> > get_bounding_box() const = 0;
  };
} // namespace rti
