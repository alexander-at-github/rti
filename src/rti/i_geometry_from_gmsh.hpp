#pragma once

namespace rti {
  // Interface
  class i_geometry_from_gmsh {
  public:
    // virtual destructor
    virtual ~i_geometry_from_gmsh() {};
    virtual void invert_surface_normals() = 0;
    // TODO: Remove the function invert_surface_normals() from this interface.
    // The function should be part of a gmsh class
    virtual std::string to_string() = 0;
    virtual std::string prim_to_string(unsigned int) = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual std::string get_input_file_path() = 0;
  };
} // namespace rti
