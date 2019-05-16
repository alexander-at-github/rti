#pragma once

#include "rti/i_geometry_from_gmsh.hpp"

namespace rti {
class absc_geometry_from_gmsh : public i_geometry_from_gmsh {
  public:
  virtual ~absc_geometry_from_gmsh() {}
  //virtual void invert_surface_normals() = 0;
  //virtual std::string to_string() = 0;
  //virtual std::string prim_to_string(unsigned int) = 0;
  RTCDevice& get_rtc_device() override {
    return mDevice;
  }
  RTCGeometry& get_rtc_geometry() override {
    return mGeometry;
  }
  protected:
  RTCDevice mDevice;
  RTCGeometry mGeometry;
  absc_geometry_from_gmsh(RTCDevice pDevice) :
    mDevice(pDevice) {}
};
} // namespace rti
