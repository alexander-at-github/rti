#pragma once

#include "rti/i_geometry.hpp"

namespace rti {
class absc_geometry_from_gmsh : public i_geometry {
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
  std::string get_input_file_path() override {
    return mGmshReader.get_mesh_file_path();
  }
  protected:
  RTCDevice mDevice;
  gmsh_reader& mGmshReader;
  RTCGeometry mGeometry;
  absc_geometry_from_gmsh(RTCDevice& pDevice, gmsh_reader& pGmshReader) :
    mDevice(pDevice),
    mGmshReader(pGmshReader) {}
};
} // namespace rti
