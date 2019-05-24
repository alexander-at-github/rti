#pragma once

#include <cmath>

#include "rti/gmsh_reader.hpp"
#include "rti/i_geometry.hpp"
#include "rti/types.hpp"

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
  gmsh_reader& mGmshReader;
  RTCDevice mDevice;
  RTCGeometry mGeometry;
  absc_geometry_from_gmsh(RTCDevice& pDevice, gmsh_reader& pGmshReader) :
    mGmshReader(pGmshReader),
    mDevice(pDevice) {}

  // Input an array of triples where each tripple are the x, y and z coordinates of a vertex
  // of a triangle
  static rti::triple_t<float> centroid(rti::triple_t<rti::triple_t<float> > pTriangle) {
    rti::triple_t<float> result;
    for (auto idx : {0, 1, 2}) {
      result[idx] = (pTriangle[0][idx] + pTriangle[1][idx] + pTriangle[2][idx]) / 3;
    }
    return result;
  }
  static float distance(rti::pair_t<rti::triple_t<float> > pPnts) {
    auto p1 = pPnts[0][0] - pPnts[1][0];
    auto p2 = pPnts[0][1] - pPnts[1][1];
    auto p3 = pPnts[0][2] - pPnts[1][2];
    return std::sqrt(p1*p1 + p2*p2 + p3*p3);
  }
};
} // namespace rti
