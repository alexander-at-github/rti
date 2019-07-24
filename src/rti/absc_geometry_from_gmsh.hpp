#pragma once

#include <cmath>

#include "rti/gmsh_reader.hpp"
#include "rti/i_geometry.hpp"
#include "rti/utils.hpp"

namespace rti {
  template<typename Ty>
  class absc_geometry_from_gmsh : public i_geometry<Ty> {
  public:

  virtual ~absc_geometry_from_gmsh() {}

  RTCDevice& get_rtc_device() override final {
    return mDevice;
  }
  RTCGeometry& get_rtc_geometry() override final {
    return mGeometry;
  }
  std::string get_input_file_path() override final {
    return mGmshReader.get_mesh_file_path();
  }
  protected:
  gmsh_reader& mGmshReader;
  RTCDevice mDevice;
  RTCGeometry mGeometry;
  absc_geometry_from_gmsh(RTCDevice& pDevice, gmsh_reader& pGmshReader) :
    mGmshReader(pGmshReader),
    mDevice(pDevice) {}

  // Input: a triple of triples where each inner triple holds the x, y and z coordinates
  // of a vertex of a triangle.
  static rti::triple<Ty> centroid(rti::triple<rti::triple<Ty> > pTriangle) {
    rti::triple<Ty> result;
    result[0] = (pTriangle[0][0] + pTriangle[1][0] + pTriangle[2][0]) / 3;
    result[1] = (pTriangle[0][1] + pTriangle[1][1] + pTriangle[2][1]) / 3;
    result[2] = (pTriangle[0][2] + pTriangle[1][2] + pTriangle[2][2]) / 3;
    return result;
  }
  static Ty distance(rti::pair<rti::triple<Ty> > pPnts) {
    auto p1 = pPnts[0][0] - pPnts[1][0];
    auto p2 = pPnts[0][1] - pPnts[1][1];
    auto p3 = pPnts[0][2] - pPnts[1][2];
    return std::sqrt(p1*p1 + p2*p2 + p3*p3);
  }
};
} // namespace rti
