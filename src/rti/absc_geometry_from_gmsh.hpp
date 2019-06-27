#pragma once

#include <cmath>

#include "rti/gmsh_reader.hpp"
#include "rti/i_geometry.hpp"
#include "rti/types.hpp"

namespace rti {
class absc_geometry_from_gmsh : public i_geometry {
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
  static rti::triple<float> centroid(rti::triple<rti::triple<float> > pTriangle) {
    rti::triple<float> result;
    result.frst = (pTriangle.frst.frst + pTriangle.scnd.frst + pTriangle.thrd.frst) / 3;
    result.scnd = (pTriangle.frst.scnd + pTriangle.scnd.scnd + pTriangle.thrd.scnd) / 3;
    result.thrd = (pTriangle.frst.thrd + pTriangle.scnd.thrd + pTriangle.thrd.thrd) / 3;
    return result;
  }
  static float distance(rti::pair<rti::triple<float> > pPnts) {
    auto p1 = pPnts.frst.frst - pPnts.scnd.frst;
    auto p2 = pPnts.frst.scnd - pPnts.scnd.scnd;
    auto p3 = pPnts.frst.thrd - pPnts.scnd.thrd;
    return std::sqrt(p1*p1 + p2*p2 + p3*p3);
  }
};
} // namespace rti
