#pragma once

#include <cmath>

#include "rti/geo/i_geometry.hpp"
#include "rti/io/gmsh_reader.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class absc_geometry_from_gmsh : public rti::geo::i_geometry<Ty> {
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
  rti::io::gmsh_reader& mGmshReader;
  RTCDevice mDevice;
  RTCGeometry mGeometry;
  absc_geometry_from_gmsh(RTCDevice& pDevice, rti::io::gmsh_reader& pGmshReader) :
    mGmshReader(pGmshReader),
    mDevice(pDevice) {}

  // Input: a triple of triples where each inner triple holds the x, y and z coordinates
  // of a vertex of a triangle.
  static rti::util::triple<Ty> centroid(rti::util::triple<rti::util::triple<Ty> > pTriangle) {
    rti::util::triple<Ty> result;
    result[0] = (pTriangle[0][0] + pTriangle[1][0] + pTriangle[2][0]) / 3;
    result[1] = (pTriangle[0][1] + pTriangle[1][1] + pTriangle[2][1]) / 3;
    result[2] = (pTriangle[0][2] + pTriangle[1][2] + pTriangle[2][2]) / 3;
    return result;
  }
  static Ty distance(rti::util::pair<rti::util::triple<Ty> > pPnts) {
    auto p1 = pPnts[0][0] - pPnts[1][0];
    auto p2 = pPnts[0][1] - pPnts[1][1];
    auto p3 = pPnts[0][2] - pPnts[1][2];
    return std::sqrt(p1*p1 + p2*p2 + p3*p3);
  }
};
}} // namespace
