#pragma once

//#include <ostream>

#include <embree3/rtcore.h>

#include "rti/geo/i_boundary.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class boundary_x_y : public rti::geo::i_boundary<Ty> {
  public:

    boundary_x_y (RTCDevice& pDevice, rti::util::pair<rti::util::triple<Ty> > pBdBox) :
      mDevice(pDevice),
      mBdBox(pBdBox) {
      // For the time beeing the following assumption is not needed.
      // assert(mBdBox[0][0] <= mBdBox[1][0] &&
      //        mBdBox[0][1] <= mBdBox[1][1] &&
      //        mBdBox[0][2] <= mBdBox[1][2] &&
      //        "ordering of coordinates in bounding box");
      init_this();
    }

    RTCDevice& get_rtc_device() override final {
      return mDevice;
    }

    RTCGeometry& get_rtc_geometry() override final {
      return mGeometry;
    }

    void print(std::ostream&pOs) const override final {
      assert(false && "Not implemented");
    }

    rti::util::triple<Ty> get_normal(unsigned int pPrimID) const override final {
      return mNormals[pPrimID];
    }

    rti::util::triple<Ty> get_new_origin(unsigned int primID) const override final {
      assert(false && "Not implemented");
    }

  private:

    // Member variables
    RTCDevice& mDevice;
    RTCGeometry mGeometry;
    std::vector<rti::util::triple<Ty> > mNormals;
    rti::util::pair<rti::util::triple<Ty> > mBdBox;

    // Member functions
    void init_this() {
      this->mGeometry = rtcNewGeometry(mDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
      struct vertex_f3_t { // vertex is the nomenclature of Embree
        // The triangle geometry has a vertex buffer which uses x, y, and z
        // in single precision floating point types.
        float xx, yy, zz;
      };
      struct triangle_t {
        // The triangle geometry uses an index buffer that contains an array
        // of three 32-bit indices per triangle.
        uint32_t v0, v1, v2;
      };
      auto numVertices = 8;
      auto numTriangles = 8;
      auto vertBuff = (vertex_f3_t*)
        rtcSetNewGeometryBuffer(mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // the slot
                                RTC_FORMAT_FLOAT3,
                                sizeof(vertex_f3_t),
                                numVertices);
      auto triBuff = (triangle_t*)
        rtcSetNewGeometryBuffer(mGeometry,
                                RTC_BUFFER_TYPE_INDEX,
                                0, //slot
                                RTC_FORMAT_UINT3,
                                sizeof(triangle_t),
                                numTriangles);
      // Fill the vertiex
      auto xmin = std::min(mBdBox[0][0], mBdBox[1][0]);
      auto xmax = std::max(mBdBox[0][0], mBdBox[1][0]);
      auto ymin = std::min(mBdBox[0][1], mBdBox[1][1]);
      auto ymax = std::max(mBdBox[0][1], mBdBox[1][1]);
      auto zmin = std::min(mBdBox[0][2], mBdBox[1][2]);
      auto zmax = std::max(mBdBox[0][2], mBdBox[1][2]);
      vertBuff[0].xx = xmin; vertBuff[0].yy = ymin; vertBuff[0].zz = zmin;
      vertBuff[1].xx = xmax; vertBuff[1].yy = ymin; vertBuff[1].zz = zmin;
      vertBuff[2].xx = xmin; vertBuff[2].yy = ymax; vertBuff[2].zz = zmin;
      vertBuff[3].xx = xmax; vertBuff[3].yy = ymax; vertBuff[3].zz = zmin;
      vertBuff[4].xx = xmin; vertBuff[4].yy = ymin; vertBuff[4].zz = zmax;
      vertBuff[5].xx = xmax; vertBuff[5].yy = ymin; vertBuff[5].zz = zmax;
      vertBuff[6].xx = xmin; vertBuff[6].yy = ymax; vertBuff[6].zz = zmax;
      vertBuff[7].xx = xmax; vertBuff[7].yy = ymax; vertBuff[7].zz = zmax;
      // Fill the triangles
      triBuff[0].v0 = 0; triBuff[0].v1 = 1; triBuff[0].v2 = 2;
      triBuff[1].v0 = 3; triBuff[1].v1 = 2; triBuff[1].v2 = 1;
      //
      triBuff[2].v0 = 4; triBuff[2].v1 = 5; triBuff[2].v2 = 0;
      triBuff[3].v0 = 1; triBuff[3].v1 = 0; triBuff[3].v2 = 5;
      //
      triBuff[4].v0 = 6; triBuff[4].v1 = 7; triBuff[4].v2 = 4;
      triBuff[5].v0 = 5; triBuff[5].v1 = 4; triBuff[5].v2 = 7;
      //
      triBuff[6].v0 = 2; triBuff[6].v1 = 3; triBuff[6].v2 = 6;
      triBuff[7].v0 = 7; triBuff[7].v1 = 6; triBuff[7].v2 = 3;

      // TODO: fill normals
      for (size_t idx = 0; idx < numTriangles; ++idx) {
        auto tt = triBuff[idx];
        auto triangle = rti::util::triple<rti::util::triple<Ty> >
          { vertBuff[tt.v0].xx, vertBuff[tt.v0].yy, vertBuff[tt.v0].zz,
            vertBuff[tt.v1].xx, vertBuff[tt.v1].yy, vertBuff[tt.v1].zz,
            vertBuff[tt.v2].xx, vertBuff[tt.v2].yy, vertBuff[tt.v2].zz};
        auto normal = rti::util::compute_normal(triangle);
        rti::util::normalize(normal);
        mNormals.push_back(normal);
      }
      mNormals.shrink_to_fit();

      rtcCommitGeometry(mGeometry);
      assert (RTC_ERROR_NONE == rtcGetDeviceError(mDevice) &&
              "Embree device error after rtcSetNewGeometryBuffer()");
    }
  };
}} // namespace
