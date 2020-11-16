#pragma once

//#include <ostream>

#include <embree3/rtcore.h>

#include "absc_boundary.hpp"
#include "../rng/i_rng.hpp"
#include "../util/logger.hpp"
#include "../util/utils.hpp"

namespace rti { namespace geo {
    
  enum class bound_cond {
    REFLECTIVE, PERIODIC
  };

  template<typename Ty>
  class boundary_x_y : public absc_boundary<Ty> {
  public:

    boundary_x_y
    (RTCDevice& pDevice,
     util::pair<util::triple<Ty> >& pBdBox,
     bound_cond pXCond,
     bound_cond pYCond) :
      mDevice(pDevice),
      mBdBox(pBdBox),
      mXCond(pXCond),
      mYCond(pYCond) {
      init_this();
    }

    boundary_x_y
    (RTCDevice& pDevice,
     util::pair<util::triple<Ty> >& pBdBox) :
      boundary_x_y(pDevice, pBdBox, bound_cond::REFLECTIVE, bound_cond::REFLECTIVE) {}


    
    RTCDevice& get_rtc_device() override final
    {
      return mDevice;
    }

    RTCGeometry& get_rtc_geometry() override final
    {
      return mGeometry;
    }

    void print(std::ostream&pOs) override final
    {
      assert(false && "Not implemented");
    }

    util::triple<Ty>& get_normal_ref(unsigned int pPrimID) override final
    {
      return mNormals[pPrimID];
    }

    util::triple<Ty> get_normal(unsigned int pPrimID) override final
    {
      return mNormals[pPrimID];
    }

    std::vector<util::triple<Ty> > get_vertices() override final
    {
      auto result = std::vector<util::triple<Ty> > {mNumVertices};
      //std::cerr << "Debug: get_vertices() result.size() == " << result.size() << std::endl;
      for (size_t idx = 0; idx < mNumVertices; ++idx) {
        auto vv = mVertBuff[idx];
        result[idx] = {vv.xx, vv.yy, vv.zz};
      }
      return result;
    }

    std::vector<util::triple<size_t> > get_triangles() override final
    {
      auto result = std::vector<util::triple<size_t> > {mNumTriangles};
      //std::cerr << "Debug: get_triangles() result.size() == " << result.size() << std::endl;
      for (size_t idx = 0; idx < mNumTriangles; ++idx) {
        result[idx] = get_triangle(idx);
      }
      return result;
    }

    std::vector<util::triple<Ty> > get_triangle_normals() override final
    {
      return mNormals;
    }

    util::pair<util::triple<Ty> >
    process_hit(RTCRay& rayin, RTCHit& hitin, rng::i_rng& rng, rng::i_rng::i_state& rngstate)
    {
      assert(false && "not implemented; TODO");
      return {};
    }

  private:

    void init_this()
    {
      this->mGeometry = rtcNewGeometry(mDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
      mNumVertices = 8;
      mNumTriangles = 8;
      mVertBuff = (vertex_f3_t*) rtcSetNewGeometryBuffer
        (mGeometry,
         RTC_BUFFER_TYPE_VERTEX,
         0, // the slot
         RTC_FORMAT_FLOAT3,
         sizeof(vertex_f3_t),
         mNumVertices);
      assert(mVertBuff != nullptr && "Error in acquiring new buffer from Embree");
      auto deverr1 = rtcGetDeviceError(mDevice);
      RLOG_DEBUG << "rtc device error: " << deverr1 << std::endl;
      RLOG_DEBUG << "where RTC_ERROR_UNKOWN == " << RTC_ERROR_UNKNOWN << " holds."
                 << std::endl;
      RLOG_DEBUG << "where RTC_ERROR_INVALID_ARGUMENT == " << RTC_ERROR_INVALID_ARGUMENT
                 << " holds." << std::endl;
      mTriBuff = (triangle_t*) rtcSetNewGeometryBuffer
        (mGeometry,
         RTC_BUFFER_TYPE_INDEX,
         0, //slot
         RTC_FORMAT_UINT3,
         sizeof(triangle_t),
         mNumTriangles);
      assert(mTriBuff != nullptr && "Error in acquiring new buffer from Embree");
      //assert(false && "test");
      auto deverr2 = rtcGetDeviceError(mDevice);
      RLOG_DEBUG << "rtc device error: " << deverr2 << std::endl;
      RLOG_DEBUG << "where RTC_ERROR_UNKOWN == " << RTC_ERROR_UNKNOWN
                 << " holds." << std::endl;
      RLOG_DEBUG << "where RTC_ERROR_INVALID_ARGUMENT == " << RTC_ERROR_INVALID_ARGUMENT
                 << " holds." << std::endl;

      // Fill the vertiex
      auto xmin = std::min(mBdBox[0][0], mBdBox[1][0]);
      auto xmax = std::max(mBdBox[0][0], mBdBox[1][0]);
      auto ymin = std::min(mBdBox[0][1], mBdBox[1][1]);
      auto ymax = std::max(mBdBox[0][1], mBdBox[1][1]);
      auto zmin = std::min(mBdBox[0][2], mBdBox[1][2]);
      auto zmax = std::max(mBdBox[0][2], mBdBox[1][2]);
      mVertBuff[0].xx = xmax; mVertBuff[0].yy = ymin; mVertBuff[0].zz = zmin;
      mVertBuff[1].xx = xmax; mVertBuff[1].yy = ymin; mVertBuff[1].zz = zmax;
      mVertBuff[2].xx = xmax; mVertBuff[2].yy = ymax; mVertBuff[2].zz = zmin;
      mVertBuff[3].xx = xmax; mVertBuff[3].yy = ymax; mVertBuff[3].zz = zmax;
      mVertBuff[4].xx = xmin; mVertBuff[4].yy = ymin; mVertBuff[4].zz = zmin;
      mVertBuff[5].xx = xmin; mVertBuff[5].yy = ymin; mVertBuff[5].zz = zmax;
      mVertBuff[6].xx = xmin; mVertBuff[6].yy = ymax; mVertBuff[6].zz = zmin;
      mVertBuff[7].xx = xmin; mVertBuff[7].yy = ymax; mVertBuff[7].zz = zmax;
      // Fill the triangles
      mTriBuff[0].v0 = 0; mTriBuff[0].v1 = 1; mTriBuff[0].v2 = 2;
      mTriBuff[1].v0 = 3; mTriBuff[1].v1 = 2; mTriBuff[1].v2 = 1;
      //
      mTriBuff[2].v0 = 4; mTriBuff[2].v1 = 5; mTriBuff[2].v2 = 0;
      mTriBuff[3].v0 = 1; mTriBuff[3].v1 = 0; mTriBuff[3].v2 = 5;
      //
      mTriBuff[4].v0 = 6; mTriBuff[4].v1 = 7; mTriBuff[4].v2 = 4;
      mTriBuff[5].v0 = 5; mTriBuff[5].v1 = 4; mTriBuff[5].v2 = 7;
      //
      mTriBuff[6].v0 = 2; mTriBuff[6].v1 = 3; mTriBuff[6].v2 = 6;
      mTriBuff[7].v0 = 7; mTriBuff[7].v1 = 6; mTriBuff[7].v2 = 3;

      for (size_t idx = 0; idx < mNumTriangles; ++idx) {
        auto triangle = get_triangle_with_coords(idx);
        auto normal = util::compute_normal(triangle);
        util::normalize(normal);
        mNormals.push_back(normal);
      }
      mNormals.shrink_to_fit();
      rtcCommitGeometry(mGeometry);
      assert (RTC_ERROR_NONE == rtcGetDeviceError(mDevice) &&
              "Embree device error after rtcSetNewGeometryBuffer()");
    }

    util::triple<util::triple<Ty> > get_triangle_with_coords(size_t idx)
    {
      auto tt = mTriBuff[idx];
      return {mVertBuff[tt.v0].xx, mVertBuff[tt.v0].yy, mVertBuff[tt.v0].zz,
              mVertBuff[tt.v1].xx, mVertBuff[tt.v1].yy, mVertBuff[tt.v1].zz,
              mVertBuff[tt.v2].xx, mVertBuff[tt.v2].yy, mVertBuff[tt.v2].zz};
    }

    util::triple<size_t> get_triangle(size_t idx)
    {
      auto tt = mTriBuff[idx];
      return {tt.v0, tt.v1, tt.v2};
    }
    

  private:
    
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

    RTCDevice& mDevice;
    RTCGeometry mGeometry;
    vertex_f3_t* mVertBuff = nullptr;
    triangle_t* mTriBuff = nullptr;
    std::vector<util::triple<Ty> > mNormals;
    util::pair<util::triple<Ty> > mBdBox;
    size_t mNumVertices = 0;
    size_t mNumTriangles = 0;

    const bound_cond mXCond;
    const bound_cond mYCond;
  };
}}
