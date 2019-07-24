#pragma once

#include <boost/core/demangle.hpp>

#include "rti/i_geometry.hpp"
#include "rti/i_geometry_reader.hpp"

namespace rti {
  // The type Ty is supposed to be a numeric type - float or double.
  template<typename Ty>
  class point_cloud_geometry : public i_geometry<Ty> {
  public:

    point_cloud_geometry(RTCDevice& pDevice, rti::i_geometry_reader<Ty>& pGReader) :
      mDevice(pDevice),
      mInfilename(pGReader.get_input_file_name()) {
      init_this(pDevice, pGReader);
    }

    void print(std::ostream& pOs) const override final {
      pOs << "(:class " << boost::core::demangle(typeid(this).name());
      if (mVVBuffer != nullptr)
        for (size_t idx = 0; idx < this->mNumPoints; ++idx)
          pOs << this->prim_to_string(idx);
      pOs << ")";
    }

    std::string prim_to_string(unsigned int pPrimID) const override final {
      std::stringstream strstream;
      strstream
        << "(" << mVVBuffer[pPrimID].xx << " " << mVVBuffer[pPrimID].yy
        << " " << mVVBuffer[pPrimID].zz << " " << mVVBuffer[pPrimID].radius << ")";
      return strstream.str();
    }

    rti::triple<Ty> get_normal(unsigned int pPrimID) const override final {
      return mNormals[pPrimID];
    }

    RTCDevice& get_rtc_device() override final {
      return mDevice;
    }

    RTCGeometry& get_rtc_geometry() override final {
      return mGeometry;
    }

    std::string get_input_file_path() override final {
      return mInfilename;
    }

    rti::pair<rti::triple<Ty> > get_bounding_box() const override final {
      assert(mVVBuffer != nullptr && "No data");
      if (mVVBuffer == nullptr) // no data in this instance
        return {rti::triple<Ty> {0,0,0}, rti::triple<Ty> {0,0,0}};
      Ty min = std::numeric_limits<Ty>::min();
      Ty max = std::numeric_limits<Ty>::max();
      Ty xmin=max, xmax=min, ymin=max, ymax=min, zmin=max, zmax=min;
      for (size_t idx = 0; idx < mNumPoints; ++idx) {
        xmin = std::min(xmin, mVVBuffer[idx].xx);
        xmax = std::max(xmax, mVVBuffer[idx].xx);
        ymin = std::min(ymin, mVVBuffer[idx].yy);
        ymax = std::max(ymax, mVVBuffer[idx].yy);
        zmin = std::min(zmin, mVVBuffer[idx].zz);
        zmax = std::max(zmax, mVVBuffer[idx].zz);
      }
      return {rti::triple<Ty> {xmin, ymin, zmin}, rti::triple<Ty> {xmax, ymax, zmax}};
    }

  private:
    // Local types
    struct point_4f_t {
      float xx, yy, zz, radius;
      // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
      // single precision x, y, z floating point coordinates
      // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred
      // from the size of that buffer. "
      // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
    };
    // Data members
    RTCDevice& mDevice;
    RTCGeometry mGeometry;
    point_4f_t* mVVBuffer = nullptr;
    std::vector<rti::triple<Ty> > mNormals;
    size_t mNumPoints = 0;
    std::string mInfilename;
    // Functions
    void init_this(RTCDevice& pDevice, i_geometry_reader<Ty>& pGReader) {
      // "Points with per vertex radii are supported with sphere, ray-oriented
      // discs, and normal-oriented discs geometric represetntations. Such point
      // geometries are created by passing RTC_GEOMETRY_TYPE_SPHERE_POINT,
      // RTC_GEOMETRY_TYPE_DISC_POINT, or RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT to
      // the rtcNewGeometry function."
      //
      // "RTC_GEOMETRY_TYPE_SPHERE_POINT -
      //   point geometry spheres
      // RTC_GEOMETRY_TYPE_DISC_POINT -
      //   point geometry with ray-oriented discs
      // RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT -
      //   point geometry with normal-oriented discs"
      //
      // Source: https://embree.github.io/api.html#rtc_geometry_type_point

      // "The point vertices can be specified t through a vertex buffer (RTC_BUFFER_TYPE_VERTEX).
      // For the normal oriented discs a normal buffer (RTC_BUFFER_TYPE_NORMAL) has to get specified
      // additionally. See rtcSetGeometryBuffer and rtcSetSharedGeometryBuffer for more details on
      // how to set buffers."
      // Source: https//embree.github.io/api.html#rtc_geometry_type_point

      // "Points with per vertex radii are supported with sphere, ray-oriented discs, and
      // normal-oriented discs geometric represetntations. Such point geometries are created by
      // passing RTC_GEOMETRY_TYPE_SPHERE_POINT, RTC_GEOMETRY_TYPE_DISC_POINT, or
      // RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT to the rtcNewGeometry function. The point vertices can
      // be specified t through a vertex buffer (RTC_BUFFER_TYPE_VERTEX). For the normal oriented
      // discs a normal buffer (RTC_BUFFER_TYPE_NORMAL) has to get specified additionally. See
      // rtcSetGeometryBuffer and rtcSetSharedGeometryBuffer for more details on how to set buffers.
      //
      // The vertex buffer stores each control vertex in the form of a single precision position and
      // radius stored in (x, y, z, r) order in memory (RTC_FORMAT_FLOAT4 format). The number of
      // vertices is inferred from the size of this buffer. Similarly, the normal buffer stores a
      // single precision normal per control vertex (x, y, z order and RTC_FORMAT_FLOAT3 format).
      //
      // In the RTC_GEOMETRY_TYPE_SPHERE_POINT mode, a real geometric surface is rendered for the
      // curve, which is more expensive but allows closeup views.
      //
      // The RTC_GEOMETRY_TYPE_DISC_POINT flat mode is a fast mode designed to render distant point.
      // In this mode the point is rendered as a ray facing disc.
      //
      // The RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT mode is a mode designed as a midpoint
      // geometrically between ray facing discs and spheres. In this mode the point is rendered as
      // a normal oriented disc.
      //
      // For all point types, only the hit distance and geometry normal is returned as hit
      // information, u and v are set to zero."
      // Source: https://www.embree.org/api.html#rtc_geometry_type_point

      this->mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_DISC_POINT); // using a ray facing disc

      std::vector<rti::quadruple<Ty> > points = pGReader.get_points();
      this->mNumPoints = points.size();
      // Acquire memory via Embree API
      mVVBuffer = (point_4f_t*)
        rtcSetNewGeometryBuffer(mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // slot
                                RTC_FORMAT_FLOAT4,
                                sizeof(point_4f_t),
                                this->mNumPoints);
      // Write points to Embree data structure
      for (size_t idx = 0; idx < this->mNumPoints; ++idx) {
        rti::quadruple<Ty>& qudtrpl = points[idx];
        // Here we have to cast to float because:
        // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
        // single precision x, y, z floating point coordinates ..."
        // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
        mVVBuffer[idx].xx = (float) qudtrpl[0];
        mVVBuffer[idx].yy = (float) qudtrpl[1];
        mVVBuffer[idx].zz = (float) qudtrpl[2];
        mVVBuffer[idx].radius = (float) qudtrpl[3]; // TODO: Might need adjustment
      }
      this->mNormals = pGReader.get_normals(); // TODO: Deep or shallow copy
      rtcCommitGeometry(mGeometry);
      assert (RTC_ERROR_NONE == rtcGetDeviceError(pDevice) &&
              "Embree device error after rtcSetNewGeometryBuffer()");
    }
  };
} // namespace rti
