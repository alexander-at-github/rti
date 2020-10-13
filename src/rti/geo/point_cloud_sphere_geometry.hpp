#pragma once

#include "absc_point_cloud_geometry.hpp"
#include "../io/i_point_cloud_reader.hpp"
#include "../util/utils.hpp"

namespace rti { namespace geo {
  // The type Ty is supposed to be a numeric type - float or double.
  template<typename Ty>
  class point_cloud_sphere_geometry : public rti::geo::absc_point_cloud_geometry<Ty> {
  public:

    point_cloud_sphere_geometry(RTCDevice& pDevice, rti::io::i_point_cloud_reader<Ty>& pGReader, Ty pStickingC) :
      rti::geo::absc_point_cloud_geometry<Ty>(pDevice, pGReader, pStickingC) {
      init_this(pDevice, pGReader);
    }

    // point_cloud_sphere_geometry(point_cloud_sphere_geometry& pOther) = default;
    // point_cloud_sphere_geometry(point_cloud_sphere_geometry&& pOther) = default;
    // point_cloud_sphere_geometry& operator= (point_cloud_sphere_geometry& pOther) = default;
    // point_cloud_sphere_geometry& operator= (point_cloud_sphere_geometry&& pOther) = default;

    // needs explicit destructor? No
    //~point_cloud_sphere_geometry() override final {}

    std::string prim_to_string(unsigned int pPrimID) const override final {
      std::stringstream strstream;
      strstream
        << "(" << this->mVVBuffer[pPrimID].xx << " " << this->mVVBuffer[pPrimID].yy
        << " " << this->mVVBuffer[pPrimID].zz << " " << this->mVVBuffer[pPrimID].radius << ")";
      return strstream.str();
    }

    rti::util::quadruple<Ty> get_prim(unsigned int pPrimID) const override final {
     auto pnt = this->mVVBuffer[pPrimID];
     return {pnt.xx, pnt.yy, pnt.zz, pnt.radius};
    }

    rti::util::triple<Ty> get_normal(unsigned int pPrimID) const override final {
      return mNormals[pPrimID];
    }

    rti::util::triple<Ty> get_new_origin(RTCRay& pRay, unsigned int pPrimID) const override final {
      // take the origin of the sphere
      auto sphere = this->mVVBuffer[pPrimID];
      // ... and move in the direction of the surface normal
      auto normal = mNormals[pPrimID];
      assert(rti::util::is_normalized(normal) && "Condition: Surface normal is normalized");
      // ... such that the new origin is just above the sphere
      auto scale = sphere.radius;
      auto xx = sphere.xx + normal[0] * scale;
      auto yy = sphere.yy + normal[1] * scale;
      auto zz = sphere.zz + normal[2] * scale;
      return {xx, yy, zz};
    }

  private:
    std::vector<rti::util::triple<Ty> > mNormals;

    // Functions
    void init_this(RTCDevice& pDevice, rti::io::i_point_cloud_reader<Ty>& pGReader) {
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

      std::vector<rti::util::quadruple<Ty> > points = pGReader.get_points();
      this->mNumPoints = points.size();
      // Acquire memory via Embree API
      this->mVVBuffer = (rti::geo::point_4f_t*)
        rtcSetNewGeometryBuffer(this->mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // slot
                                RTC_FORMAT_FLOAT4,
                                sizeof(rti::geo::point_4f_t),
                                this->mNumPoints);
      // Write points to Embree data structure
      for (size_t idx = 0; idx < this->mNumPoints; ++idx) {
        rti::util::quadruple<Ty>& qudtrpl = points[idx];
        // Here we have to cast to float because:
        // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
        // single precision x, y, z floating point coordinates ..."
        // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
        this->mVVBuffer[idx].xx = (float) qudtrpl[0];
        this->mVVBuffer[idx].yy = (float) qudtrpl[1];
        this->mVVBuffer[idx].zz = (float) qudtrpl[2];
        this->mVVBuffer[idx].radius = (float) qudtrpl[3]; // TODO: Might need adjustment
      }
      this->mNormals = pGReader.get_normals(); // TODO: Deep or shallow copy
      assert(rti::util::each_normalized<Ty>(mNormals) && "Condition: surface normals are normalized");
      rtcCommitGeometry(this->mGeometry);
      assert (RTC_ERROR_NONE == rtcGetDeviceError(pDevice) &&
              "Embree device error after rtcSetNewGeometryBuffer()");
    }
  };
}} // namespace
