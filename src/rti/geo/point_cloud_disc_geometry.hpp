#pragma once

#include <cassert>
#include <sstream>

#include "absc_point_cloud_geometry.hpp"
#include "../io/i_point_cloud_reader.hpp"
#include "../util/utils.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class point_cloud_disc_geometry : public rti::geo::absc_point_cloud_geometry<Ty> {
  public:

    point_cloud_disc_geometry(RTCDevice& pDevice, rti::io::i_point_cloud_reader<Ty>& pGReader) :
      rti::geo::absc_point_cloud_geometry<Ty>(pDevice, pGReader) {
      init_this(pDevice, pGReader);
    }

    point_cloud_disc_geometry(RTCDevice& device,
                              std::vector<rti::util::quadruple<Ty> > points,
                              std::vector<rti::util::triple<Ty> > normals) :
      rti::geo::absc_point_cloud_geometry<Ty>(device) {
      init_this(device, points, normals);
    }

    std::string prim_to_string(unsigned int pPrimID) override final
    {
      auto strs = std::stringstream {};
      assert(false && "Not implemented");
      return "prim_to_string() not implemented";
    }

    rti::util::quadruple<Ty>& get_prim_ref(unsigned int pPrimID) override final
    {
      return *reinterpret_cast<rti::util::quadruple<Ty>* > (&this->mVVBuffer[pPrimID]);
    }

    rti::util::quadruple<Ty> get_prim(unsigned int pPrimID) override final
    {
      auto pnt = this->mVVBuffer[pPrimID];
      return {pnt.xx, pnt.yy, pnt.zz, pnt.radius};
    }

    rti::util::triple<Ty>& get_normal_ref(unsigned int pPrimID) override final
    {
      return *reinterpret_cast<rti::util::triple<Ty>* > (&this->mNNBuffer[pPrimID]);
    }

    rti::util::triple<Ty> get_normal(unsigned int pPrimID) override final
    {
      auto nml = this->mNNBuffer[pPrimID];
      return {(Ty) nml.xx, (Ty) nml.yy, (Ty) nml.zz};
    }

    rti::util::triple<Ty> get_new_origin(RTCRay& pRay, unsigned int pPrimID) override final {
      //auto epsilon = 1e-12; // magic number; less than 1e-3 does definitely not work
      auto xx = pRay.org_x + pRay.dir_x * pRay.tfar;
      auto yy = pRay.org_y + pRay.dir_y * pRay.tfar;
      auto zz = pRay.org_z + pRay.dir_z * pRay.tfar;
      // // add a small epsilon from the surface to be sure to be above the surface
      // auto normal = get_normal(pPrimID);
      // xx += normal[0] * epsilon;
      // yy += normal[1] * epsilon;
      // zz += normal[2] * epsilon;
      return {(Ty) xx, (Ty) yy, (Ty) zz};
    }

    Ty get_area(unsigned int primID) override final
    {
      auto radius = this->mVVBuffer[primID].radius;
      return radius * radius * rti::util::pi();
    }

  private:
    // "RTC_GEOMETRY_TYPE_POINT:
    // [...] the normal buffer stores a single precision normal per control
    // vertex (x, y, z order and RTC_FORMAT_FLOAT3 format)."
    // Source: https://embree.github.io/api.html#rtc_geometry_type_point
    struct normal_vec_3f_t {
      float xx, yy, zz;
    };
    normal_vec_3f_t* mNNBuffer = nullptr;

    void init_this(RTCDevice& device, rti::io::i_point_cloud_reader<Ty>& preader)
    {
      auto points = preader.get_points();
      auto normals = preader.get_normals();
      init_this(device, points, normals);
    }

    void init_this(RTCDevice& device,
                   std::vector<rti::util::quadruple<Ty> > points,
                   std::vector<rti::util::triple<Ty> > normals)
    {
      assert(rti::util::each_normalized<Ty>(normals) && "Condition: surface normals are normalized violated");
      assert(normals.size() == points.size() && "Number of surface normals not consistent with number of points");
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

      this->mGeometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT); // using a ray facing disc

      this->mNumPoints = points.size();

      this->mVVBuffer = (rti::geo::point_4f_t*) // type specified in the file absc_point_cloud_geometry.hpp
        rtcSetNewGeometryBuffer(this->mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // slot
                                RTC_FORMAT_FLOAT4,
                                sizeof(rti::geo::point_4f_t),
                                this->mNumPoints);

      for (size_t idx = 0; idx < this->mNumPoints; ++idx) {
        rti::util::quadruple<Ty>& qudtrpl = points[idx];
        // Here we have to cast to float because:
        // "RTC_GEOMETRY_TYPE_POINT:
        // [...] the normal buffer stores a single precision normal per control
        // vertex (x, y, z order and RTC_FORMAT_FLOAT3 format)."
        // Source: https://embree.github.io/api.html#rtc_geometry_type_point
        this->mVVBuffer[idx].xx = (float) qudtrpl[0];
        this->mVVBuffer[idx].yy = (float) qudtrpl[1];
        this->mVVBuffer[idx].zz = (float) qudtrpl[2];
        this->mVVBuffer[idx].radius = (float) qudtrpl[3];
      }
      this->mNNBuffer = (normal_vec_3f_t*)
        rtcSetNewGeometryBuffer(this->mGeometry,
                                RTC_BUFFER_TYPE_NORMAL,
                                0, // slot
                                RTC_FORMAT_FLOAT3,
                                sizeof(normal_vec_3f_t),
                                this->mNumPoints);
      for (size_t idx = 0; idx < this->mNumPoints; ++idx) {
        this->mNNBuffer[idx].xx = normals[idx][0];
        this->mNNBuffer[idx].yy = normals[idx][1];
        this->mNNBuffer[idx].zz = normals[idx][2];
      }

      rtcCommitGeometry(this->mGeometry);
      assert (RTC_ERROR_NONE == rtcGetDeviceError(device) &&
              "Embree device error after rtcSetNewGeometryBuffer()");
    }
  };
}}
