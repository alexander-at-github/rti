#pragma once

#include <cassert>
#include <sstream>

#include "disc_neighborhood.hpp"
#include "i_abs_geometry.hpp"
#include "../io/i_point_cloud_reader.hpp"
#include "../util/timer.hpp"
#include "../util/utils.hpp"

namespace rti { namespace geo {
  template<typename numeric_type>
  class point_cloud_disc_geometry : public i_abs_geometry<numeric_type> {
    
  public:

    point_cloud_disc_geometry
    (RTCDevice& pDevice, io::i_point_cloud_reader<numeric_type>& pGReader) :
      mDevice(pDevice),
      mInfilename(pGReader.get_input_file_name()) {
      init_this(mDevice, pGReader);
    }

    point_cloud_disc_geometry
    (RTCDevice& pDevice,
     std::vector<util::quadruple<numeric_type> > points,
     std::vector<util::triple<numeric_type> > normals) :
      mDevice(pDevice),
      mInfilename("") {
      init_this(mDevice, points, normals);
    }

    std::string prim_to_string(unsigned int pPrimID)
    {
      assert(false && "Not implemented");
      return "prim_to_string() not implemented";
    }

    util::quadruple<numeric_type>& get_prim_ref(unsigned int pPrimID)
    {
      static_assert(std::is_same<numeric_type, float>::value,
                    "Error: Trying to use refernces of incompatible types.");
      return *reinterpret_cast<util::quadruple<numeric_type>* > (&mVVBuffer[pPrimID]);
    }

    util::quadruple<numeric_type> get_prim(unsigned int pPrimID)
    {
      auto const& pnt = mVVBuffer[pPrimID];
      return {(numeric_type) pnt.xx, (numeric_type) pnt.yy, (numeric_type) pnt.zz, (numeric_type) pnt.radius};
    }

    util::triple<numeric_type>& get_normal_ref(unsigned int pPrimID) override final
    {
      static_assert(std::is_same<numeric_type, float>::value,
                    "Error: Trying to use refernces of incompatible types.");
      return *reinterpret_cast<util::triple<numeric_type>* > (&mNNBuffer[pPrimID]);
    }

    util::triple<numeric_type> get_normal(unsigned int pPrimID) override final
    {
      auto nml = mNNBuffer[pPrimID];
      return {(numeric_type) nml.xx, (numeric_type) nml.yy, (numeric_type) nml.zz};
    }

    util::triple<numeric_type> get_new_origin(RTCRay& pRay, unsigned int pPrimID) override final
    {
      //auto epsilon = 1e-12; // magic number; less than 1e-3 does definitely not work
      auto xx = pRay.org_x + pRay.dir_x * pRay.tfar;
      auto yy = pRay.org_y + pRay.dir_y * pRay.tfar;
      auto zz = pRay.org_z + pRay.dir_z * pRay.tfar;
      // // add a small epsilon from the surface to be sure to be above the surface
      // auto normal = get_normal(pPrimID);
      // xx += normal[0] * epsilon;
      // yy += normal[1] * epsilon;
      // zz += normal[2] * epsilon;
      return {(numeric_type) xx, (numeric_type) yy, (numeric_type) zz};
    }

    RTCDevice& get_rtc_device() override final
    {
      return mDevice;
    }

    RTCGeometry& get_rtc_geometry() override final
    {
      return mGeometry;
    }

    std::string get_input_file_path()
    {
      return mInfilename;
    }

    util::pair<util::triple<numeric_type> >
    get_bounding_box()
    {
      return {mincoords, maxcoords};
    }

    size_t get_num_primitives()
    {
      return mNumPoints;
    }

    void print(std::ostream& pOs) override final
    {
      pOs << "(:class " << typeid(this).name();
      if (mVVBuffer != nullptr)
        for (size_t idx = 0; idx < mNumPoints; ++idx)
          pOs << prim_to_string(idx);
      pOs << ")";
    }

    std::vector<size_t>& get_neighbors(unsigned int id)
    {
      return discnbhd.get_neighbors(id);
    }

  private:

    void init_this
    (RTCDevice& device,
     io::i_point_cloud_reader<numeric_type>& preader)
    {
      auto points = preader.get_points();
      auto normals = preader.get_normals();
      init_this(device, points, normals);
    }

    void init_this
    (RTCDevice& device,
     std::vector<util::quadruple<numeric_type> > points,
     std::vector<util::triple<numeric_type> > normals)
    {
      //std::cout << "initializing geometry with points.size() == " << points.size() << std::endl;
      assert(util::each_normalized<numeric_type>(normals) &&
             "Condition: surface normals are normalized violated");
      assert(normals.size() == points.size() &&
             "Number of surface normals not consistent with number of points");
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

      // using a ray facing discs
      mGeometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);
      mNumPoints = points.size();

      mVVBuffer = (point_4f_t*) rtcSetNewGeometryBuffer
        (mGeometry,
         RTC_BUFFER_TYPE_VERTEX,
         0, // slot
         RTC_FORMAT_FLOAT4,
         sizeof(point_4f_t),
         mNumPoints);
      
      for (size_t idx = 0; idx < mNumPoints; ++idx) {
        util::quadruple<numeric_type> const& qudtrpl = points[idx];
        // Here we have to cast to float because:
        // "RTC_GEOMETRY_TYPE_POINT:
        // [...] the normal buffer stores a single precision normal per control
        // vertex (x, y, z order and RTC_FORMAT_FLOAT3 format)."
        // Source: https://embree.github.io/api.html#rtc_geometry_type_point
        mVVBuffer[idx].xx = (float) qudtrpl[0];
        mVVBuffer[idx].yy = (float) qudtrpl[1];
        mVVBuffer[idx].zz = (float) qudtrpl[2];
        mVVBuffer[idx].radius = (float) qudtrpl[3];
        if (qudtrpl[0] < mincoords[0]){ mincoords[0] = qudtrpl[0]; }
        if (qudtrpl[1] < mincoords[1]){ mincoords[1] = qudtrpl[1]; }
        if (qudtrpl[2] < mincoords[2]){ mincoords[2] = qudtrpl[2]; }
        if (qudtrpl[0] > maxcoords[0]){ maxcoords[0] = qudtrpl[0]; }
        if (qudtrpl[1] > maxcoords[1]){ maxcoords[1] = qudtrpl[1]; }
        if (qudtrpl[2] > maxcoords[2]){ maxcoords[2] = qudtrpl[2]; }
      }
      mNNBuffer = (normal_vec_3f_t*) rtcSetNewGeometryBuffer
        (mGeometry,
         RTC_BUFFER_TYPE_NORMAL,
         0, // slot
         RTC_FORMAT_FLOAT3,
         sizeof(normal_vec_3f_t),
         mNumPoints);
      
      for (size_t idx = 0; idx < mNumPoints; ++idx) {
        mNNBuffer[idx].xx = normals[idx][0];
        mNNBuffer[idx].yy = normals[idx][1];
        mNNBuffer[idx].zz = normals[idx][2];
      }

      rtcCommitGeometry(mGeometry);
      assert (RTC_ERROR_NONE == rtcGetDeviceError(device) &&
              "Embree device error after rtcSetNewGeometryBuffer()");

      std::cout << "Creating neighborhood ... " << std::flush;
      auto timer = util::timer {};
      // discnbhd.setup_neighborhood_naive(points);
      discnbhd.setup_neighborhood(points, mincoords, maxcoords);
      auto elapsed = timer.elapsed_seconds();
      std::cout << " took " << elapsed << " seconds" << std::endl;
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
    // "RTC_GEOMETRY_TYPE_POINT:
    // The vertex buffer stores each control vertex in the form of a single
    // precision position and radius stored in (x, y, z, r) order in memory
    // (RTC_FORMAT_FLOAT4 format). The number of vertices is inferred from the
    // size of this buffer. Similarly, the normal buffer stores a single
    // precision normal per control vertex (x, y, z order and RTC_FORMAT_FLOAT3 format)."
    // Source: https://embree.github.io/api.html#rtc_geometry_type_point
    struct point_4f_t {
      float xx, yy, zz, radius;
    };
    point_4f_t* mVVBuffer = nullptr;

    RTCDevice& mDevice;
    RTCGeometry mGeometry;
    size_t mNumPoints = 0;
    std::string mInfilename;
    geo::disc_neighborhood<numeric_type> discnbhd;

    constexpr static numeric_type nummax = std::numeric_limits<numeric_type>::max();
    constexpr static numeric_type nummin = std::numeric_limits<numeric_type>::lowest();
    util::triple<numeric_type> mincoords = util::triple<numeric_type> {nummax, nummax, nummax};
    util::triple<numeric_type> maxcoords = util::triple<numeric_type> {nummin, nummin, nummin};
  };
}}
