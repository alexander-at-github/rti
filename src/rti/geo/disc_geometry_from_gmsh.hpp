#pragma once

#include <sstream>

#include <embree3/rtcore.h>
#include <gmsh.h>

#include "rti/geo/absc_geometry_from_gmsh.hpp"
#include "rti/io/gmsh_reader.hpp"
#include "rti/util/logger.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  // One might want to templatize this function and remove all float type specifiers
  class disc_geometry_from_gmsh : public rti::geo::absc_geometry_from_gmsh<float> {
  public:
    disc_geometry_from_gmsh(RTCDevice& pDevice, rti::io::gmsh_reader& pGmshReader) :
      absc_geometry_from_gmsh(pDevice, pGmshReader) {
      init_this(pDevice, pGmshReader);
    }

    void print(std::ostream& pOs) const override final {
      pOs << "(:class disc_geometry_from_gmsh ";
      if (mVVBuffer != nullptr) {
        for (size_t idx = 0; idx < this->mNumVertices; ++idx) {
          pOs << this->prim_to_string(idx);
        }
      }
      pOs << ")";
    }

    std::string prim_to_string(unsigned int pPrimID) const override final {
      //assert(false && "not implemented");
      std::stringstream strstream;
      strstream
        << "(" << mVVBuffer[pPrimID].xx << "," << mVVBuffer[pPrimID].yy
        << "," << mVVBuffer[pPrimID].zz << "," << mVVBuffer[pPrimID].radius << ")";
      return strstream.str();
    }

    rti::util::triple<float> get_normal(unsigned int primID) const {
      assert(false && "Not implemented");
      rti::util::triple<float> rr;
      return rr;
    }

  private:
    ////////////////////////
    // Local algebraic types
    ////////////////////////
    struct vertex_f4_t {
      float xx, yy, zz, radius;
      // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
      // single precision x, y, z floating point coordinates
      // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred
      // from the size of that buffer. "
      // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
    };
    ///////////////
    // Data members
    ///////////////
    vertex_f4_t* mVVBuffer = nullptr;
    size_t mNumTriangles = 0;
    size_t mNumVertices = 0;
    ////////////
    // Functions
    ////////////
    void init_this(RTCDevice& pDevice, rti::io::gmsh_reader& pGmshReader) {
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

      mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_DISC_POINT); // using a ray facing disc

      // Read input from gmsh
      std::vector<rti::util::triple<double> > vertices = pGmshReader.get_vertices();

      this->mNumVertices = vertices.size();
      // Acquire memory from Embree
      mVVBuffer = (vertex_f4_t*)
        rtcSetNewGeometryBuffer(
                                mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // slot
                                RTC_FORMAT_FLOAT4,
                                sizeof(vertex_f4_t),
                                this->mNumVertices);

      // Write vertices to Embree
      for (size_t idx = 0; idx < this->mNumVertices; ++idx) {
        auto& triple = vertices[idx];
        mVVBuffer[idx].xx = triple[0];
        mVVBuffer[idx].yy = triple[1];
        mVVBuffer[idx].zz = triple[2];
        mVVBuffer[idx].radius = 0; // Update to sound value later
      }

      std::vector<rti::util::triple<std::size_t> > triangles = pGmshReader.get_triangles();

      this->mNumTriangles = triangles.size();
      // Set radii of discs
      for (size_t idx = 0; idx < this->mNumTriangles; ++idx) {
        auto& triangle = triangles[idx];
        rti::util::triple<rti::util::triple<float> > trnglCoords
          {rti::util::triple<float> {mVVBuffer[triangle[0]].xx, mVVBuffer[triangle[0]].yy, mVVBuffer[triangle[0]].zz},
           rti::util::triple<float> {mVVBuffer[triangle[1]].xx, mVVBuffer[triangle[1]].yy, mVVBuffer[triangle[1]].zz},
           rti::util::triple<float> {mVVBuffer[triangle[2]].xx, mVVBuffer[triangle[2]].yy, mVVBuffer[triangle[2]].zz}};
        rti::util::triple<float> centroid = this->centroid(trnglCoords);

        for (auto& vertex : triangle) {
          // Do not reuse the coordinates from above, because here we would depend on the ordering introduced
          // above, which could lead to suddle bugs when changing the code.
          vertex_f4_t& vbv = mVVBuffer[vertex];
          rti::util::triple<float> crds = {vbv.xx, vbv.yy, vbv.zz};
          float tmp = this->distance(rti::util::pair<rti::util::triple<float> > {crds, centroid});
          // RLOG_TRACE
          //   << "Centroid distance: " << tmp << std::endl;
          if (tmp > vbv.radius) {
            vbv.radius = tmp;
            RLOG_TRACE
              << "Setting disc radius of " << vertex << " " << this->prim_to_string(vertex) << std::endl;
          }
        }
      }
      rtcCommitGeometry(mGeometry);
    }
  };
}} // namespace