#pragma once

#include <embree3/rtcore.h>
#include <gmsh.h>

#include "rti/gmsh_reader.hpp"
#include "rti/logger.hpp"
#include "rti/absc_geometry_from_gmsh.hpp"

namespace rti {
  class sphere_geometry_from_gmsh : public absc_geometry_from_gmsh {
  public:
    sphere_geometry_from_gmsh(RTCDevice& pDevice, gmsh_reader& pGmshReader) :
      absc_geometry_from_gmsh(pDevice) {
      init_this(pDevice, pGmshReader);
    }
    void invert_surface_normals() override {
      BOOST_LOG_SEV(rti::mRLogger, blt::info)
        <<  "Function rti::sphere_geometry_from_gmsh::invert_surface_normals() does not have any functionality";
      //assert(false && "not implemented");
    }
    std::string to_string() override {
      std::stringstream strstream;
      strstream << "[";
      if (mVVBuffer != nullptr) {
        for (size_t idx = 0; idx < mNumVertices; ++idx) {
          strstream << "(" << mVVBuffer[idx].xx << "," << mVVBuffer[idx].yy
                    << "," << mVVBuffer[idx].zz << "," << mVVBuffer[idx].radius << ")";
        }
      }
      strstream << "]";
      return strstream.str();
    }
    std::string prim_to_string(unsigned int pPrimID) override {
      //assert(false && "not implemented");
      std::stringstream strstream;
      strstream << "foo";
      return strstream.str();
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
    void init_this(RTCDevice& pDevice) {
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

      mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_SPHERE_POINT);

      // Read input from gmsh
      std::vector<rti::triple_t<double> > vertices = pGmshReader.get_vertices();

      this->mNumVertices = vertices.size();
      // Acquire memory from Embree
      mVVBuffer = (vertex_f4_t*)
        rtcSetNewGeometryBuffer(
                                mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // slot
                                RTC_FORMAT_FLOAT4,
                                sizeof(vertex_f4_t),
                                mNumVertices);

      // Write vertices to Embree
      BOOST_LOG_SEV(rti::mRLogger, blt::error) << "ERROR: radii of spheres set to some default value"; // TODO
      for (size_t idx = 0; idx < this->mNumVertices; ++idx) {
        auto& triple = vertices[idx];
        mVVBuffer[idx].xx = std::get<0>(triple);
        mVVBuffer[idx].yy = std::get<1>(triple);
        mVVBuffer[idx].zz = std::get<2>(triple);
        mVVBuffer[idx].radius = 0.1; // Set sphere radius to some default value FIXME
      }

      std::vector<rti::triple_t<std::size_t> > triangles = pGmshReader.get_triangles();

      this->mNumTriangles = triangles.size();
      // Set radii of spheres
      BOOST_LOG_SEV(rti::mRLogger, blt::error) << "ERROR: input data not used compute radii of spheres"; // TODO
      for (size_t idx = 0; idx < this->mNumTriangles; ++idx) {
        //assert(false && "TODO: process trinagles");
        /*  Compute centroid of triangle
         *  For each vertex of the triangle
         *    compute distance vertex<->centroid
         *    vertex.radius = max(vertex.radius, distance vertex<->controid)
         */
        //assert(vvtagsminval <= mTTBuffer[idx].v0 && mTTBuffer[idx].v0 <= vvtagsmaxval && "Invalid Vertex");
        //assert(vvtagsminval <= mTTBuffer[idx].v1 && mTTBuffer[idx].v1 <= vvtagsmaxval && "Invalid Vertex");
        //assert(vvtagsminval <= mTTBuffer[idx].v2 && mTTBuffer[idx].v2 <= vvtagsmaxval && "Invalid Vertex");
        /***************************************
         * Subtracting one to fix the indices. *
         ***************************************/
        //mTTBuffer[idx].v0 -= 1;
        //mTTBuffer[idx].v1 -= 1;
        //mTTBuffer[idx].v2 -= 1;
      }
    }
  };
}
