#pragma once

#include <embree3/rtcore.h>
#include <gmsh.h>

#include "rti/absc_geometry_from_gmsh.hpp"
#include "rti/gmsh_reader.hpp"
#include "rti/logger.hpp"
#include "rti/types.hpp"

namespace rti {
  class oriented_disc_geometry_from_gmsh : public absc_geometry_from_gmsh {
  public:
    oriented_disc_geometry_from_gmsh(RTCDevice& pDevice, gmsh_reader& pGmshReader) :
      absc_geometry_from_gmsh(pDevice, pGmshReader) {
      init_this(pDevice, pGmshReader);
    }
    std::string to_string() override {
      std::stringstream strstream;
      strstream << "(:class oriented_disc_geometry_from_gmsh ";
      if (mVVBuffer != nullptr) {
        for (size_t idx = 0; idx < this->mNumVertices; ++idx) {
          strstream << this->prim_to_string(idx);
        }
      }
      strstream << ")";
      return strstream.str();
    }
    std::string prim_to_string(unsigned int pPrimID) override {
      //assert(false && "not implemented");
      std::stringstream strstream;
      strstream << "(" << mVVBuffer[pPrimID].xx << "," << mVVBuffer[pPrimID].yy
                << "," << mVVBuffer[pPrimID].zz << ", radius: " << mVVBuffer[pPrimID].radius
                << ", normal_vector: "
                << "(" << mNNBuffer[pPrimID].xx
                << "," << mNNBuffer[pPrimID].yy
                << "," << mNNBuffer[pPrimID].zz << "))";
      return strstream.str();
    }
  private:
    ////////////////////////
    // Local algebraic types
    ////////////////////////
    // https://embree.github.io/api.html#rtc_geometry_type_point
    struct vertex_f4_t {
      float xx, yy, zz, radius;
    };
    struct normal_vec_f3_t {
      float xx, yy, zz;
    };
    ///////////////
    // Data members
    ///////////////
    // Buffer for vertices / disc centers
    vertex_f4_t* mVVBuffer = nullptr;
    // Buffer for disc normals
    normal_vec_f3_t* mNNBuffer = nullptr;
    size_t mNumTriangles = 0;
    size_t mNumVertices = 0;
    ////////////
    // Functions
    ////////////
    void init_this(RTCDevice& pDevice,gmsh_reader& pGmshReader) {
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

      mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT); // using a disc with surface normal

      // Read input from gmsh
      std::vector<rti::triple<double> > vertices = pGmshReader.get_vertices();

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
        mVVBuffer[idx].xx = triple.frst;
        mVVBuffer[idx].yy = triple.scnd;
        mVVBuffer[idx].zz = triple.thrd;
        mVVBuffer[idx].radius = 0; // Update to sound value later
      }

      std::vector<rti::triple<size_t> > triangles = pGmshReader.get_triangles();

      mNNBuffer = (normal_vec_f3_t*)
        rtcSetNewGeometryBuffer(
                                mGeometry,
                                RTC_BUFFER_TYPE_NORMAL,
                                0, //slot
                                RTC_FORMAT_FLOAT3,
                                sizeof(normal_vec_f3_t),
                                this->mNumVertices);
      // Initialize all disc normals to (0, 0, 0)
      for (size_t idx = 0; idx < this->mNumVertices; ++idx) {
        mNNBuffer[idx].xx = mNNBuffer[idx].yy = mNNBuffer[idx].zz = 0;
      }

      this->mNumTriangles = triangles.size();
      // Set radii and normals of discs
      for (size_t idx = 0; idx < this->mNumTriangles; ++idx) {
        auto& triangle = triangles[idx];

        // Set radius
        rti::triple<rti::triple<float> > trnglCoords
        {{mVVBuffer[triangle.frst].xx, mVVBuffer[triangle.frst].yy, mVVBuffer[triangle.frst].zz},
         {mVVBuffer[triangle.scnd].xx, mVVBuffer[triangle.scnd].yy, mVVBuffer[triangle.scnd].zz},
         {mVVBuffer[triangle.thrd].xx, mVVBuffer[triangle.thrd].yy, mVVBuffer[triangle.thrd].zz}};
        rti::triple<float> centroid = this->centroid(trnglCoords);
        for (auto& vertex : triangle.get_iterable()) {
          // Do not reuse the coordinates from above, because here we would depend on the ordering introduced
          // above, which could lead to suddle bugs when changing the code.
          vertex_f4_t& vbv = mVVBuffer[vertex];
          rti::triple<float> crds = {vbv.xx, vbv.yy, vbv.zz};
          float tmp = this->distance({crds, centroid});
          // BOOST_LOG_SEV(rti::mRLogger, blt::trace)
          //   << "Centroid distance: " << tmp;
          if (tmp > vbv.radius) {
            /********************************
             * Make discs 68% bigger
             ********************************/
            vbv.radius = 1.68f * tmp;
            BOOST_LOG_SEV(rti::mRLogger, blt::trace)
              << "Setting disc radius of " << vertex << " " << this->prim_to_string(vertex);
          }
        }

        // Set normals
        rti::triple<float> normal = this->compute_triangle_normal(triangle);
        // Every disc will get the averaged normal over all adjacent triangles.
        for (auto& vertex : triangle.get_iterable()) {
          mNNBuffer[vertex].xx += normal.frst;
          mNNBuffer[vertex].yy += normal.scnd;
          mNNBuffer[vertex].zz += normal.thrd;
        }
      }
    }

    rti::triple<float> compute_triangle_normal(rti::triple<size_t> ptriangle) {
      size_t p1 = ptriangle.frst;
      size_t p2 = ptriangle.scnd;
      size_t p3 = ptriangle.thrd;
      rti::triple<float> p2m1 {
        mVVBuffer[p2].xx - mVVBuffer[p1].xx,
        mVVBuffer[p2].yy - mVVBuffer[p1].yy,
        mVVBuffer[p2].zz - mVVBuffer[p1].zz};
      rti::triple<float> p3m1 {
        mVVBuffer[p3].xx - mVVBuffer[p1].xx,
        mVVBuffer[p3].yy - mVVBuffer[p1].yy,
        mVVBuffer[p3].zz - mVVBuffer[p1].zz};
      // vector cross product of the two edges p2m1, p3m1
      return {
        (p2m1.scnd * p3m1.thrd) - (p2m1.thrd * p3m1.scnd),
        (p2m1.thrd * p3m1.frst) - (p2m1.frst * p3m1.thrd),
        (p2m1.frst * p3m1.scnd) - (p2m1.scnd * p3m1.frst)};
    }
  };
}
