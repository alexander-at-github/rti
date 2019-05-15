#pragma once

#include <embree3/rtcore.h>
#include <gmsh.h>

#include "rti/gmsh_reader.hpp"
#include "rti/logger.hpp"
#include "rti/absc_geometry_from_gmsh.hpp"

namespace rti {
  class triangle_geometry_from_gmsh : public absc_geometry_from_gmsh {
  public:
    // triangle_geometry_from_gmsh(RTCDevice& pDevice) {
    //   init_this(pDevice);
    // }
    triangle_geometry_from_gmsh(RTCDevice& pDevice, gmsh_reader& pGmshReader) :
      absc_geometry_from_gmsh(pDevice) {
      init_this(pDevice, pGmshReader);
    }
    void invert_surface_normals() override {
      for (size_t idx = 0; idx < mNumTriangles; ++idx) {
        std::swap<uint32_t>(mTTBuffer[idx].v1, mTTBuffer[idx].v2);
      }
    }
    std::string to_string() override {
      std::stringstream strstream;
      strstream << "(:class triangle_geometry_from_gmsh ";
      for (size_t idxA = 0; idxA < mNumTriangles; ++idxA) {
        strstream << this->prim_to_string(idxA);
        if (idxA < mNumTriangles-1) {
          strstream << " ";
        }
      }
      strstream << ")";
      return strstream.str();
    }
    std::string prim_to_string(unsigned int pPrimID) override {
      std::stringstream strstream;
      strstream << "(" << mVVBuffer[mTTBuffer[pPrimID].v0].xx
                << " " << mVVBuffer[mTTBuffer[pPrimID].v0].yy
                << " " << mVVBuffer[mTTBuffer[pPrimID].v0].zz << ")"
                << "(" << mVVBuffer[mTTBuffer[pPrimID].v1].xx
                << " " << mVVBuffer[mTTBuffer[pPrimID].v1].yy
                << " " << mVVBuffer[mTTBuffer[pPrimID].v1].zz << ")"
                << "(" << mVVBuffer[mTTBuffer[pPrimID].v2].xx
                << " " << mVVBuffer[mTTBuffer[pPrimID].v2].yy
                << " " << mVVBuffer[mTTBuffer[pPrimID].v2].zz << ")";
      return strstream.str();
    }
  private:
    ////////////////////////
    // Local algebraic types
    ////////////////////////
    struct vertex_f3_t {
      float xx, yy, zz; // No padding here!
      // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
      // single precision x, y, z floating point coordinates
      // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred
      // from the size of that buffer. "
      // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
    };
    // Does Embree need aligned memory?
    struct alignas(4) triangle_t {
      uint32_t v0, v1, v2;
      //int v0, v1, v2;
      // "RTC_GEOMETRY_TYPE_TRIANGLE: The index buffer contains an array of three
      // 32-bit indices per triangle (RTC_FORMAT_UINT format)"
      // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
    };
    ///////////////
    // Data members
    ///////////////
    triangle_t* mTTBuffer = nullptr;
    vertex_f3_t* mVVBuffer = nullptr;
    size_t mNumTriangles = 0;
    size_t mNumVertices = 0;
    ////////////
    // Functions
    ////////////
    void init_this(RTCDevice& pDevice, gmsh_reader& pGmshReader) {
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

      // RTC_GEOMETRY_TYPE_TRIANGLE: https://embree.github.io/api.html#rtc_geometry_type_triangle
      // "The index buffer contains an array of three 32-bit indices per triangle (RTC_FORMAT_UINT
      // format) and the number of primitives is inferred from the size of that buffer. The vertex
      // buffer contains an array of single precision x, y, z floating point coordinates
      // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred from the size of that
      // buffer."
      mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_TRIANGLE);

      std::vector<rti::triple_t<double> > vertices = pGmshReader.get_vertices();

      // The vertex buffer contains an array of single precision x, y, z floating
      // point coordinates (RTC_FORMAT_FLOAT3 format).
      mNumVertices = vertices.size();
      // Acquire memory from Embree
      mVVBuffer = (vertex_f3_t*)
        rtcSetNewGeometryBuffer(
                                mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, /* the slot */
                                RTC_FORMAT_FLOAT3,
                                sizeof(vertex_f3_t),
                                mNumVertices);

      // Write vertices to Embree
      for (size_t idx = 0; idx < vertices.size(); ++idx) {
        auto& triple = vertices[idx];
        mVVBuffer[idx].xx = std::get<0>(triple);
        mVVBuffer[idx].yy = std::get<1>(triple);
        mVVBuffer[idx].zz = std::get<2>(triple);
      }

      std::vector<rti::triple_t<size_t> > triangles = pGmshReader.get_triangles();
      mNumTriangles = triangles.size();
      // Acquire memory from Embree
      mTTBuffer = (triangle_t*)
        rtcSetNewGeometryBuffer(
                                mGeometry,
                                RTC_BUFFER_TYPE_INDEX,
                                0, // slot
                                RTC_FORMAT_UINT3,
                                sizeof(triangle_t),
                                mNumTriangles);

      // Check Embree device error
      if (RTC_ERROR_NONE != rtcGetDeviceError(pDevice)) {
        BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Embree device error after rtcSetNewGeometryBuffer()";
      }

      // Write triangle to Embree
      for (size_t idx = 0; idx < triangles.size(); ++idx) {
        auto& triple = triangles[idx];
        mTTBuffer[idx].v0 = std::get<0>(triple);
        mTTBuffer[idx].v1 = std::get<1>(triple);
        mTTBuffer[idx].v2 = std::get<2>(triple);
        assert(0 <= mTTBuffer[idx].v0 && mTTBuffer[idx].v0 < mNumVertices && "Invalid Vertex");
        assert(0 <= mTTBuffer[idx].v1 && mTTBuffer[idx].v1 < mNumVertices && "Invalid Vertex");
        assert(0 <= mTTBuffer[idx].v2 && mTTBuffer[idx].v2 < mNumVertices && "Invalid Vertex");
        //BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "(" << mVVBuffer[mTTBuffer[idx].v0].xx
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v0].yy
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v0].zz << ")"
        //                                  << "(" << mVVBuffer[mTTBuffer[idx].v1].xx
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v1].yy
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v1].zz << ")"
        //                                  << "(" << mVVBuffer[mTTBuffer[idx].v2].xx
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v2].yy
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v2].zz << ")";
      }
    }





















    void init_this__OLD_DELETE(RTCDevice& pDevice) {
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
    
      // RTC_GEOMETRY_TYPE_TRIANGLE: https://embree.github.io/api.html#rtc_geometry_type_triangle
      // "The index buffer contains an array of three 32-bit indices per triangle (RTC_FORMAT_UINT
      // format) and the number of primitives is inferred from the size of that buffer. The vertex
      // buffer contains an array of single precision x, y, z floating point coordinates
      // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred from the size of that
      // buffer."
      mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_TRIANGLE);

      /////////////////////////////////////////////////
      // - vvtags holds the vertices ("nodes in Gmsh")
      // - eetags[someIndex] contains the triangles
      /////////////////////////////////////////////////

      std::vector<std::size_t> vvtags;
      std::vector<double> vvxyz;
      std::vector<double> vvuvw;
      gmsh::model::mesh::getNodes(vvtags, vvxyz, vvuvw); // the vvtags seem to be the Nodes in Gmsh.
      // The vertex buffer contains an array of single precision x, y, z floating
      // point coordinates (RTC_FORMAT_FLOAT3 format).
      mNumVertices = vvtags.size();
      // Acquire memory from Embree
      mVVBuffer = (vertex_f3_t*) rtcSetNewGeometryBuffer(
                                                         mGeometry,
                                                         RTC_BUFFER_TYPE_VERTEX,
                                                         0, // slot
                                                         RTC_FORMAT_FLOAT3,
                                                         sizeof(vertex_f3_t),
                                                         mNumVertices);

      size_t vvtagsmaxval = *std::max_element(std::begin(vvtags), std::end(vvtags));
      // The tag with maximal value in vvtags has a number at most equal to the number
      // of elements in vvtags. That is, tags are numbered consecutively.
      assert(vvtagsmaxval <= vvtags.size() && "Non-consecutive node tags in the vector vvtags");
      // tags in vvtags start from the integer 1.
      size_t vvtagsminval = *std::min_element(std::begin(vvtags), std::end(vvtags));
      assert(vvtagsminval >= 1 && "vvtags contains a tag with value less than 1");

      // In Gmsh the node tags are strictly positive. In Embree vertices do not have explicite
      // tags. They are just members of an array. The index of the array starts from zero.
      // In order to fix this discrepancy we subtract one from each tag (below). When reading
      // the triangles late, one has to perform this subtraction, too!

      // Write vertices ("nodes" in Gmsh) from Gmsh to Embree
      for (size_t idx = 0; idx < vvtags.size(); ++idx) {
        size_t vcidx = 3 * idx;
        //
        size_t vvtag = vvtags[idx];
        // Gmsh uses tags starting from 1.
        assert(vvtag >= 1);
        // Arrays are indexed starting from 0. Hence, we subtract one.
        /*********************************************
         * This can easily cause a bug further down. *
         *********************************************/
        vvtag -= 1; 
        // vvxyz[0]..vvxyz[2] contains the coordinates for gmsh-tag 1.
        mVVBuffer[vvtag].xx = vvxyz[vcidx];
        mVVBuffer[vvtag].yy = vvxyz[vcidx+1];
        mVVBuffer[vvtag].zz = vvxyz[vcidx+2];
      }

      std::vector<int> eetypes;
      std::vector<std::vector<std::size_t>> eetags;
      std::vector<std::vector<std::size_t>> nntags;
      // element types are selected by the number of dimensions of that elements. E.g., integer 2 for
      // triangles.
      int selecttriangles = 2; // dimensions
      gmsh::model::mesh::getElements(
                                     eetypes,
                                     eetags,
                                     nntags,
                                     selecttriangles, // dimension
                                     -1); // select all elements with respect to their tag
      // When calling gmsh::getElements() with a dimension argument, then the vectors
      // eetypes, eetags and nntags are of size 1.
      assert(eetypes.size() == 1 && eetags.size() == 1 && nntags.size() == 1 && "Assumptions not met");
      int selectresult = 0;
      // Sanity check
      assert(eetags[selectresult].size() * 3 == nntags[selectresult].size() &&
             "Size missmatch in triangle data");

      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Reading "
                                               << eetags[selectresult].size() << " triangles";

      mNumTriangles = eetags[selectresult].size();
      mTTBuffer = (triangle_t*)rtcSetNewGeometryBuffer(
                                                       mGeometry,
                                                       RTC_BUFFER_TYPE_INDEX,
                                                       0, // slot
                                                       RTC_FORMAT_UINT3,
                                                       sizeof(triangle_t),
                                                       mNumTriangles);

      // Check Embree device error
      if (RTC_ERROR_NONE != rtcGetDeviceError(pDevice)) {
        BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Embree device error after rtcSetNewGeometryBuffer()";
      }

      // Write triangle from Gmsh to Embree
      for (size_t idx = 0; idx < mNumTriangles; ++idx) {
        size_t ntidx = 3 * idx;
        mTTBuffer[idx].v0 = nntags[selectresult][ntidx];
        mTTBuffer[idx].v1 = nntags[selectresult][ntidx+1];
        mTTBuffer[idx].v2 = nntags[selectresult][ntidx+2];
        assert(vvtagsminval <= mTTBuffer[idx].v0 && mTTBuffer[idx].v0 <= vvtagsmaxval && "Invalid Vertex");
        assert(vvtagsminval <= mTTBuffer[idx].v1 && mTTBuffer[idx].v1 <= vvtagsmaxval && "Invalid Vertex");
        assert(vvtagsminval <= mTTBuffer[idx].v2 && mTTBuffer[idx].v2 <= vvtagsmaxval && "Invalid Vertex");
        /***************************************
         * Subtracting one to fix the indices. *
         ***************************************/
        mTTBuffer[idx].v0 -= 1;
        mTTBuffer[idx].v1 -= 1;
        mTTBuffer[idx].v2 -= 1;
        //BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "(" << mVVBuffer[mTTBuffer[idx].v0].xx
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v0].yy
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v0].zz << ")"
        //                                  << "(" << mVVBuffer[mTTBuffer[idx].v1].xx
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v1].yy
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v1].zz << ")"
        //                                  << "(" << mVVBuffer[mTTBuffer[idx].v2].xx
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v2].yy
        //                                  << "," << mVVBuffer[mTTBuffer[idx].v2].zz << ")";
      }
    }
  };
}
