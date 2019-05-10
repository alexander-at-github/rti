#pragma once

#include <embree3/rtcore.h>
#include <gmsh.h>

#include "rti/logger.hpp"
#include "rti/absc_geometry_from_gmsh.hpp"

namespace rti {
class disc_geometry_from_gmsh : public absc_geometry_from_gmsh {
  public:
    disc_geometry_from_gmsh(RTCDevice& pDevice) {
      init_this(pDevice);
    }
    void invert_surface_normals() override {
      assert(false && "not implemented");
    }
    std::string prim_to_string(unsigned int pPrimID) override {
      assert(false && "not implemented");
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

    mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_DISC_POINT); // using a ray facing disc


    /////////////////////////////////////////////////
    // - vvtags holds the vertices ("nodes in Gmsh")
    // - eetags[someIndex] contains the triangles
    /////////////////////////////////////////////////

    // TODO: Refactor a gmsh_reader into its own class

    std::vector<std::size_t> vvtags;
    std::vector<double> vvxyz;
    std::vector<double> vvuvw;
    gmsh::model::mesh::getNodes(vvtags, vvxyz, vvuvw); // the vvtags seem to be the Nodes in Gmsh.
    // The vertex buffer contains an array of single precision x, y, z floating
    // point coordinates (RTC_FORMAT_FLOAT3 format).
    mNumVertices = vvtags.size();
    // Acquire memory from Embree
    mVVBuffer = (vertex_f4_t*) rtcSetNewGeometryBuffer(
        mGeometry,
        RTC_BUFFER_TYPE_VERTEX,
        0, // slot
        RTC_FORMAT_FLOAT4,
        sizeof(vertex_f4_t),
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
			assert(false && "FIXME: set disc radius");
			mVVBuffer[vvtag].radius == 1;
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
    // Write triangle from Gmsh to Embree
    for (size_t idx = 0; idx < mNumTriangles; ++idx) {
      size_t ntidx = 3 * idx;
			assert(false && "TODO: process trinagles");
      assert(vvtagsminval <= mTTBuffer[idx].v0 && mTTBuffer[idx].v0 <= vvtagsmaxval && "Invalid Vertex");
      assert(vvtagsminval <= mTTBuffer[idx].v1 && mTTBuffer[idx].v1 <= vvtagsmaxval && "Invalid Vertex");
      assert(vvtagsminval <= mTTBuffer[idx].v2 && mTTBuffer[idx].v2 <= vvtagsmaxval && "Invalid Vertex");
      /***************************************
       * Subtracting one to fix the indices. *
       ***************************************/
      //mTTBuffer[idx].v0 -= 1;
      //mTTBuffer[idx].v1 -= 1;
      //mTTBuffer[idx].v2 -= 1;
    }
  }

	/** Pseudo Code **
    * void new_function() {
    *   get-points-from-gmsh
    *   get-triangles-from-gmsh
    *   for each triangle
    *     compute-centroid-for-triangle // also called barycenter
    *     for each corner (point) of the triangle
    *       calculate distance corner<->centroid
    *       set minimum disc radius to distance of corner<->centroid
    * }
	  */
};
}
