#include "rti/Utils.h"

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_feature.hpp> // Do we need that header?
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>

#include <embree3/rtcore.h>
#include <gmsh.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

#include <algorithm>
#include <iostream>
#include <string>

#include <assert.h>

namespace blt = boost::log::trivial;
namespace bls = boost::log::sources;

class Manager {
  public:
    // Set logging with levels trace, debug, info, warning, error, and fatal.
    bls::severity_logger<blt::severity_level> lg;
  private:
  Manager() {
    // Program Initialization in constructor.

    // Enable Flush-to-Zero and Denormals-are-Zero for the MXSCR status and
    // control registers for performance reasons.
    // See: https://embree.github.io/api.html#mxcsr-control-and-status-register
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

    // Log to the console
    boost::log::add_console_log(std::cout);
    // optional: set boost log format
    //
    // Set logging severity level
    boost::log::core::get()->set_filter(blt::severity >= blt::trace);
  }
  friend int main(int, char*[]);
};

int main(int argc, char* argv[]) {
  Manager man;

  // Set-up gmsh
  gmsh::initialize();
  // Set Gmsh to print information to the terminal
  gmsh::option::setNumber("General.Terminal", 1);
  std::string inpath;
  if (argc > 1) {
    inpath = argv[1];
  } else {
    // Try default file name
    //inpath = "../resources/test.cyl.msh";
    // Gmsh says 14 nodes and 44 elements for box.msh
    // Gmsh says 314 nodes, 704 elements, 8 points, 72 lines, and 624 triangles for box.fine.msh
    inpath = "/home/alexander/cdr/rti/resources/box.fine.msh";
  }
  std::cout << "Reading input file " << inpath << std::endl;
  // Read gmsh fille from arguments
  // Equivalent to the File->Open menu in the gmsh gui application.
  // This file should be a mesh.
  gmsh::open(inpath);

  //////////////////////////////////////////////////////////////////////////////
  // - vvtags seems to hold the vertices ("nodes in Gmsh") since the number of
  // elements is equal to the number of nodes in the gmsh gui.
  // - eetags[2] seems to contain the triangles, since the size fits the
  // number of triangles shown in the gmsh gui.
  //////////////////////////////////////////////////////////////////////////////

  { // Debug
    BOOST_LOG_SEV(man.lg, blt::debug) << "gmsh::model::mesh::getElements():";
    // Do we need a new model with gmsh::model::add("<name>") first.
    // Do we need gmsh::model::geo::synchronize()?
    std::vector<std::size_t> vvtags;
    std::vector<double> vvcoord;
    std::vector<double> vvpacoord;
    // the vvtags seem to be the Nodes in Gmsh.
    gmsh::model::mesh::getNodes(vvtags, vvcoord, vvpacoord);
    BOOST_LOG_SEV(man.lg, blt::debug) << "vvtags.size()=" << vvtags.size();
    std::vector<int> eetypes;
    std::vector<std::vector<std::size_t>> eetags;
    std::vector<std::vector<std::size_t>> nntags;
    gmsh::model::mesh::getElements(eetypes, eetags, nntags);
    // Later we can query elements of certain dimension, because we actually
    // need points only (dimension 0) and triangles (dimension 2).
    // gmsh::model::mesh::getElements() does not return the parametric
    // coordinates. We need gmsh::model::mesh::getNodes().


    /*********************************************************************************************
     * When querying Gmash with getElements(eetypes, eetags, nntags) then the content of eetypes *
     * is a vector of ints. The node tags vector nntags must be indexed by the same index as     *
     * eetypes. That is, the content of eetypes[idx] is NOT the correct index for nntags. The    *
     * correct index is the variable idx. eetypes[idx] corresponds to nntags[idx].               *
     *********************************************************************************************/

     
    BOOST_LOG_SEV(man.lg, blt::trace) << "eetypes.size()=" << eetypes.size();


    BOOST_LOG_SEV(man.lg, blt::debug) << "##################################################";
    BOOST_LOG_SEV(man.lg, blt::debug) << "### Using loop with 0 <= idx <= eetypes.size() ###";
    BOOST_LOG_SEV(man.lg, blt::debug) << "##################################################";
    for (size_t idx = 0; idx < eetypes.size(); ++idx) {
      int eetype = eetypes[idx];
      std::string eename;
      int eedim;
      int eeorder;
      int eenumnodes;
      std::vector<double> eeparamcoord;
      gmsh::model::mesh::getElementProperties(eetype, eename, eedim, eeorder, eenumnodes, eeparamcoord);
      BOOST_LOG_SEV(man.lg, blt::trace) << "idx=" << idx  << " eename=" << eename << " eetype=" << eetype << " eedim=" << eedim << " eenumnodes=" << eenumnodes << "; number of elements of that type (eetags[idx].size()) = " << eetags[idx].size() << " nntags[idx].size()=" << nntags[idx].size();
    }
    BOOST_LOG_SEV(man.lg, blt::debug) << "###############################################";
    BOOST_LOG_SEV(man.lg, blt::debug) << "### Using loop for (auto& eetype : eetypes) ###";
    BOOST_LOG_SEV(man.lg, blt::debug) << "###############################################";
    for (auto& eetype : eetypes) {
      std::string eename;
      int eedim;
      int eeorder;
      int eenumnodes;
      std::vector<double> eeparamcoord;
      gmsh::model::mesh::getElementProperties(eetype, eename, eedim, eeorder, eenumnodes, eeparamcoord);
      BOOST_LOG_SEV(man.lg, blt::trace) << "eename=" << eename << " eetype=" << eetype << " eedim=" << eedim << " eenumnodes=" << eenumnodes << "; number of elements of that type (eetags[eetype].size()) = " << eetags[eetype].size() << " nntags[eetype].size()=" << nntags[eetype].size();
    }

    {
      BOOST_LOG_SEV(man.lg, blt::debug) << "##################################################";
      BOOST_LOG_SEV(man.lg, blt::debug) << "### Using getElements(\"select triangles\")      ###";
      BOOST_LOG_SEV(man.lg, blt::debug) << "### Using loop with 0 <= idx <= eetypes.size() ###";
      BOOST_LOG_SEV(man.lg, blt::debug) << "##################################################";
      std::vector<int> eetypes;
      std::vector<std::vector<std::size_t>> eetags;
      std::vector<std::vector<std::size_t>> nntags;
      const int selecttriangles = 2; // triangles correspond to two dimensions
      gmsh::model::mesh::getElements(eetypes, eetags, nntags, selecttriangles, -1); // -1: all tags
      for (size_t idx = 0; idx < eetypes.size(); ++idx) {
        int eetype = eetypes[idx];
        std::string eename;
        int eedim;
        int eeorder;
        int eenumnodes;
        std::vector<double> eeparamcoord;
        gmsh::model::mesh::getElementProperties(eetype, eename, eedim, eeorder, eenumnodes, eeparamcoord);
        BOOST_LOG_SEV(man.lg, blt::trace) << "idx=" << idx  << " eename=" << eename << " eetype=" << eetype << " eedim=" << eedim << " eenumnodes=" << eenumnodes << "; number of elements of that type (eetags[idx].size()) = " << eetags[idx].size() << " nntags[idx].size()=" << nntags[idx].size();
      }
    }
  } // End Debug
  BOOST_LOG_SEV(man.lg, blt::debug) << "##############################";
  BOOST_LOG_SEV(man.lg, blt::debug) << "### End of Debug section 1 ###";
  BOOST_LOG_SEV(man.lg, blt::debug) << "##############################";


  // TODO: Verify the input.

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  RTCScene scene = rtcNewScene(device);

  // No scene flags for now.
  rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
  // Selecting higher build quality results in better rendering performance but slower
  // scene commit times. The default build quality for a scene is RTC_BUILD_QUALITY_MEDIUM.
  RTCBuildQuality bbquality = RTC_BUILD_QUALITY_HIGH;

  // TODO: which type of geometry to use?
  //
  // "Points with per vertex radii are supported with sphere, ray-oriented
  // discs, and normal-oriented discs geometric represetntations. Such point
  // geometries are created by passing RTC_GEOMETRY_TYPE_SPHERE_POINT,
  // RTC_GEOMETRY_TYPE_DISC_POINT, or RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT to
  // the rtcNewGeometry function."
  // Source: https://embree.github.io/api.html#rtc_geometry_type_point
  RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);


  // RTC_GEOMETRY_TYPE_TRIANGLE: https://embree.github.io/api.html#rtc_geometry_type_triangle
  // "The index buffer contains an array of three 32-bit indices per triangle (RTC_FORMAT_UINT
  // format) and the number of primitives is inferred from the size of that buffer. The vertex
  // buffer contains an array of single precision x, y, z floating point coordinates
  // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred from the size of that
  // buffer."

  std::vector<std::size_t> vvtags;
  std::vector<double> vvxyz;
  std::vector<double> vvuvw;
  gmsh::model::mesh::getNodes(vvtags, vvxyz, vvuvw); // the vvtags seem to be the Nodes in Gmsh.
  // The vertex buffer contains an array of single precision x, y, z floating
  // point coordinates (RTC_FORMAT_FLOAT3 format).
  int numvertices = vvtags.size();
  // Acquire memory from Embree
  rti::Utils::Vertex* vvbuffer = (rti::Utils::Vertex*) rtcSetNewGeometryBuffer(
      geometry,
      RTC_BUFFER_TYPE_VERTEX,
      0, // slot
      RTC_FORMAT_FLOAT3,
      sizeof(rti::Utils::Vertex),
      numvertices);

  size_t vvtagsmaxval = *std::max_element(std::begin(vvtags), std::end(vvtags));
  // The tag with maximal value in vvtags has a number at most equal to the number
  // of elements in vvtags. That is, tags are numbered consecutively.
  assert(vvtagsmaxval <= vvtags.size() && "Non-consecutive node tags in the vector vvtags");
  // tags in vvtags start from the integer 1.
  size_t vvtagsminval = *std::min_element(std::begin(vvtags), std::end(vvtags));
  assert(vvtagsminval >= 1 && "vvtags contains a tag with value less than 1");

  // TODO: CONTINUE here
  // In Gmsh the node tags are strictly positive. In Embree vertices do not have explicite
  // tags. They are just members of an array. The index of the array starts from zero.
  // In order to fix this discrepancy we adopt the 

  // Write vertices ("nodes" in Gmsh) from Gmsh to Embree
  for (size_t idx = 0; idx < vvtags.size(); ++idx) {
    size_t vcidx = 3 * idx;
    //
    size_t vvtag = vvtags[idx];
    // Gmsh uses tags starting from 1.
    assert(vvtag >= 1);
    // Arrays are indexed starting from 0. Hence, we subtract one.
    vvtag -= 1; // FIXME: This is part of a bug!
    // Note: In Gmsh the tags are strictly positive (starting from 1); at the
    // same time the array of coordinates (vvxyz) still contains useful data
    // starting from index 0. That is vvxyz[0]..vvxyz[2] contains the
    // coordinates for tag 1.
    vvbuffer[vvtag].xx = vvxyz[vcidx];
    vvbuffer[vvtag].yy = vvxyz[vcidx+1];
    vvbuffer[vvtag].zz = vvxyz[vcidx+2];
    //BOOST_LOG_SEV(man.lg, blt::debug) << "(" << vvbuffer[vvtag].xx << " " << vvbuffer[vvtag].yy << " " << vvbuffer[vvtag].zz << ")" << std::endl;
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
  assert(eetags[selectresult].size() * 3 == nntags[selectresult].size() && "Size missmatch in triangle data");

  BOOST_LOG_SEV(man.lg, blt::debug) << "Number of triangles: " << eetags[selectresult].size();

  size_t numtriangles = eetags[selectresult].size();
  rti::Utils::Triangle* ttbuffer = (rti::Utils::Triangle*)rtcSetNewGeometryBuffer(
      geometry,
      RTC_BUFFER_TYPE_INDEX,
      0, // slot
      RTC_FORMAT_UINT3,
      sizeof(rti::Utils::Triangle),
      numtriangles);

  // Check Embree device error
  if (RTC_ERROR_NONE != rtcGetDeviceError(device)) {
    BOOST_LOG_SEV(man.lg, blt::debug) << "Embree device error after rtcSetNewGeometryBuffer()";
  }

  // Write triangle from Gmsh to Embree
  for (size_t idx = 0; idx < numtriangles; ++idx) {
    size_t ntidx = 3 * idx;
    ttbuffer[idx].v0 = nntags[selectresult][ntidx];
    ttbuffer[idx].v1 = nntags[selectresult][ntidx+1];
    ttbuffer[idx].v2 = nntags[selectresult][ntidx+2];
    assert(vvtagsminval <= ttbuffer[idx].v0 && ttbuffer[idx].v0 <= vvtagsmaxval && "Invalid Vertex");
    assert(vvtagsminval <= ttbuffer[idx].v1 && ttbuffer[idx].v1 <= vvtagsmaxval && "Invalid Vertex");
    assert(vvtagsminval <= ttbuffer[idx].v2 && ttbuffer[idx].v2 <= vvtagsmaxval && "Invalid Vertex");
    ttbuffer[idx].v0 -= 1; // FIXME
    ttbuffer[idx].v1 -= 1; // FIXME
    ttbuffer[idx].v2 -= 1; // FIXME
    BOOST_LOG_SEV(man.lg, blt::debug) << "(" << vvbuffer[ttbuffer[idx].v0].xx
                                      << "," << vvbuffer[ttbuffer[idx].v0].yy
                                      << "," << vvbuffer[ttbuffer[idx].v0].zz << ")"
                                      << "(" << vvbuffer[ttbuffer[idx].v1].xx
                                      << "," << vvbuffer[ttbuffer[idx].v1].yy
                                      << "," << vvbuffer[ttbuffer[idx].v1].zz << ")"
                                      << "(" << vvbuffer[ttbuffer[idx].v2].xx
                                      << "," << vvbuffer[ttbuffer[idx].v2].yy
                                      << "," << vvbuffer[ttbuffer[idx].v2].zz << ")";
  }
    
  // Invert surface normals.
  rti::Utils::invertnormals(ttbuffer, numtriangles); //TODO

  rtcCommitGeometry(geometry);
  auto geomID = rtcAttachGeometry(scene, geometry);
  rtcCommitScene(scene);
  // Release geomtery.
  // The geometry object will be destructed when the scene object is
  // destructed, because the geometry object is attached to the scene
  // object
  rtcReleaseGeometry(geometry);

  // Ray queries
  // The Embree tutorial uses one context per ray.
  RTCIntersectContext context;
  rtcInitIntersectContext(&context);

  //RTCRay ray = rti::Utils::constructRay(0,0,1, 1,1,-1); // origin (0,0,1) direction (1,1,-1)
  //RTCRay ray = rti::Utils::constructRay(0,0,1, 1,0,-1); // origin (0,0,1) direction (1,0,-1)
  //RTCRay ray = rti::Utils::constructRay(0,0,0.5, 1,0,0);
  RTCRay ray = rti::Utils::constructRay(0,0,0.5, 1,0.1,0); // origin (0,0,1) direction (1,1,-1)
  RTCRayHit rayhit = rti::Utils::constructRayHit(ray, rti::Utils::constructHit()); // TODO
  BOOST_LOG_SEV(man.lg, blt::trace) << "Before rtcIntersect1()";
  BOOST_LOG_SEV(man.lg, blt::trace) <<  "rayhit.ray.tfar=" << rayhit.ray.tfar;
  // TODO: According to gdb there is an error in the following call.
  rtcIntersect1(scene, &context, &rayhit);
  BOOST_LOG_SEV(man.lg, blt::trace) << "After rtcIntersect1()";
  BOOST_LOG_SEV(man.lg, blt::trace) <<  "rayhit.ray.tnear=" << rayhit.ray.tnear
                                    << " rayhit.ray.tfar=" << rayhit.ray.tfar
                                    << " rayhit.hit.Ng_x=" << rayhit.hit.Ng_x
                                    << " rayhit.hit.Ng_y=" << rayhit.hit.Ng_y
                                    << " rayhit.hit.Ng_z=" << rayhit.hit.Ng_z
                                    << " rayhit.hit.primID=" << rayhit.hit.primID;
  unsigned int primID = rayhit.hit.primID;
  BOOST_LOG_SEV(man.lg, blt::trace) << "hit-primID=" << primID
                                    << " (" << vvbuffer[ttbuffer[primID].v0].xx
                                    << " ," << vvbuffer[ttbuffer[primID].v0].yy
                                    << " ," << vvbuffer[ttbuffer[primID].v0].zz << ")"
                                    <<  "(" << vvbuffer[ttbuffer[primID].v1].xx
                                    << " ," << vvbuffer[ttbuffer[primID].v1].yy
                                    << " ," << vvbuffer[ttbuffer[primID].v1].zz << ")"
                                    <<  "(" << vvbuffer[ttbuffer[primID].v2].xx
                                    << " ," << vvbuffer[ttbuffer[primID].v2].yy
                                    << " ," << vvbuffer[ttbuffer[primID].v2].zz << ")";

  gmsh::fltk::run();
  gmsh::finalize();
}

