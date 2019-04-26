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
#include <string>
#include <xmmintrin.h>

#include <assert.h>
#include <iostream>

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
    inpath = "../resources/box.fine.msh";
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
    BOOST_LOG_SEV(man.lg, blt::debug) << "gmsh::model::mesh::getElemets():";
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
     
    BOOST_LOG_SEV(man.lg, blt::trace) << "eetypes.size()=" << eetypes.size();

    size_t idx = 0;
    for (auto& eetype : eetypes) {
      std::string eename;
      int eedim;
      int eeorder;
      int eenumnodes;
      std::vector<double> eeparamcoord;
      gmsh::model::mesh::getElementProperties(eetype, eename, eedim, eeorder, eenumnodes, eeparamcoord);
      BOOST_LOG_SEV(man.lg, blt::trace) << "eename=" << eename << ", eedim=" << eedim << ", eenumnodes=" << eenumnodes << "; number of elements of that type (eetags[eetype].size()) = " << eetags[idx].size();
      ++idx;
    }

  } // End Debug

  // TODO: Verify the input.

  // Enable huge page support.
  const std::string device_config = "hugepages=1";
  RTCDevice device = rtcNewDevice(device_config.c_str());
  RTCScene scene = rtcNewScene(device);

  // No scene flags for now.
  rtcSetSceneFlags(scene, RTC_SENE_FLAG_NONE);
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

  // Write vertices ("nodes" in Gmsh) from Gmsh to Embree
  for (size_t idx = 0; idx < numvertices; ++idx) {
    size_t npidx = idx * 3;
    assert(npidx <= vvxyz.size() && "npidx out of bounds");
    vvbuffer[idx].xx = vvxyz[npidx];
    vvbuffer[idx].yy = vvxyz[npidx+1];
    vvbuffer[idx].zz = vvxyz[npidx+2];
  }

  std::vector<int> eetypes;
  std::vector<std::vector<std::size_t>> eetags;
  std::vector<std::vector<std::size_t>> nntags;
  // TODO: Can one use getElementsByType() to query for triangles only?
  gmsh::model::mesh::getElements(eetypes, eetags, nntags);
  int selecttriangles = 2;

  size_t numtriangles = eetags[selecttriangles].size();
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
    ttbuffer[idx].v0 = nntags[selecttriangles][3*idx];
    ttbuffer[idx].v1 = nntags[selecttriangles][3*idx+1];
    ttbuffer[idx].v2 = nntags[selecttriangles][3*idx+2];
  }
    
  // Invert surface normals.
  rti::Utils::invertnormals(ttbuffer, numtriangles);

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

  RTCRay ray = rti::Utils::getRay(0,0,0,1,1,1); //TODO
  RTCRayHit rayhit = rti::Utils::getRayHit(ray, RTCHit{}); // TODO
  // TODO: According to gdb there is an error in the following call.
  rtcIntersect1(scene, &context, &rayhit);

  //gmsh::fltk::run();
  gmsh::finalize();
}

