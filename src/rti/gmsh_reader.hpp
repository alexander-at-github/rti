#pragma once

#include <gmsh.h>

#include "rti/logger.hpp"
#include "rti/types.hpp"

namespace rti {

  class gmsh_reader {
    ////////////////////////////
    // This is a singleton class
    ////////////////////////////
    // For this implementation see suggestions at
    // https://stackoverflow.com/questions/1008019/c-singleton-design-pattern
  public:
    static gmsh_reader& getInstance() {
      // the static keyword ensures that there is only one instance.
      // (At least per translation unit).
      static gmsh_reader instance;
      return instance;
    }

    // Delete default constructor, copy constructor, and copy assignment operator
    //gmsh_reader() = delete;
    gmsh_reader(gmsh_reader const&) = delete;
    gmsh_reader& operator=(gmsh_reader const&) = delete;

  private:
    // Declare constructor private
    gmsh_reader() {
      mMshFilePath = rti::command_line_options::get_instance().
        get_option_value(rti::command_line_options::option_type::MESH_FILE);
      if (mMshFilePath.empty()) {
        // Try default path
        mMshFilePath = "../resources/cylinder/c.msh";
      }
      gmsh::initialize();
      gmsh::option::setNumber("General.Terminal", 1);
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Reading input file " << mMshFilePath;
      gmsh::open(mMshFilePath);
      // BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Will read vertices";
      mVertices = read_vertices();
      // BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Will read triangles";
      mTriangles = read_triangles();
    }
    // Destructor calls gmsh::finalize(); RAII
    ~gmsh_reader() {
      gmsh::finalize();
    }

    ////////////////////////////////////////////////////////
    // Code which is not related to the singleton behaviour.
    ////////////////////////////////////////////////////////
  public:
    std::vector<rti::triple<double> > get_vertices() {
      return this->mVertices;
    }

    std::vector<rti::triple<std::size_t> > get_triangles() {
      return this->mTriangles;
    }

    std::string get_mesh_file_path() {
      return this->mMshFilePath;
    }
  private:
    ///////////////
    // Data Members
    ///////////////
    std::string mMshFilePath;
    std::vector<rti::triple<double> > mVertices;
    std::vector<rti::triple<std::size_t> > mTriangles;
    ////////////
    // Functions
    ////////////
    // std::string getMshFilePath(int argc, char* argv[]) {
    //   std::string optStr{"--msh-file"};
    //   BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Reading command line";
    //   for(int idx = 0; idx < argc; ++idx) {
    //     BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "argv[" << idx << "] == " << argv[idx];
    //   }
    //   for(int idx = 0; idx < argc; ++idx) {
    //     // BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "idx == " << idx << " argv[idx] == " << argv[idx];
    //     if (optStr.compare(argv[idx]) == 0 && idx < (argc-1)) {
    //       std::string filePath(argv[idx+1]);
    //       BOOST_LOG_SEV(rti::mRLogger, blt::debug)
    //         << "Found Mesh file option string: '" << argv[idx] << " " << filePath << "' at index " << idx;
    //       // Mesh file path was found in argv.
    //       return filePath;
    //     }
    //   }
    //   // No mesh file path in argv. Return empty string.
    //   return std::string();
    // }

    std::vector<rti::triple<double> > read_vertices() {
      std::vector<std::size_t> vvtags;
      std::vector<double> vvxyz;
      std::vector<double> vvuvw;
      gmsh::model::mesh::getNodes(vvtags, vvxyz, vvuvw);
      assert(vvxyz.size() == 3 * vvtags.size() && "Vertex data missmatch");

      // BOOST_LOG_SEV(rti::mRLogger, blt::debug)
      //   << "min element is " << *std::min_element(vvtags.begin(), vvtags.end()) << " should be 1";
      // In Gmsh's world node tags start from 1. In our world vertex tags start from 0.
      // We adjust it here. Subtract one from each vertex identifier.
      std::for_each(vvtags.begin(), vvtags.end(), [](auto &tt) {--tt;});
      // assert(*std::min_element(std::begin(vvtags), std::end(vvtags)));
      // BOOST_LOG_SEV(rti::mRLogger, blt::debug)
      //   << "min element is " << *std::min_element(vvtags.begin(), vvtags.end()) << " should be 0";
      assert(*std::min_element(vvtags.begin(), vvtags.end()) == 0 && "Vertex tag assumption not met");
      assert(*std::max_element(vvtags.begin(), vvtags.end()) == vvtags.size()-1
             && "Vertex tag assumption not met");

      //std::vector<double> result(vvxyz.size());
      std::vector<rti::triple<double> > result(vvtags.size());
      // BOOST_LOG_SEV(rti::mRLogger, blt::debug)
      //   << "result vector created";
      for (size_t idx = 0; idx < vvtags.size(); ++idx) {
        size_t xyzidx = 3 * idx;
        //assert(std::is_unsigned<decltype(xyzidx)>::value); // not necessary; unsigned type
        assert(xyzidx < vvxyz.size() && "Error in reading spatial data");
        size_t vvtag = vvtags[idx];
        //assert(std::is_unsigned<decltype(vvtag)>::value); // not necessary; unsigned type
        assert(vvtag < vvtags.size() && "Error in tag/index computation");
        rti::triple<double> rr {vvxyz[xyzidx], vvxyz[xyzidx+1], vvxyz[xyzidx+2]};
        // Would this statement use move semantics without explicit call to std::move()?
        result[vvtag] = std::move(rr);
      }
      return result;
    }

    std::vector<rti::triple<size_t> > read_triangles() {
      std::vector<int> eetypes;
      std::vector<std::vector<std::size_t> > eetags;
      std::vector<std::vector<std::size_t> > nntags;
			// element types are selected by the number of dimensions of that elements.
      // E.g., integer 2 for triangles.
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
      size_t numTriangles = eetags[selectresult].size();
      assert(nntags[selectresult].size() == 3 * numTriangles  && "Size missmatch in triangle data");
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Reading " << eetags[selectresult].size() << " triangles";

      std::vector<std::size_t> selected = nntags[selectresult];
      // Again, like in the get_vertices function, adjust the tags of the vertices to start
      // from 0 instead of 1.
			std::for_each(selected.begin(), selected.end(), [](auto &nn) {--nn;});
      // Some sanity checks
      // Not needed because of unsigned type
      //assert(*std::min_element(selected.begin(), selected.end()) >= 0 && "Vertex tag assumption not met");
      if (this->mVertices.size() > 0) {
        // We can verify this property only if the vertices are set in mVertices.
        assert(*std::max_element(selected.begin(), selected.end()) < mVertices.size() &&
               "Vertex tag assumption not met");
      }

      // Note: we do not consider the element tags (eetags) from Gmsh here. That is, the tags/ids of
      // the triangels may be different than in Gmsh.

      std::vector<rti::triple<size_t> > result(numTriangles);
      for (size_t idx = 0; idx < numTriangles; ++idx) {
        size_t ntidx = 3 * idx;
        //assert(0 <= ntidx); // not needed; unsigned type
        assert(ntidx <= selected.size() && "Index out of bounds");
        rti::triple<size_t> rr {selected[ntidx], selected[ntidx+1], selected[ntidx+2]};
        result[idx] = std::move(rr); // Do we need std::move() for move semantics?
      }
      return result;
    }
  };
}
