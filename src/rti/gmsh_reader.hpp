#pragma once

#include <gmsh.h>

namespace rti {

  // Alias template for triples
  template<typename T>
  using triple_t = std::tuple<T, T, T>;

  class gmsh_reader {
    ////////////////////////////
    // This is a singleton class
    ////////////////////////////
    // For this implementation see suggestions at
    // https://stackoverflow.com/questions/1008019/c-singleton-design-pattern
  public:
    static gmsh_reader& getInstance(int argc, char* argv[]) {
      // the static keyword ensures that there is only one instance.
      // (At least per translation unit).
      static gmsh_reader instance(argc, argv);
      return instance;
    }

    // Delete default constructor, copy constructor, and copy assignment operator
    gmsh_reader() = delete;
    gmsh_reader(gmsh_reader const&) = delete;
    void operator=(gmsh_reader const&) = delete;

  private:
    // Declare constructor private
    gmsh_reader(int argc, char* argv[]) {
      mMshFilePath = getMshFilePath(argc, argv);
      gmsh::initialize();
      gmsh::option::setNumber("General.Terminal", 1);
      if (mMshFilePath.empty()) {
        // Try default path
        mMshFilePath = "../resources/box.fine.msh";
      }
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Reading input file " << mMshFilePath;
      gmsh::open(mMshFilePath);
      mVertices = read_vertices();
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

    std::vector<triple_t<double> > get_vertices() {
      return this->mVertices;
    }

    std::vector<triple_t<std::size_t> > get_triangles() {
      return this->mTriangles;
    }

  private:
    ///////////////
    // Data Members
    ///////////////
    std::string mMshFilePath;
    std::vector<triple_t<double> > mVertices;
    std::vector<triple_t<std::size_t> > mTriangles;
    ////////////
    // Functions
    ////////////
    std::string getMshFilePath(int argc, char* argv[]) {
      std::string optStr{"--msh-file"};
      for(size_t idx = 0; idx < argc; ++idx) {
        if (optStr.compare(argv[idx]) && idx < (argc-1)) {
          BOOST_LOG_SEV(rti::mRLogger, blt::debug)
            << "Mesh file option string: " << argv[idx] << " " << argv[idx + 1];
          // Mesh file path was found in argv.
          return std::string(argv[idx]);
        }
      }
      // No mesh file path in argv. Return empty string.
      return std::string();
    }

    std::vector<triple_t<double> > read_vertices() {
      std::vector<std::size_t> vvtags;
      std::vector<double> vvxyz;
      std::vector<double> vvuvw;
      gmsh::model::mesh::getNodes(vvtags, vvxyz, vvuvw);
      assert(vvxyz.size() == 3 * vvtags.size() && "Vertex data missmatch");

      // BOOST_LOG_SEV(rti::mRLogger, blt::debug)
      //   << "min element is " << *std::min_element(vvtags.begin(), vvtags.end()) << " should be 0";
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
      std::vector<triple_t<double> > result(vvtags.size());
      for (size_t idx; idx < vvtags.size(); ++idx) {
        size_t xyzidx = 3 * idx;
        assert(0 <= xyzidx && xyzidx < vvxyz.size() && "Error in reading spatial data");
        size_t vvtag = vvtags[idx];
        assert(0 <= vvtag && vvtag < vvtags.size() && "Error in tag/index computation");
        triple_t<double> rr {vvxyz[xyzidx], vvxyz[xyzidx+1], vvxyz[xyzidx+2]};
        result[vvtag] = std::move(rr); // Does this use move semantics without explicit call to std::move()?
      }
      return result;
    }

    std::vector<triple_t<size_t> > read_triangles() {
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
      assert(*std::min_element(selected.begin(), selected.end()) >= 0 && "Vertex tag assumption not met");
      if (this->mVertices.size() > 0) {
        // We can verify this property only if the vertices are set in mVertices.
        assert(*std::max_element(selected.begin(), selected.end()) < mVertices.size() &&
               "Vertex tag assumption not met");
      }

      // Note: we do not consider the element tags (eetags) from Gmsh here. That is, the tags/ids of
      // the triangels may be different than in Gmsh.

      std::vector<triple_t<size_t> > result(numTriangles);
      for (size_t idx = 0; idx < numTriangles; ++idx) {
        size_t ntidx = 3 * idx;
        assert(0 <= ntidx && ntidx <= selected.size() && "Index out of bounds");
        triple_t<size_t> rr {selected[ntidx], selected[ntidx+1], selected[ntidx+2]};
        result[idx] = std::move(rr); // Do we need std::move() for move semantics?
      }
      return result;
    }
  };
}
