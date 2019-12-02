#pragma once

namespace rti { namespace io {
  template<typename Ty>
  class i_triangle_reader {
  public:
    i_triangle_reader(const std::string& pFilename) :
      mInfilename(pFilename) {}

    std::vector<rti::util::triple<Ty> > get_points() const {
      return mPoints;
    }

    std::vector<rti::util::triple<size_t> > get_triangles() const {
      return mTriangles;
    }

    std::string get_input_file_name() {
      return mInfilename;
    }

    virtual std::vector<rti::util::triple<Ty> > get_normals() const {
      // compute normals from triangles
      auto result = std::vector<rti::util::triple<Ty> > {};
      result.reserve(mTriangles.size());
      for (auto const& triangle: mTriangles) {
        auto& v0 = this->mPoints[triangle[0]];
        auto& v1 = this->mPoints[triangle[1]];
        auto& v2 = this->mPoints[triangle[2]];
        // auto triangle = rti::util::triple<rti::util::triple<Ty> > {v0[0], v0[1], v0[2],
        //                                                            v1[0], v1[1], v1[2],
        //                                                            v2[0], v2[1], v2[2]};
        auto trianglepoints = rti::util::triple<rti::util::triple<Ty> > {v0, v1, v2};
        auto normal = rti::util::compute_normal(trianglepoints);
        rti::util::normalize(normal);
        result.push_back(normal);
      }
      result.shrink_to_fit();
      return result;
    }

  protected:
    std::string mInfilename;
    std::vector<rti::util::triple<Ty> > mPoints;
    std::vector<rti::util::triple<size_t> > mTriangles;
  };
}} // namespace
