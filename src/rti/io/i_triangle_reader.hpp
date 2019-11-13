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
  protected:
    std::string mInfilename;
    std::vector<rti::util::triple<Ty> > mPoints;
    std::vector<rti::util::triple<size_t> > mTriangles;

  };
}} // namespace
