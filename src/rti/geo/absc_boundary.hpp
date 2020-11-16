#pragma once

#include "meta_geometry.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class absc_boundary : public rti::geo::meta_geometry<Ty> {
  public:
    virtual ~absc_boundary() {}
    virtual std::vector<rti::util::triple<Ty> > get_vertices() = 0;
    virtual std::vector<rti::util::triple<size_t> > get_triangles() = 0;
    virtual std::vector<rti::util::triple<Ty> > get_triangle_normals() = 0;
  };
}} // namespace
