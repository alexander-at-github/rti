#pragma once

#include "rti/geo/i_abs_geometry.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class i_boundary : public rti::geo::i_abs_geometry<Ty> {
  public:
    virtual ~i_boundary() {}
    virtual std::vector<rti::util::triple<Ty> > get_vertices() const = 0;
    virtual std::vector<rti::util::triple<size_t> > get_triangles() const = 0;
    virtual std::vector<rti::util::triple<Ty> > get_triangle_normals() const = 0;
  };
}} // namespace
