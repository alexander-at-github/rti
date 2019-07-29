#pragma once

#include "rti/geo/i_abs_geometry.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class i_boundary : public rti::geo::i_abs_geometry<Ty> {
  public:
    virtual ~i_boundary() {}
  };
}} // namespace
