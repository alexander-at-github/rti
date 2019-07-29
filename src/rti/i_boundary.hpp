#pragma once

#include "rti/i_abs_geometry.hpp"

namespace rti {
  template<typename Ty>
  class i_boundary : public i_abs_geometry<Ty> {
  public:
    virtual ~i_boundary() {}
  };
} // namespace rti
