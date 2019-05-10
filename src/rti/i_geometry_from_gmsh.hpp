#pragma once

#include "rti/utils.hpp"

namespace rti {
// Interface ...
class i_geometry_from_gmsh {
  public:
    // ... one should not be able to instanciate it.
    //virtual ~i_geometry_from_gmsh() = default; // = 0; // FIXME
    virtual void invert_surface_normals() = 0;
    virtual std::string prim_to_string(unsigned int) = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
};
} // namespace rti
