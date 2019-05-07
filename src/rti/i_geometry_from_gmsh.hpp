#pragma once

#include "rti/utils.hpp"

namespace rti {
class i_geometry_from_gmsh {
  public:
    // This is an abstract class. One should not be able to instanciate it.
    virtual ~i_geometry_from_gmsh() {};
    //virtual void create_geometry();
    virtual void invert_surface_normals() = 0;
    virtual std::string prim_to_string(unsigned int) = 0;
    RTCGeometry& get_rtc_geometry() {
      return mGeometry;
    }
  protected:
    RTCGeometry mGeometry;
};
} // namespace rti
