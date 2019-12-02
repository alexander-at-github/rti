#pragma once

#include "rti/geo/i_abs_geometry.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  // Interface
  template<typename Ty>
  class i_geometry : public rti::geo::i_abs_geometry<Ty> {
  public:
    // virtual destructor
    virtual ~i_geometry() {}
    //virtual void print(std::ostream& pOs) const = 0;
    //virtual RTCDevice& get_rtc_device() = 0;
    //virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual std::string get_input_file_path() = 0;
    virtual size_t get_num_primitives() const = 0;
    virtual std::string prim_to_string(unsigned int primID) const = 0;
    virtual rti::util::triple<Ty> get_normal(unsigned int primID) const = 0;
    virtual rti::util::pair<rti::util::triple<Ty> > get_bounding_box() const = 0;
    //virtual rti::util::triple<Ty> get_new_origin(RTCRay& pRay, unsigned int primID) const = 0;
    virtual Ty get_sticking_coefficient() const = 0;
    // get_primitive returns the values from the underlying Embree primitive (a disc, a sphere, or something else)
  };
}} // namespace
