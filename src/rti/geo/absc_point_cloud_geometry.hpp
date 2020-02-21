#pragma once

#include <limits>

#include <boost/core/demangle.hpp>

#include "rti/geo/i_geometry.hpp"
#include "rti/io/i_point_cloud_reader.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace geo {

    // "RTC_GEOMETRY_TYPE_POINT:
    // The vertex buffer stores each control vertex in the form of a single
    // precision position and radius stored in (x, y, z, r) order in memory
    // (RTC_FORMAT_FLOAT4 format). The number of vertices is inferred from the
    // size of this buffer. Similarly, the normal buffer stores a single
    // precision normal per control vertex (x, y, z order and RTC_FORMAT_FLOAT3 format)."
    // Source: https://embree.github.io/api.html#rtc_geometry_type_point
    struct point_4f_t {
      float xx, yy, zz, radius;
    };

  template<typename numeric_type>
  class absc_point_cloud_geometry : public rti::geo::i_geometry<numeric_type> {
  public:
    // abstract class
    virtual ~absc_point_cloud_geometry() {}
    // inherits also some virtual function declarations from i_geometry

    virtual rti::util::quadruple<numeric_type> get_prim(unsigned int pPrimID) const = 0;

    absc_point_cloud_geometry(RTCDevice& pDevice,
                              rti::io::i_point_cloud_reader<numeric_type>& pGReader,
                              numeric_type pStickingC) :
      mDevice(pDevice),
      mStickingC(pStickingC),
      mInfilename(pGReader.get_input_file_name()) {}

    absc_point_cloud_geometry(RTCDevice& device, numeric_type stickingC) :
      mDevice(device),
      mStickingC(stickingC) {}

    void print(std::ostream& pOs) const override final
    {
      pOs << "(:class " << boost::core::demangle(typeid(this).name());
      if (mVVBuffer != nullptr)
        for (size_t idx = 0; idx < this->mNumPoints; ++idx)
          pOs << this->prim_to_string(idx);
      pOs << ")";
    }

    RTCDevice& get_rtc_device() override final
    {
      return mDevice;
    }

    RTCGeometry& get_rtc_geometry() override final
    {
      return mGeometry;
    }

    std::string get_input_file_path() override final
    {
      return mInfilename;
    }

    rti::util::pair<rti::util::triple<numeric_type> > get_bounding_box() const override final
    {
      assert(mVVBuffer != nullptr && "No data");
      if (mVVBuffer == nullptr) // no data in this instance
        return {rti::util::triple<numeric_type> {0,0,0}, rti::util::triple<numeric_type> {0,0,0}};
      numeric_type min = std::numeric_limits<numeric_type>::lowest();
      numeric_type max = std::numeric_limits<numeric_type>::max();
      numeric_type xmin=max, xmax=min, ymin=max, ymax=min, zmin=max, zmax=min;
      for (size_t idx = 0; idx < mNumPoints; ++idx) {
        xmin = std::min(xmin, (numeric_type) mVVBuffer[idx].xx);
        xmax = std::max(xmax, (numeric_type) mVVBuffer[idx].xx);
        ymin = std::min(ymin, (numeric_type) mVVBuffer[idx].yy);
        ymax = std::max(ymax, (numeric_type) mVVBuffer[idx].yy);
        zmin = std::min(zmin, (numeric_type) mVVBuffer[idx].zz);
        zmax = std::max(zmax, (numeric_type) mVVBuffer[idx].zz);
      }
      return {rti::util::triple<numeric_type> {xmin, ymin, zmin},
              rti::util::triple<numeric_type> {xmax, ymax, zmax}};
    }

    size_t get_num_primitives() const override final
    {
      return this->mNumPoints;
    }

    numeric_type get_sticking_coefficient() const override final
    {
      return mStickingC;
    }

  protected:
    // Data members
    RTCDevice& mDevice;
    numeric_type mStickingC = 1; // initialize to some value
    RTCGeometry mGeometry;
    point_4f_t* mVVBuffer = nullptr;
    size_t mNumPoints = 0;
    std::string mInfilename;

  };
}}
