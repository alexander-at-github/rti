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

    virtual rti::util::quadruple<numeric_type> get_prim(unsigned int pPrimID) = 0;

    absc_point_cloud_geometry(RTCDevice& pDevice,
                              rti::io::i_point_cloud_reader<numeric_type>& pGReader) :
      mDevice(pDevice),
      mInfilename(pGReader.get_input_file_name()) {}

    absc_point_cloud_geometry(RTCDevice& device) :
      mDevice(device) {}

    void print(std::ostream& pOs) override final
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

    rti::util::pair<rti::util::triple<numeric_type> > get_bounding_box() override final
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

    size_t get_num_primitives() override final
    {
      return this->mNumPoints;
    }

    bool get_relevance(unsigned int primID) override final
    {
      assert(primID <= mNumPoints && "Correctness Assertion");
      auto relevance = (mVVBuffer[primID].xx <= 0.12 && mVVBuffer[primID].zz <= 0.27);
      // std::cout
      //   << "absc_point_cloud_geometry::get_relevanceI() will return "
      //   << (relevance ? "TRUE" : "FALSE") << " "
      //   << "primID == " << primID << " "
      //   << "x == " << mVVBuffer[primID].xx << " "
      //   << "y == " << mVVBuffer[primID].yy << " "
      //   << "z == " << mVVBuffer[primID].zz
      //   << std::endl;
      return relevance;
      // assert(false && "Not implemented");
      // return false;
    }

  protected:
    // Data members
    RTCDevice& mDevice;
    RTCGeometry mGeometry;
    point_4f_t* mVVBuffer = nullptr;
    size_t mNumPoints = 0;
    std::string mInfilename;

  };
}}
