#pragma once

#include <cassert>

#include <embree3/rtcore.h>

#include "../io/i_triangle_reader.hpp"
#include "i_geometry.hpp"
#include "../util/logger.hpp"
#include "../util/utils.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class triangle_geometry : public rti::geo::i_geometry<Ty> {
  public:

    triangle_geometry(RTCDevice& pRTCDevice,
                      rti::io::i_triangle_reader<Ty>& pReader) :
      mRTCDevice(pRTCDevice),
      mReader(pReader) {
      init_this(pRTCDevice, pReader);
    }

    Ty get_area(unsigned int primID) override final
    {
      auto tri = mTTBuffer[primID];
      auto v0 = mVVBuffer[tri.v0];
      auto v1 = mVVBuffer[tri.v1];
      auto v2 = mVVBuffer[tri.v2];
      auto triangle = rti::util::triple<rti::util::triple<Ty> >
        {v0.xx, v0.yy, v0.zz,
         v1.xx, v1.yy, v1.zz,
         v2.xx, v2.yy, v2.zz};
      return rti::util::area_of_triangle(triangle);
      // assert(false && "not implemented");
      // return 0;
    }

  private:
    // "RTC_GEOMETRY_TYPE_TRIANGLE: The vertex buffer contains an array of
    // single precision x, y, z floating point coordinates
    // (RTC_FORMAT_FLOAT3 format), and the number of vertices are inferred
    // from the size of that buffer. "
    // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
    struct vertex_f3_t {
      float xx, yy, zz;
    };
    // using vertex_f3_t = rti::util::triple<float>;
    // "RTC_GEOMETRY_TYPE_TRIANGLE: The index buffer contains an array of three
    // 32-bit indices per triangle (RTC_FORMAT_UINT format)"
    // Source: https://embree.github.io/api.html#rtc_geometry_type_triangle
    struct triangle_t {
      uint32_t v0, v1, v2;
    };

    void init_this(RTCDevice& pDevice, rti::io::i_triangle_reader<Ty>& pReader)
    {
      this->mGeometry = rtcNewGeometry(pDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
      auto vertices = pReader.get_points();
      this->mNumVertices = vertices.size();
      this->mVVBuffer = (vertex_f3_t*)
        rtcSetNewGeometryBuffer(this->mGeometry,
                                RTC_BUFFER_TYPE_VERTEX,
                                0, // slot
                                RTC_FORMAT_FLOAT3,
                                sizeof(vertex_f3_t),
                                this->mNumVertices);
      for (size_t idx = 0; idx < this->mNumVertices; ++idx) {
        auto& vv = vertices[idx];
        this->mVVBuffer[idx].xx = vv[0];
        this->mVVBuffer[idx].yy = vv[1];
        this->mVVBuffer[idx].zz = vv[2];
      }
      auto triangles = pReader.get_triangles();
      this->mNumTriangles = triangles.size();
      this->mTTBuffer = (triangle_t*)
        rtcSetNewGeometryBuffer(this->mGeometry,
                                RTC_BUFFER_TYPE_INDEX,
                                0, // slot
                                RTC_FORMAT_UINT3,
                                sizeof(triangle_t),
                                this->mNumTriangles);

      assert (RTC_ERROR_NONE == rtcGetDeviceError(pDevice) && "Error");

      for (size_t idx = 0; idx < this->mNumTriangles; ++idx) {
        auto& tt = triangles[idx];
        this->mTTBuffer[idx].v0 = tt[0];
        this->mTTBuffer[idx].v1 = tt[1];
        this->mTTBuffer[idx].v2 = tt[2];
        // assert(0 <= mTTBuffer[idx][0]); // not necessary; unsigned
        assert(mTTBuffer[idx].v0 < (long long) mNumVertices && "Invalid Vertex");
        assert(mTTBuffer[idx].v1 < (long long) mNumVertices && "Invalid Vertex");
        assert(mTTBuffer[idx].v2 < (long long) mNumVertices && "Invalid Vertex");
      }
      rtcCommitGeometry(this->mGeometry);

      this->mNNBuffer = pReader.get_normals();
    }

  public:

    void print(std::ostream& pOs) override final
    {
      pOs << "(:class " << typeid(this).name();
      if (this->mVVBuffer != nullptr)
        for (size_t idx = 0; idx < this->mNumTriangles; ++idx) {
          pOs << this->prim_to_string(idx);
        }
      pOs << ")";
    }

    RTCDevice& get_rtc_device() override final
    {
      return this->mRTCDevice;
    }

    RTCGeometry& get_rtc_geometry() override final
    {
      return this->mGeometry;
    }

    rti::util::triple<Ty>& get_normal_ref(unsigned int pPrimID) override final
    {
      assert(this->mNNBuffer.size() == this->mNumTriangles &&
             pPrimID <= this->mNumTriangles && "Assumption");
      return reinterpret_cast<rti::util::triple<Ty> > (this->mNNBuffer[pPrimID]);
    }
    
    rti::util::triple<Ty> get_normal(unsigned int pPrimID) override final
    {
      // auto& tri = this->mTTBuffer[pPrimID];
      // auto& v0 = this->mVVBuffer[tri.v0];
      // auto& v1 = this->mVVBuffer[tri.v1];
      // auto& v2 = this->mVVBuffer[tri.v2];
      // auto triangle = rti::util::triple<rti::util::triple<Ty> > {v0.xx, v0.yy, v0.zz,
      //                                                            v1.xx, v1.yy, v1.zz,
      //                                                            v2.xx, v2.yy, v2.zz};
      // return rti::util::compute_normal(triangle);
      //std::cerr << this->mNNBuffer.size() << " " << this->mNumTriangles << " " << pPrimID << std::endl;
      assert(this->mNNBuffer.size() == this->mNumTriangles && pPrimID <= this->mNumTriangles && "Assumption");
      return this->mNNBuffer[pPrimID];
    }

    rti::util::triple<Ty> get_new_origin(RTCRay& pRay, unsigned int primID) override final
    {
      //auto epsilon = 1e-12; // magic number; less than 1e-3 does definitely not work
      auto xx = pRay.org_x + pRay.dir_x * pRay.tfar;
      auto yy = pRay.org_y + pRay.dir_y * pRay.tfar;
      auto zz = pRay.org_z + pRay.dir_z * pRay.tfar;
      // // add a small epsilon from the surface to be sure to be above the surface
      // auto normal = get_normal(pPrimID);
      // xx += normal[0] * epsilon;
      // yy += normal[1] * epsilon;
      // zz += normal[2] * epsilon;
      return {(Ty) xx, (Ty) yy, (Ty) zz};
    }

    std::string get_input_file_path() override final
    {
      return this->mReader.get_input_file_name();
    }

    size_t get_num_primitives() override final
    {
      return this->mNumTriangles;
    }

    std::string prim_to_string(unsigned int pPrimID) override final
    {
      assert(false && "Not implemented");
      return {};
    }

    rti::util::pair<rti::util::triple<Ty> > get_bounding_box() override final
    {
      assert(mVVBuffer != nullptr && "No data");
      if (mVVBuffer == nullptr) // no data in this instance
        return {rti::util::triple<Ty> {0,0,0}, rti::util::triple<Ty> {0,0,0}};
      Ty min = std::numeric_limits<Ty>::lowest();
      Ty max = std::numeric_limits<Ty>::max();
      Ty xmin=max, xmax=min, ymin=max, ymax=min, zmin=max, zmax=min;
      for (size_t idx = 0; idx < mNumVertices; ++idx) {
        xmin = std::min(xmin, mVVBuffer[idx].xx);
        xmax = std::max(xmax, mVVBuffer[idx].xx);
        ymin = std::min(ymin, mVVBuffer[idx].yy);
        ymax = std::max(ymax, mVVBuffer[idx].yy);
        zmin = std::min(zmin, mVVBuffer[idx].zz);
        zmax = std::max(zmax, mVVBuffer[idx].zz);
      }
      return {rti::util::triple<Ty> {xmin, ymin, zmin}, rti::util::triple<Ty> {xmax, ymax, zmax}};
    }

    // Note: this function does not read the data from this->mVVBuffer (the Embree buffer).
    // Instead it reads the data from the input data structure of type
    // rti::io::i_triangle_reader<Ty>.
    std::vector<rti::util::triple<Ty> > get_vertices() const
    {
      return this->mReader.get_points();
    }

    // Note: this function does not read the data from this->mVVBuffer (the Embree buffer).
    // Instead it reads the data from the input data structure of type
    // rti::io::i_triangle_reader<Ty>.
    std::vector<rti::util::triple<size_t> > get_triangles() const
    {
      return this->mReader.get_triangles();
    }

    rti::util::triple<rti::util::triple<Ty> > get_prim(unsigned int pPrimID) const
    {
      auto& tt = this->mTTBuffer[pPrimID];
      return {rti::util::triple<Ty> {this->mVVBuffer}};
    }

  private:
    RTCDevice& mRTCDevice;
    rti::io::i_triangle_reader<Ty> mReader;

    RTCGeometry mGeometry;
    vertex_f3_t* mVVBuffer = nullptr;
    size_t mNumVertices = 0;
    triangle_t* mTTBuffer = nullptr;
    size_t mNumTriangles = 0;
    std::vector<rti::util::triple<Ty> > mNNBuffer;
  };
}}
