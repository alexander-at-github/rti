#pragma once

#include "rti/geo/i_factory.hpp"
#include "rti/geo/point_cloud_disc_geometry.hpp"
#include "rti/io/vtp_point_cloud_reader.hpp"
#include "rti/io/vtp_writer.hpp"
#include "rti/trace/point_cloud_context.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class point_cloud_disc_factory : public rti::geo::i_factory<Ty> {
  public:
    point_cloud_disc_factory(RTCDevice& pDevice,
                             rti::io::i_point_cloud_reader<Ty>& pReader,
                             Ty pStickingC) :
      mGeometry(pDevice, pReader, pStickingC) {}

    rti::geo::i_geometry<Ty>& get_geometry() override final {
      return mGeometry;
    }

    void register_intersect_filter_funs(rti::geo::i_boundary<Ty>& pBoundary) {
      rti::trace::point_cloud_context<Ty>::register_intersect_filter_funs(mGeometry, pBoundary);
    }

    std::unique_ptr<rti::trace::absc_context<Ty> > get_new_context(
      unsigned int pGeometryID,
      rti::geo::i_geometry<Ty>& pGeometry,
      rti::reflection::i_reflection_model<Ty>& pReflectionModel,
      rti::trace::i_hit_accumulator<Ty>& pHitAccumulator,
      unsigned int pBoundaryID,
      rti::geo::i_boundary<Ty>& pBoundary,
      rti::reflection::i_reflection_model<Ty>& pBoundaryReflectionModel,
      rti::rng::i_rng& pRng,
      rti::rng::i_rng::i_state& pRngState) override final {
      auto cntxt = std::make_unique<rti::trace::point_cloud_context<Ty> >
        (pGeometryID,
         // the cast characterizes a precondition to this function
         *dynamic_cast<rti::geo::absc_point_cloud_geometry<Ty>*>(&pGeometry),
         pReflectionModel,
         pHitAccumulator, pBoundaryID, pBoundary,
         pBoundaryReflectionModel, pRng, pRngState);
      return cntxt;
    }

    void write_to_file(rti::trace::i_hit_accumulator<Ty>& pHA,
                       std::string pOutfilename,
                       std::vector<rti::util::pair<std::string> > pMetadata) override final {
      rti::io::vtp_writer<Ty>::write(this->mGeometry, pHA, pOutfilename, pMetadata);
    }

  private:
    rti::geo::point_cloud_disc_geometry<Ty> mGeometry;
  };
}} // namespace
