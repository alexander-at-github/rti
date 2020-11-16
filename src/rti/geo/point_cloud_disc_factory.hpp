#pragma once

#include <type_traits>

#include "i_factory.hpp"
#include "point_cloud_disc_geometry.hpp"
//#include "../io/vtp_point_cloud_reader.hpp"
#include "../io/vtp_writer.hpp"
#include "../trace/point_cloud_context.hpp"

namespace rti { namespace geo {
  template<typename numeric_type, typename context_type>
  class point_cloud_disc_factory : public rti::geo::i_factory<numeric_type> {

    // Precodition:
    static_assert(std::is_base_of<rti::trace::absc_context<numeric_type>, context_type>::value,
                  "Precondition");

  public:
    point_cloud_disc_factory(RTCDevice& device,
                             rti::io::i_point_cloud_reader<numeric_type>& reader) :
      mGeometry(device, reader) {}

    point_cloud_disc_factory(RTCDevice& device,
                             std::vector<rti::util::quadruple<numeric_type> > points,
                             std::vector<rti::util::triple<numeric_type> > normals) :
      mGeometry(device, points, normals) {}

    rti::geo::absc_geometry<numeric_type>& get_geometry() override final
    {
      return mGeometry;
    }

    void register_intersect_filter_funs(rti::geo::absc_boundary<numeric_type>& pBoundary) override final
    {
      context_type::register_intersect_filter_funs(mGeometry, pBoundary);
    }

    std::unique_ptr<rti::trace::absc_context<numeric_type> >
    get_new_context(
      unsigned int pGeometryID,
      rti::geo::absc_geometry<numeric_type>& pGeometry,
      rti::reflection::i_reflection_model<numeric_type>& pReflectionModel,
      rti::trace::i_hit_accumulator<numeric_type>& pHitAccumulator,
      unsigned int pBoundaryID,
      rti::geo::absc_boundary<numeric_type>& pBoundary,
      rti::reflection::i_reflection_model<numeric_type>& pBoundaryReflectionModel,
      rti::rng::i_rng& pRng,
      rti::rng::i_rng::i_state& pRngState,
      rti::particle::i_particle<numeric_type>& particle) override final
    {
      auto cntxt = std::make_unique<context_type>
        (pGeometryID,
         // the cast characterizes a precondition to this function
         *dynamic_cast<rti::geo::absc_point_cloud_geometry<numeric_type>*>(&pGeometry),
         pReflectionModel,
         pHitAccumulator, pBoundaryID, pBoundary,
         pBoundaryReflectionModel, pRng, pRngState, particle);
      return std::move(cntxt);
    }

    void write_to_file(rti::trace::i_hit_accumulator<numeric_type>& pHA,
                       std::string pOutfilename,
                       std::vector<rti::util::pair<std::string> > pMetadata) override final
    {
      rti::io::vtp_writer<numeric_type>::write(this->mGeometry, pHA, pOutfilename, pMetadata);
    }

  private:
    rti::geo::point_cloud_disc_geometry<numeric_type> mGeometry;
  };
}}
