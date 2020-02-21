#pragma once

#include "rti/geo/i_boundary.hpp"
#include "rti/geo/i_geometry.hpp"
#include "rti/reflection/i_reflection_model.hpp"
#include "rti/trace/absc_context.hpp"
#include "rti/trace/i_hit_accumulator.hpp"

namespace rti { namespace geo {
  template<typename Ty>
  class i_factory {
  public:
    virtual ~i_factory() {}
    virtual rti::geo::i_geometry<Ty>& get_geometry() = 0;
    virtual void register_intersect_filter_funs(rti::geo::i_boundary<Ty>&) = 0;
    virtual std::unique_ptr<rti::trace::absc_context<Ty> > get_new_context(
      unsigned int pGeometryID,
      rti::geo::i_geometry<Ty>& pGeometry,
      rti::reflection::i_reflection_model<Ty>& pReflectionModel,
      rti::trace::i_hit_accumulator<Ty>& pHitAccumulator,
      unsigned int pBoundaryID,
      rti::geo::i_boundary<Ty>& pBoundary,
      rti::reflection::i_reflection_model<Ty>& pBoundaryReflectionModel,
      rti::rng::i_rng& pRng,
      rti::rng::i_rng::i_state& pRngState) = 0;
    virtual void write_to_file(
      rti::trace::i_hit_accumulator<Ty>& pHA,
      std::string pOutfilename,
      std::vector<rti::util::pair<std::string> > pMetadata) = 0;
  };
}}
