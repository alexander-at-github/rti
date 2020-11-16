#pragma once

#include "absc_boundary.hpp"
#include "absc_geometry.hpp"
#include "../particle/i_particle.hpp"
#include "../reflection/i_reflection_model.hpp"
#include "../trace/absc_context.hpp"
#include "../trace/i_hit_accumulator.hpp"

namespace rti { namespace geo {
  template<typename numeric_type>
  class i_factory {
  public:
    virtual ~i_factory() {}
    virtual rti::geo::absc_geometry<numeric_type>& get_geometry() = 0;
    virtual void register_intersect_filter_funs(rti::geo::absc_boundary<numeric_type>&) = 0;
    virtual std::unique_ptr<rti::trace::absc_context<numeric_type> > get_new_context(
      unsigned int pGeometryID,
      rti::geo::absc_geometry<numeric_type>& pGeometry,
      rti::reflection::i_reflection_model<numeric_type>& pReflectionModel,
      rti::trace::i_hit_accumulator<numeric_type>& pHitAccumulator,
      unsigned int pBoundaryID,
      rti::geo::absc_boundary<numeric_type>& pBoundary,
      rti::reflection::i_reflection_model<numeric_type>& pBoundaryReflectionModel,
      rti::rng::i_rng& pRng,
      rti::rng::i_rng::i_state& pRngState,
      rti::particle::i_particle<numeric_type>& particle) = 0;
    virtual void write_to_file(
      rti::trace::i_hit_accumulator<numeric_type>& pHA,
      std::string pOutfilename,
      std::vector<rti::util::pair<std::string> > pMetadata) = 0;
  };
}}
