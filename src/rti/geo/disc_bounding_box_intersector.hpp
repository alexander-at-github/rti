#pragma once

#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  class disc_bounding_box_intersector {
    using numeric_type = float;
  public:
    disc_bounding_box_intersector(rti::util::pair<rti::util::pair<numeric_type> > boundingbox) :
      boundingbox(boundingbox) {}

    numeric_type area_inside(rti::util::quadruple<numeric_type> disc, rti::util::triple<numeric_type> normal)
    {
      return rti::util::pi();
    }

  private:
    rti::util::pair<rti::util::pair<numeric_type> > boundingbox;
  };
}}
