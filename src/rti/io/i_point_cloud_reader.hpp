#pragma once

#include "rti/util/utils.hpp"

namespace rti { namespace io {
  template <typename Ty>
  class i_point_cloud_reader {
  public:
    // pure virtual class
    virtual ~i_point_cloud_reader() {}
    virtual std::vector<rti::util::quadruple<Ty> > get_points() = 0;
    virtual std::vector<rti::util::triple<Ty> > get_normals() = 0;
    virtual std::string get_input_file_name() const = 0;
  };
}} // namespace rti
