#include <cassert>
#include <vector>

#include <embree3/rtcore.h>

#include "rti/geo/boundary_x_y.hpp"
#include "rti/geo/point_cloud_disc_factory.hpp"
#include "rti/ray/cosine_direction.hpp"
#include "rti/ray/rectangle_origin_z.hpp"
#include "rti/ray/source.hpp"
#include "rti/trace/tracer.hpp"
#include "rti/util/utils.hpp"

#include "device.hpp"

namespace rti {

  //// Auxiliary functions
  template<typename numeric_type>
  std::vector<rti::util::quadruple<numeric_type> >
  combine_points_with_grid_spacing(deviceImpl<numeric_type>& obj)
  {
    auto result = std::vector<rti::util::quadruple<numeric_type> > {};
    assert(obj.points.size() == obj.spacing.size() && "Assumption");
    result.reserve(obj.points.size());
    for (size_t idx = 0; idx < obj.points.size(); ++idx) {
      auto tri = obj.points[idx];
      auto sca = obj.spacing[idx];
      result.push_back({tri[0], tri[1], tri[2], sca});
    }
    return result;
  }

  template<typename numeric_type>
  rti::util::pair<rti::util::triple<numeric_type> >
  increase_size_of_bounding_box_by_eps_on_z_axis
  (rti::util::pair<rti::util::triple<numeric_type> > bdbox, double epsT)
  {
    auto eps = (numeric_type) epsT;
    if (bdbox[0][2] > bdbox[1][2]) {
      bdbox[0][2] += eps;
    } else {
      bdbox[1][2] += eps;
    }
    return bdbox;
  }

  template<typename numeric_type>
  rti::ray::rectangle_origin_z<numeric_type>
  prepare_source_from_bounding_box(rti::util::pair<rti::util::triple<numeric_type> > bdbox)
  {
    auto zmax = std::max(bdbox[0][2], bdbox[1][2]);
    auto originC1 = rti::util::pair<numeric_type> {bdbox[0][0], bdbox[0][1]};
    auto originC2 = rti::util::pair<numeric_type> {bdbox[1][0], bdbox[1][1]};
    // ASSUMPTION: the source is on a plain above (positive values) the structure
    // such that z == c for some constant c. (That is in accordance with the silvaco
    // verification instances.)
    return rti::ray::rectangle_origin_z<numeric_type> {zmax, originC1, originC2};
  }

  template<typename numeric_type>
  std::vector<numeric_type>
  extract_mc_estimates(rti::trace::result<numeric_type>& traceresult)
  {
    auto& hitacc = traceresult.hitAccumulator;
    auto sums = hitacc->get_values();
    auto cnts = hitacc->get_cnts();
    assert (sums.size() == cnts.size() && "Assumption");
    auto mcestimates = std::vector<numeric_type> {};
    mcestimates.reserve(sums.size());
    for (size_t idx = 0; idx < sums.size(); ++idx) {
      mcestimates.push_back(sums[idx] / cnts[idx]);
    }
    return mcestimates;
  }

  //// Class Implementation

  template<typename numeric_type>
  device<numeric_type>::~device() {}

  template<typename numeric_type>
  device<numeric_type>::device(device<numeric_type>&& rhs) {
    assert(false && "TODO");
  }

  template<typename numeric_type>
  device<numeric_type>& device<numeric_type>::operator=(device<numeric_type>&& rhs) {
    assert(false && "TODO");




  }

  template<typename numeric_type>
  device<numeric_type>::device() :
    pimpl(std::make_unique<rti::deviceImpl<numeric_type> > ()) {}

  // template<typename numeric_type>
  // device<numeric_type>::

  template<typename numeric_type>
  void device<numeric_type>::set_points(std::vector<std::array<numeric_type, 3> > points)
  {
    pimpl->points = points;
  }

  template<typename numeric_type>
  void device<numeric_type>::set_normals(std::vector<std::array<numeric_type, 3> > normals)
  {
    pimpl->normals = normals;
  }

  template<typename numeric_type>
  void device<numeric_type>::set_grid_spacing(std::vector<numeric_type> spacing)
  {
    pimpl->spacing = spacing;
  }

  template<typename numeric_type>
  void device<numeric_type>::set_sticking_coefficient(numeric_type stickingCoefficient)
  {
    assert (0 <= stickingCoefficient && stickingCoefficient <= 1 && "Precondition");
    pimpl->stickingCoefficient = stickingCoefficient;
  }

  template<typename numeric_type>
  void device<numeric_type>::set_number_of_rays(size_t numberOfRays)
  {
    pimpl->numberOfRays = numberOfRays;
  }

  template<typename numeric_type>
  void device<numeric_type>::set_particle(std::unique_ptr<rti::i_particle<numeric_type> > particle)
  {
    pimpl->particle = std::move(particle);
  }

  template<typename numeric_type>
  void device<numeric_type>::run()
  {
    auto device_config = "hugepages=1";
    auto device = rtcNewDevice(device_config);
    auto pointsandradii = combine_points_with_grid_spacing(*pimpl);
    auto geometryFactory = rti::geo::point_cloud_disc_factory<numeric_type>
      (device, pointsandradii, pimpl->normals, pimpl->stickingCoefficient);
    auto bdbox = geometryFactory.get_geometry().get_bounding_box();
    bdbox = increase_size_of_bounding_box_by_eps_on_z_axis(bdbox, 0.1);
    auto boundary = rti::geo::boundary_x_y<numeric_type> {device, bdbox};
    auto origin = prepare_source_from_bounding_box(bdbox);
    auto direction = rti::ray::cosine_direction<numeric_type>::construct_in_opposite_direction_of_z_axis();
    auto source = rti::ray::source<numeric_type> {origin, direction};
    auto tracer = rti::trace::tracer<numeric_type> {geometryFactory, boundary, source, pimpl->numberOfRays};
    auto traceresult = tracer.run();
    pimpl->mcestimates = extract_mc_estimates(traceresult);
    rtcReleaseDevice(device);
  }

  template<typename numeric_type>
  std::vector<numeric_type> device<numeric_type>::get_mc_estimates()
  {
    return pimpl->mcestimates;
  }

  // POD struct, hence we take the freedom to make everything public
  template<typename numeric_type> struct deviceImpl
  {
  public:
    std::vector<rti::util::triple<numeric_type> > points;
    std::vector<rti::util::triple<numeric_type> > normals;
    std::vector<numeric_type> spacing;
    std::vector<numeric_type> mcestimates;
    std::unique_ptr<rti::i_particle<numeric_type> > particle;
    numeric_type stickingCoefficient = (numeric_type) 0.5;
    size_t numberOfRays = 1024ul;
  };

}
