#pragma once

#include <array>
#include <memory>
#include <vector>

#include <pmmintrin.h>
#include <xmmintrin.h>

#include <embree3/rtcore.h>


#include "rti/geo/boundary_x_y.hpp"
#include "rti/geo/point_cloud_disc_factory.hpp"
//#include "rti/particle/i_particle.hpp"
#include "rti/particle/i_particle_factory.hpp"
#include "rti/ray/cosine_direction.hpp"
#include "rti/ray/disc_origin_z.hpp"
#include "rti/ray/rectangle_origin_z.hpp"
#include "rti/ray/source.hpp"
#include "rti/trace/point_cloud_context.hpp"
#include "rti/trace/tracer.hpp"
#include "rti/util/utils.hpp"

namespace rti {
  template<typename numeric_type>
  class device final {

  public:

    device()
    {
      init_memory_flags();
    }

    void set_points(std::vector<std::array<numeric_type, 3> > points)
    {
      points = points;
    }

    void set_normals(std::vector<std::array<numeric_type, 3> > normals)
    {
      normals = normals;
    }

    void set_grid_spacing(std::vector<numeric_type> spacing)
    {
      spacing = spacing;
    }

    void set_number_of_rays(size_t numberOfRays)
    {
      numberOfRays = numberOfRays;
    }

    // The register_particle_factory() function takes ownership of the particle.
    // That is, you cannot use the unique_ptr after passing it into this function.
    // Call this function in one of two ways:
    //   (a) auto factory = std::make_unique<concrete::factory_class> (\*constructor-arguments*\)
    //       instance.register_particel_factory(std::move(factory))
    //   (b) instance.register_particle_factory(std::unique_ptr<rti::particle::i_particle_factory>
    //         (new concrete::factory_class(\*constructor-arguments*\)))
    void register_particle_factory
    (std::unique_ptr<rti::particle::i_particle_factory<numeric_type> > particlefactory)
    {
      particlefactory = std::move(particlefactory);
    }

    void run()
    {
      if (particle_factory_is_not_set()) {
        return;
      }
      auto device_config = "hugepages=1";
      auto device = rtcNewDevice(device_config);
      auto pointsandradii = combine_points_with_grid_spacing();
      auto geometryFactory =
        rti::geo::point_cloud_disc_factory<numeric_type, rti::trace::point_cloud_context<numeric_type> >
        (device, pointsandradii, normals);
      auto bdbox = geometryFactory.get_geometry().get_bounding_box();
      bdbox = increase_size_of_bounding_box_by_eps_on_z_axis(bdbox, 0.1);
      //bdbox = increase_size_of_bounding_box_on_x_and_y_axes(bdbox, 8);
      auto origin = create_rectangular_source_from_bounding_box(bdbox);
      // auto origin = create_circular_source_from_bounding_box(bdbox);
      //bdbox = increase_size_of_bounding_box_on_x_and_y_axes(bdbox, 0.5);
      auto boundary = rti::geo::boundary_x_y<numeric_type> {device, bdbox};
      auto direction = rti::ray::cosine_direction<numeric_type>::construct_in_opposite_direction_of_z_axis();
      auto source = rti::ray::source<numeric_type> {origin, direction};
      auto tracer = rti::trace::tracer<numeric_type>
        {geometryFactory, boundary, source, numberOfRays, *(particlefactory)};
      auto traceresult = tracer.run();
      mcestimates = extract_mc_estimates(traceresult);
      hitcnts = extract_hit_cnts(traceresult);
      assert(mcestimates.size() == hitcnts.size() && "Correctness Assumption");
      rtcReleaseDevice(device);
      { // Debug
        //rti::io::vtp_writer<numeric_type>::write(boundary, "bounding-box.vtp");
      }
    }

    std::vector<numeric_type> get_mc_estimates()
    {
      return mcestimates;
    }

    std::vector<size_t> get_hit_cnts()
    {
      return hitcnts;
    }

  private:
    //// Auxiliary functions

    void init_memory_flags()
    {
      _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
      _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    }

    bool particle_factory_is_not_set()
    {
      return (! particlefactory);
    }

    std::vector<size_t>
    extract_hit_cnts(rti::trace::result<numeric_type>& traceresult)
    {
      return traceresult.hitAccumulator->get_cnts();
    }

    std::vector<numeric_type>
    extract_mc_estimates(rti::trace::result<numeric_type>& traceresult)
    {
      auto& hitacc = traceresult.hitAccumulator;
      auto sums = hitacc->get_values();
      auto cntssum = hitacc->get_cnts_sum();
      auto mcestimates = std::vector<numeric_type> {};
      mcestimates.reserve(sums.size());
      for (size_t idx = 0; idx < sums.size(); ++idx) {
        mcestimates.push_back(sums[idx] / cntssum);
      }
      return mcestimates;
    }

    std::vector<rti::util::quadruple<numeric_type> >
    combine_points_with_grid_spacing()
    {
      auto result = std::vector<rti::util::quadruple<numeric_type> > {};
      assert(points.size() == spacing.size() && "Assumption");
      result.reserve(points.size());
      for (size_t idx = 0; idx < points.size(); ++idx) {
        auto tri = points[idx];
        auto sca = spacing[idx];
        result.push_back({tri[0], tri[1], tri[2], sca});
      }
      return result;
    }

    rti::util::pair<rti::util::triple<numeric_type> >
    increase_size_of_bounding_box_by_eps_on_z_axis
    (rti::util::pair<rti::util::triple<numeric_type> > bdbox, double epsT)
    {
      auto result = rti::util::pair<rti::util::triple<numeric_type> >
        {bdbox[0][0], bdbox[0][1], bdbox[0][2], bdbox[1][0], bdbox[1][1], bdbox[1][2]};
      auto eps = (numeric_type) epsT;
      if (result[0][2] > result[1][2]) {
        result[0][2] += eps;
      } else {
        result[1][2] += eps;
      }
      return result;
    }

    rti::util::pair<rti::util::triple<numeric_type> >
    increase_size_of_bounding_box_on_x_and_y_axes
    (rti::util::pair<rti::util::triple<numeric_type> > bdbox, double epsT)
    {
      auto result = rti::util::pair<rti::util::triple<numeric_type> >
        {bdbox[0][0], bdbox[0][1], bdbox[0][2], bdbox[1][0], bdbox[1][1], bdbox[1][2]};
      auto eps = (numeric_type) epsT;
      if (result[0][0] > result[1][0])
        rti::util::swap(result[0][0], result[1][0]);
      if (result[0][1] > result[1][1])
        rti::util::swap(result[0][1], result[1][1]);
      result[0][0] -= eps;
      result[1][0] += eps;
      result[0][1] -= eps;
      result[1][1] += eps;
      return result;
    }

    rti::ray::rectangle_origin_z<numeric_type>
    create_rectangular_source_from_bounding_box(rti::util::pair<rti::util::triple<numeric_type> > bdbox)
    {
      // std::cout
      //   << "################################" << std::endl
      //   << "### Using rectangular source ###" << std::endl
      //   << "################################" << std::endl;
      auto zmax = std::max(bdbox[0][2], bdbox[1][2]);
      auto originC1 = rti::util::pair<numeric_type> {bdbox[0][0], bdbox[0][1]};
      auto originC2 = rti::util::pair<numeric_type> {bdbox[1][0], bdbox[1][1]};
      // ASSUMPTION: the source is on a plain above (positive values) the structure
      // such that z == c for some constant c. (That is in accordance with the silvaco
      // verification instances.)
      return rti::ray::rectangle_origin_z<numeric_type> {zmax, originC1, originC2};
    }

    rti::ray::disc_origin_z<numeric_type>
    create_circular_source_from_bounding_box(rti::util::pair<rti::util::triple<numeric_type> > bdbox)
    {
      // std::cout
      //   << "#############################" << std::endl
      //   << "### Using circular source ###" << std::endl
      //   << "#############################" << std::endl;
      auto zmax = std::max(bdbox[0][2], bdbox[1][2]);
      auto originC1 = rti::util::pair<numeric_type> {bdbox[0][0], bdbox[0][1]};
      auto originC2 = rti::util::pair<numeric_type> {bdbox[1][0], bdbox[1][1]};
      // ASSUMPTION: the source is on a plain above (positive values) the structure
      // such that z == c for some constant c. (That is in accordance with the silvaco
      // verification instances.)
      return rti::ray::disc_origin_z<numeric_type>
        {(originC1[0] + originC2[0]) / 2, (originC1[1] + originC2[1]) / 2, zmax, (originC2[0] - originC1[0])/2};
    }

  private:
    std::vector<rti::util::triple<numeric_type> > points;
    std::vector<rti::util::triple<numeric_type> > normals;
    std::vector<numeric_type> spacing;
    std::vector<numeric_type> mcestimates;
    std::vector<size_t> hitcnts;
    std::unique_ptr<rti::particle::i_particle_factory<numeric_type> > particlefactory;
    size_t numberOfRays = 1024;
  };
}
