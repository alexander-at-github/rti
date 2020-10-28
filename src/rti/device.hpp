#pragma once

#include <array>
#include <memory>
#include <vector>

#include <xmmintrin.h>
#include <pmmintrin.h>

#include <embree3/rtcore.h>


#include "geo/boundary_x_y.hpp"
#include "geo/point_cloud_disc_factory.hpp"
#include "io/vtp_writer.hpp"
//#include "particle/i_particle.hpp"
#include "particle/i_particle_factory.hpp"
#include "ray/cosine_direction.hpp"
#include "ray/disc_origin_z.hpp"
#include "ray/rectangle_origin_z.hpp"
#include "ray/source.hpp"
#include "trace/point_cloud_context.hpp"
#include "trace/tracer.hpp"
#include "util/utils.hpp"

namespace rti {
  template<typename numeric_type>
  class device final {

  public:

    device()
    {
      init_memory_flags();
    }

    void set_points(std::vector<std::array<numeric_type, 3> > points_)
    {
      points = points_;
    }

    void set_normals(std::vector<std::array<numeric_type, 3> > normals_)
    {
      normals = normals_;
    }

    void set_grid_spacing(std::vector<numeric_type> spacing_)
    {
      spacing = spacing_;
    }

    void set_number_of_rays(size_t numberOfRays_)
    {
      numberOfRays = numberOfRays_;
    }

    // The register_particle_factory() function takes ownership of the particle.
    // That is, you cannot use the unique_ptr after passing it into this function.
    // Call this function in one of two ways:
    //   (a) auto factory = std::make_unique<concrete::factory_class> (\*constructor-arguments*\)
    //       instance.register_particel_factory(std::move(factory))
    //   (b) instance.register_particle_factory(std::unique_ptr<rti::particle::i_particle_factory>
    //         (new concrete::factory_class(\*constructor-arguments*\)))
    void register_particle_factory
    (std::unique_ptr<rti::particle::i_particle_factory<numeric_type> > pf)
    {
      particlefactory = std::move(pf);
    }

    void run()
    {
      if (particle_factory_is_not_set()) {
        return;
      }
      auto device_config = "hugepages=1";
      auto device = rtcNewDevice(device_config);
      auto pointsandradii = combine_points_with_grid_spacing_and_compute_max_disc_radius();
      auto geometryFactory =
        rti::geo::point_cloud_disc_factory<numeric_type, rti::trace::point_cloud_context<numeric_type> >
        (device, pointsandradii, normals);
      auto bdbox = geometryFactory.get_geometry().get_bounding_box();
      // increase by 5% of the z-achsis
      // auto bdboxEps =  0.01 * std::abs(bdbox[1][2] - bdbox[0][2]);
      // auto bdboxEps = maxDscRad * 2;
      auto bdboxEps = maxDscRad * 1;
      std::cout << "[Alex] Using bdboxEpx == " << bdboxEps << std::endl;
      bdbox = increase_size_of_bounding_box_by_eps_on_z_axis(bdbox, bdboxEps);
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
      // mcestimates = extract_mc_estimates(traceresult);
      // normalize_mc_estimates();
      mcestimates = extract_mc_estimates_normalized(traceresult);
      hitcnts = extract_hit_cnts(traceresult);
      assert(mcestimates.size() == hitcnts.size() && "Correctness Assumption");
      rtcReleaseDevice(device);
      // { // Debug
      //   auto path = "/home/alexanders/vtk/outputs/bounding-box.vtp";
      //   std::cout << "Writing bounding box to " << path << std::endl;
      //   rti::io::vtp_writer<numeric_type>::write(boundary, path);
      // }
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

    void normalize_mc_estimates()
    {
      std::cout << "[Alex] normalizing_mc_estimates()" << std::endl;
      auto sum = 0.0;
      auto max = 0.0;
      for (auto const& ee : mcestimates) {
        sum += ee;
        if (ee > max) {
          max = ee;
        }
      }
      for (auto& ee : mcestimates) {
        // ee /= sum;
        ee /= max;
      }
    }
    
    void init_memory_flags()
    {
      #pragma omp parallel
      {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
      }
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
    extract_mc_estimates_normalized(rti::trace::result<numeric_type>& traceresult)
    {
      auto& hitacc = traceresult.hitAccumulator;
      auto values = hitacc->get_values();
      auto areas = hitacc->get_exposed_areas();

      auto mcestimates = std::vector<numeric_type> {};
      mcestimates.reserve(values.size());
      
      auto maxv = 0.0;
      assert (values.size() == areas.size() && "Correctness Assertion");
      // Account for area and find max value
      for (size_t idx = 0; idx < values.size(); ++idx) {
        auto vv = values[idx] / areas[idx];
        mcestimates.push_back(vv);
        if (maxv < vv) {
          maxv = vv;
        }
      }
      std::cout << "[Alex] normalizing_mc_estimates()" << std::endl;
      for (auto& vv : mcestimates) {
        vv /= maxv;
      }
      return mcestimates;
    }
    
    std::vector<numeric_type>
    extract_mc_estimates(rti::trace::result<numeric_type>& traceresult)
    {
      auto& hitacc = traceresult.hitAccumulator;
      auto values = hitacc->get_values();
      auto areas = hitacc->get_exposed_areas();
      auto mcestimates = std::vector<numeric_type> {};
      mcestimates.reserve(values.size());
      assert (values.size() == areas.size() && "Correctness Assertion");
      // Account for area and find max value
      for (size_t idx = 0; idx < values.size(); ++idx) {
        auto vv = values[idx] / areas[idx];
        mcestimates.push_back(vv);
      }
      return mcestimates;
    }

    std::vector<rti::util::quadruple<numeric_type> >
    combine_points_with_grid_spacing_and_compute_max_disc_radius()
    {
      auto result = std::vector<rti::util::quadruple<numeric_type> > {};
      assert(points.size() == spacing.size() && "Assumption");
      maxDscRad = 0.0;
      result.reserve(points.size());
      for (size_t idx = 0; idx < points.size(); ++idx) {
        auto tri = points[idx];
        auto sca = spacing[idx];
        result.push_back({tri[0], tri[1], tri[2], sca});
        if (maxDscRad < sca) {
          maxDscRad = sca;
        }
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
    float maxDscRad = 0.0;
  };
}
