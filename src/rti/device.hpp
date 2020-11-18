#pragma once

#include <array>
#include <memory>
#include <vector>

#include <xmmintrin.h>
#include <pmmintrin.h>

#include <embree3/rtcore.h>


#include "geo/boundary_x_y.hpp"
#include "geo/bound_condition.hpp"
#include "geo/point_cloud_disc_factory.hpp"
#include "io/vtp_writer.hpp"
#include "particle/i_particle.hpp"
#include "ray/i_direction.hpp"
#include "ray/cosine_direction_z.hpp"
#include "ray/power_cosine_direction_z.hpp"
#include "ray/disc_origin_z.hpp"
#include "ray/rectangle_origin_z.hpp"
#include "ray/source.hpp"
#include "reflection/i_reflection.hpp"
#include "reflection/diffuse.hpp"
#include "trace/point_cloud_context.hpp"
#include "trace/tracer.hpp"
#include "util/utils.hpp"

namespace rti {
  using ::rti::geo::bound_condition;
  using ::rti::particle::i_particle;

  template<typename numeric_type, typename particle_type>
  class device final {

    static_assert(std::is_base_of<i_particle<numeric_type>, particle_type>::value,
                  "Precondition");

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

    void set_number_of_rays(size_t numofrays_)
    {
      numofrays = numofrays_;
    }

    void set_x(bound_condition cond)
    {
      xCond = cond;
    }

    void set_y(bound_condition cond)
    {
      yCond = cond;
    }

    void set_cosine_source()
    {
      srcDirection = std::make_unique<ray::cosine_direction_z<numeric_type> > ();
    }

    void set_power_cosine_source(numeric_type exponent)
    {
      srcDirection = std::make_unique<ray::power_cosine_direction_z<numeric_type> > (exponent);
    }

    void set(reflection::i_reflection<numeric_type>& reflection_)
    {
      reflection = reflection_;
    }

    void run()
    {
      auto device_config = "hugepages=1";
      auto device = rtcNewDevice(device_config);
      auto pointsandradii = combine_points_with_grid_spacing_and_compute_max_disc_radius();
      auto geometry = geo::point_cloud_disc_geometry<numeric_type>{device, pointsandradii, normals};
      auto bdbox = geometry.get_bounding_box();
      auto bdboxEps = maxDscRad;
      bdbox = increase_size_of_bounding_box_by_eps_on_z_axis(bdbox, bdboxEps);
      //bdbox = increase_size_of_bounding_box_on_x_and_y_axes(bdbox, 8);
      auto origin = create_rectangular_source_from_bounding_box(bdbox);
      auto boundary = geo::boundary_x_y<numeric_type> {device, bdbox, xCond, yCond};
      auto source = ray::source<numeric_type> {origin, *srcDirection};
      auto tracer = trace::tracer<numeric_type, particle_type>
        {geometry, boundary, source, reflection, numofrays};
      auto traceresult = tracer.run();
      mcestimates = extract_mc_estimates_normalized(traceresult);
      hitcnts = extract_hit_cnts(traceresult);
      assert(mcestimates.size() == hitcnts.size() && "Correctness Assumption");
      rtcReleaseDevice(device);
      // { // Debug
      //   auto path = "/home/alexanders/vtk/outputs/bounding-box.vtp";
      //   std::cout << "Writing bounding box to " << path << std::endl;
      //   io::vtp_writer<numeric_type>::write(boundary, path);
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
      // std::cout << "[Alex] normalizing_mc_estimates()" << std::endl;
      auto sum = 0.0;
      auto max = 0.0;
      for (auto const& ee : mcestimates) {
        sum += ee;
        if (ee > max) {
          max = ee;
        }
      }
      for (auto& ee : mcestimates) {
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

    std::vector<size_t>
    extract_hit_cnts(trace::result<numeric_type>& traceresult)
    {
      return traceresult.hitAccumulator->get_cnts();
    }
    
    std::vector<numeric_type>
    extract_mc_estimates_normalized(trace::result<numeric_type>& traceresult)
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
      // std::cout << "[Alex] normalizing_mc_estimates()" << std::endl;
      for (auto& vv : mcestimates) {
        vv /= maxv;
      }
      return mcestimates;
    }
    
    std::vector<numeric_type>
    extract_mc_estimates(trace::result<numeric_type>& traceresult)
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

    std::vector<util::quadruple<numeric_type> >
    combine_points_with_grid_spacing_and_compute_max_disc_radius()
    {
      auto result = std::vector<util::quadruple<numeric_type> > {};
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

    util::pair<util::triple<numeric_type> >
    increase_size_of_bounding_box_by_eps_on_z_axis
    (util::pair<util::triple<numeric_type> > bdbox, double epsT)
    {
      auto result = util::pair<util::triple<numeric_type> >
        {bdbox[0][0], bdbox[0][1], bdbox[0][2], bdbox[1][0], bdbox[1][1], bdbox[1][2]};
      auto eps = (numeric_type) epsT;
      if (result[0][2] > result[1][2]) {
        result[0][2] += eps;
      } else {
        result[1][2] += eps;
      }
      return result;
    }

    util::pair<util::triple<numeric_type> >
    increase_size_of_bounding_box_on_x_and_y_axes
    (util::pair<util::triple<numeric_type> > bdbox, double epsT)
    {
      auto result = util::pair<util::triple<numeric_type> >
        {bdbox[0][0], bdbox[0][1], bdbox[0][2], bdbox[1][0], bdbox[1][1], bdbox[1][2]};
      auto eps = (numeric_type) epsT;
      if (result[0][0] > result[1][0])
        util::swap(result[0][0], result[1][0]);
      if (result[0][1] > result[1][1])
        util::swap(result[0][1], result[1][1]);
      result[0][0] -= eps;
      result[1][0] += eps;
      result[0][1] -= eps;
      result[1][1] += eps;
      return result;
    }

    ray::rectangle_origin_z<numeric_type>
    create_rectangular_source_from_bounding_box(util::pair<util::triple<numeric_type> > bdbox)
    {
      auto zmax = std::max(bdbox[0][2], bdbox[1][2]);
      auto originC1 = util::pair<numeric_type> {bdbox[0][0], bdbox[0][1]};
      auto originC2 = util::pair<numeric_type> {bdbox[1][0], bdbox[1][1]};
      // ASSUMPTION: the source is on a plain above (positive values) the structure
      // such that z == c for some constant c. (That is in accordance with the silvaco
      // verification instances.)
      return ray::rectangle_origin_z<numeric_type> {zmax, originC1, originC2};
    }

    ray::disc_origin_z<numeric_type>
    create_circular_source_from_bounding_box(util::pair<util::triple<numeric_type> > bdbox)
    {
      auto zmax = std::max(bdbox[0][2], bdbox[1][2]);
      auto originC1 = util::pair<numeric_type> {bdbox[0][0], bdbox[0][1]};
      auto originC2 = util::pair<numeric_type> {bdbox[1][0], bdbox[1][1]};
      // ASSUMPTION: the source is on a plain above (positive values) the structure
      // such that z == c for some constant c. (That is in accordance with the silvaco
      // verification instances.)
      return ray::disc_origin_z<numeric_type>
        {(originC1[0] + originC2[0]) / 2, (originC1[1] + originC2[1]) / 2, zmax, (originC2[0] - originC1[0])/2};
    }

  private:
    std::vector<util::triple<numeric_type> > points;
    std::vector<util::triple<numeric_type> > normals;
    std::vector<numeric_type> spacing;
    std::vector<numeric_type> mcestimates;
    std::vector<size_t> hitcnts;

    size_t numofrays = 1024;
    numeric_type maxDscRad = 0.0;

    bound_condition xCond = geo::bound_condition::REFLECTIVE;
    bound_condition yCond = geo::bound_condition::REFLECTIVE;

    reflection::diffuse<numeric_type> diffuse;
    reflection::i_reflection<numeric_type>& reflection = diffuse;

    std::unique_ptr<ray::i_direction<numeric_type> > srcDirection =
      std::make_unique<ray::cosine_direction_z<numeric_type> > ();
  };
}
