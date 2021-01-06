#pragma once

#include <embree3/rtcore.h>

#include <cmath>

#include "../util/utils.hpp"

namespace rti { namespace trace {
  class local_intersector {
    
  public:

    template<typename numeric_type>
    static bool
    intersect_v2
    (RTCRay const& ray,
     util::quadruple<numeric_type> const& disc,
     util::triple<numeric_type> const& dnormal)
    {
      //embree::avx2::OrientedDiscMiIntersector1<8, 8, true>::intersect
      assert(false && "not implemented");
      return false;
    }
    
    // See Glassner page 50

    // Given a point pp and a normal nn, a plain in 3D is defined by the
    // following equation:
    //
    // 0 = (nn[x] * x + nn[y] * y + nn[z] * z) - dd
    // where dd = - (nn[x] * pp[x] + nn[y] * pp[y] + nn[z] * pp[z])
    //
    // which is equivalent to the following formula in vector notation:
    //
    // 0 = (nn dot-prod (x, y, z)) - dd
    // where dd = - (nn dot-prod pp)
    //
    // Note that dd can be precomputed for each disk.
    
    template<typename numeric_type>
    static bool
    intersect
    (RTCRay const& ray,
     util::quadruple<numeric_type> const& disc,
     util::triple<numeric_type> const& dnormal)
    {
      static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
      // We cast from the Embree data structures which use float
      auto const& rorg = *reinterpret_cast<util::triple<numeric_type> const*> (&ray.org_x);
      auto const& rdir = *reinterpret_cast<util::triple<numeric_type> const*> (&ray.dir_x);
      auto const& dorg = *reinterpret_cast<util::triple<numeric_type> const*> (&disc);

      assert ("Precondition" && util::is_normalized(rdir) && util::is_normalized(dnormal));
      
      auto prodOfDirections = rti::util::dot_product<numeric_type>(dnormal, rdir);
      // std::cerr << "prodOfDirections == " << prodOfDirections << std::endl;
      if (prodOfDirections > 0) {
        // Disc normal is pointing away from the ray direction,
        // i.e., this might be a hit from the back or no hit at all.
        //std::cerr << "possible hit from back" << std::endl;
        return false;
      }
      assert (prodOfDirections <= 0 && "Assumption");
      auto eps = (numeric_type) 1e-9;
      if (std::fabs(prodOfDirections) < eps) {
        // Ray is parallel to disc surface
        //std::cerr << "parallel" << std::endl;
        return false;
      }

      // std::cerr << "### TODO: Memoize ddneg !!! ###" << std::endl;
      auto ddneg = util::dot_product(dorg, dnormal);
      // the nominator term of tt
      auto ttnom = (-util::dot_product(dnormal, rorg)) + ddneg;
      auto tt = ttnom / prodOfDirections;
      if (tt <= 0) {
        // Intersection point is behind or exactly on the ray origin.
        //std::cerr << "behind or exactly on" << std::endl;
        return false;
      }
        
      // copy ray direction
      auto rdirC = rti::util::triple<numeric_type> {rdir[0], rdir[1], rdir[2]};
      // std::cerr << "rdirC == " << rdirC[0] << " " << rdirC[1] << " " << rdirC[2] << std::endl;
      assert( util::is_normalized(rdirC) && "Correctness Assumption");
      util::scale(tt, rdirC);
      auto hitpoint = util::sum(rorg, rdirC);
      auto dorg2HitPoint = util::diff(hitpoint, dorg);
      auto distance = util::length_of_vec(dorg2HitPoint);
      // std::cerr << "length == " << distance << std::endl;
      auto const& radius = disc[3];
      // Be aware that Embree seems to test something like 'if (radius * 1.0001 > distance)'
      if (radius > distance) {
        return true;
      }
      return false;
    }

        template<typename numeric_type>
    static bool
    intersect_debug
    (RTCRay const& ray,
     util::quadruple<numeric_type> const& disc,
     util::triple<numeric_type> const& dnormal)
    {
      static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
      // We cast from the Embree data structures which use float
      auto const& rorg = *reinterpret_cast<util::triple<numeric_type> const*> (&ray.org_x);
      auto const& rdir = *reinterpret_cast<util::triple<numeric_type> const*> (&ray.dir_x);
      auto const& dorg = *reinterpret_cast<util::triple<numeric_type> const*> (&disc);

      assert ("Precondition" && util::is_normalized(rdir) && util::is_normalized(dnormal));
      
      auto prodOfDirections = rti::util::dot_product<numeric_type>(dnormal, rdir);
      std::cerr << "prodOfDirections == " << prodOfDirections << std::endl;
      if (prodOfDirections > 0) {
        // Disc normal is pointing away from the ray direction,
        // i.e., this might be a hit from the back or no hit at all.
        std::cerr << "possible hit from back" << std::endl;
        return false;
      }
      assert (prodOfDirections <= 0 && "Assumption");
      auto eps = (numeric_type) 1e-9;
      if (std::fabs(prodOfDirections) < eps) {
        // Ray is parallel to disc surface
        std::cerr << "parallel" << std::endl;
        return false;
      }

      // std::cerr << "### TODO: Memoize ddneg !!! ###" << std::endl;
      auto ddneg = util::dot_product(dorg, dnormal);
      // the nominator term of tt
      auto ttnom = (-util::dot_product(dnormal, rorg)) + ddneg;
      auto tt = ttnom / prodOfDirections;
      if (tt <= 0) {
        // Intersection point is behind or exactly on the ray origin.
        std::cerr << "behind or exactly on" << std::endl;
        return false;
      }
        
      // copy ray direction
      using internal_numeric_type = double;
      auto rdirC = rti::util::triple<internal_numeric_type> {rdir[0], rdir[1], rdir[2]};
      auto rorgC = rti::util::triple<internal_numeric_type> {rorg[0], rorg[1], rorg[2]};
      auto dorgC = rti::util::triple<internal_numeric_type> {dorg[0], dorg[1], dorg[2]};
      std::cerr << "rdirC == " << rdirC[0] << " " << rdirC[1] << " " << rdirC[2] << std::endl;
      assert( util::is_normalized(rdirC) && "Correctness Assumption");
      util::scale((internal_numeric_type) tt, rdirC);
      auto hitpoint = util::sum(rorgC, rdirC);
      auto dorg2HitPoint = util::diff(hitpoint, dorgC);
      auto distance = util::length_of_vec(dorg2HitPoint);
      auto const& radius = disc[3];
      std::cerr << "length == " << distance << " radius == " << radius << std::endl;
      if (radius > distance) {
        return true;
      }
      return false;
    }

    
    template<typename numeric_type>
    static bool
    intersect_old_v1
    (RTCRay const& ray,
     util::quadruple<numeric_type> const& disc,
     util::triple<numeric_type> const& dnormal)
    {
      // Note: maybe buggy. For sure does not hanlde hits from the back.
      //assert(false);
      
      // auto rorg = rti::util::triple<numeric_type> {ray.org_x, ray.org_y, ray.org_z};
      // auto rdir = rti::util::triple<numeric_type> {ray.dir_x, ray.dir_y, ray.dir_z};
      static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
      // We cast from the Embree data structures which use float
      auto const& rorg = *reinterpret_cast<util::triple<numeric_type> const*> (&ray.org_x);
      auto const& rdir = *reinterpret_cast<util::triple<numeric_type> const*> (&ray.dir_x);
      auto const& dorg = *reinterpret_cast<util::triple<numeric_type> const*> (&disc);

      // std::cerr << "rorg ==  " << rorg[0] << " " << rorg[1] << " " << rorg[2] << std::endl;
      // std::cerr << "ray.o == " << ray.org_x << " " << ray.org_y << " " << ray.org_z << std::endl;
      // std::cerr << "rdir ==  " << rdir[0] << " " << rdir[1] << " " << rdir[2] << std::endl;
      // std::cerr << "dorg == " << dorg[0] << " " << dorg[1] << " " << dorg[2] << std::endl;
      // std::cerr << "dnormal == " << dnormal[0] << " " << dnormal[1] << " " << dnormal[2] << std::endl;

      assert (util::is_normalized(rdir) && util::is_normalized(dnormal));
      auto prodOfDirections = rti::util::dot_product<numeric_type>(rdir, dnormal);
      // std::cerr << "prodOfDirections == " << prodOfDirections << std::endl;
      auto eps = (numeric_type) 1e-4;
      if (std::fabs(prodOfDirections) < eps) {
        //if (prodOfDirections < eps) {
        // Ray is parallel to disc surface
        // std::cerr << "parallel" << std::endl;
        return false;
      }
        
      auto rorg2dorg = rti::util::diff(dorg, rorg);
      // std::cerr << "rorg2dorg: " << rorg2dorg[0] << " "
      //          << rorg2dorg[1] << " " << rorg2dorg[2] << std::endl;
      auto distanceToHitPoint = rti::util::dot_product(rorg2dorg, dnormal) / prodOfDirections;
      // std::cerr << "distanceToHitPoint: " << distanceToHitPoint << std::endl;
      // copy ray direction
      auto rdirC = rti::util::triple<numeric_type> {rdir[0], rdir[1], rdir[2]};
      // std::cerr << "rdirC == " << rdirC[0] << " " << rdirC[1] << " " << rdirC[2] << std::endl;
      util::normalize(rdirC);
      rti::util::scale(distanceToHitPoint, rdirC);
      // std::cerr << "rdirC == " << rdirC[0] << " " << rdirC[1] << " " << rdirC[2] << std::endl;

      { // New
        auto hitpoint = util::sum(rorg, rdirC);
        auto dorg2HitPoint = util::diff(hitpoint, dorg);
        auto distance = util::length_of_vec(dorg2HitPoint);
        // std::cerr << "NEW length == " << distance << std::endl;
        auto const& radius = disc[3];
        if (radius > distance) {
          return true;
        }
        return false;
      }
    }
  };
}}
