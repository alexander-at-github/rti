#pragma once

// #include <embree3/rtcore.h>

#include <cmath>

#include "../util/utils.hpp"

namespace rti { namespace trace {
  class local_intersector {
    
  public:

    // // See Glassner page 40
    // originToCenterOfDiscVector = centerOfDisc - originOfRay;
    // ocsqurd = originToCenterOfDiscVector dot-product originToCenterOfDiscVector;
    // closestApproachAlongRay = originToCenterOfDiscVector dot-product rayDirection;
    // closestApproachToCenter = sqrt( ocsqurd * ocsqurd - closestApproachAlongRay * originToCenterOfDiscVector);
    // if (closestApproachToCenter < radius)
    //   hit = true;

    template<typename numeric_type>
    static bool
    does_intersect
    (RTCRay const& ray,
     util::quadruple<numeric_type> const& disc,
     util::triple<numeric_type> const& dnormal)
    {
      // auto rorg = rti::util::triple<numeric_type> {ray.org_x, ray.org_y, ray.org_z};
      // auto rdir = rti::util::triple<numeric_type> {ray.dir_x, ray.dir_y, ray.dir_z};
      static_assert(std::is_same<numeric_type, float>::value, "Casting assumption");
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
