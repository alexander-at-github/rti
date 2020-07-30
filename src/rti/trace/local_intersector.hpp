#pragma once

#include <embree3/rtcore.h>

#include "rti/util/utils.hpp"

namespace rti { namespace trace {
  template<typename numeric_type, typename disc_type, typename normal_type>
  class local_intersector {
  public:
    
    local_intersector (size_t numpoints, disc_type* discs, normal_type* normals) :
      discs(discs),
      normals(normals),
      numpoints(numpoints) {
      init();
    }

    std::vector<size_t> intersect_neighbors(RTCRayHit& rayhit)
    {
      // // See Glassner page 40
      // originToCenterOfDiscVector = centerOfDisc - originOfRay;
      // ocsqurd = originToCenterOfDiscVector dot-product originToCenterOfDiscVector;
      // closestApproachAlongRay = originToCenterOfDiscVector dot-product rayDirection;
      // closestApproachToCenter = sqrt( ocsqurd * ocsqurd - closestApproachAlongRay * originToCenterOfDiscVector);
      // if (closestApproachToCenter < radius)
      //   hit = true;

      auto ray = rayhit.ray;
      auto rorg = rti::util::triple<numeric_type> {ray.org_x, ray.org_y, ray.org_z};
      auto rdir = rti::util::triple<numeric_type> {ray.dir_x, ray.dir_y, ray.dir_z};
      auto result = std::vector<size_t> {};
      result.reserve(8);

      for (size_t idx = 0; idx < numpoints; ++idx) {

        auto dorg = rti::util::triple<numeric_type> {discs[idx].xx, discs[idx].yy, discs[idx].zz};
        auto dnorm = rti::util::triple<numeric_type> {normals[idx].xx, normals[idx].yy, normals[idx].zz};

        auto prodOfDirections = rti::util::dot_product<numeric_type>(rorg, dorg);
        auto eps = (numeric_type) 1e-4;
        if (prodOfDirections < eps) {
          // Ray is parallel to disc surface
          continue;
        }
        
        auto rorg2dorg = rti::util::diff(dorg, rorg);
        auto distanceToHitPoint = rti::util::dot_product(rorg2dorg, dnorm) / prodOfDirections;
        // copy ray direction
        auto rdirC = rti::util::triple<numeric_type> {rdir[0], rdir[1], rdir[2]};
        auto rorg2HitPoint = rti::util::sum(rorg, rti::util::scale(distanceToHitPoint, rdirC));
        auto dorg2HitPoint = rti::util::diff(rorg2HitPoint, rorg2dorg);

        auto distSqrd = dorg2HitPoint[0] * dorg2HitPoint[0] +
          dorg2HitPoint[1] * dorg2HitPoint[1] +
          dorg2HitPoint[2] * dorg2HitPoint[2];

        auto radiusSqrd = discs[idx].radius * discs[idx].radius;
        if (radiusSqrd > distSqrd) {
          result.push_back(idx);
        }
      }
      return result;
    }

  private:

    void init() {}

  private:
    
    disc_type* discs;
    normal_type* normals;
    size_t numpoints;
  };
}}


      // for (size_t idx = 0; idx < numpoints; ++idx) {

      //   auto cd = rti::util::triple<numeric_type> {discs[idx].xx, discs[idx].yy, discs[idx].zz};

      //   auto org2cd = rti::util::diff(cd, rorg);
      //   auto org2cdSqrd = rti::util::dot_product(org2cd, org2cd);
      //   auto closestApproachOnRay = rti::util::dot_product(org2cd, rdir);
      //   // Note: this does not include the normals. It assumes the disc normals are alligned with the ray.
      //   auto closestApproachToCenterSqrd = org2cdSqrd * org2cdSqrd - closestApproachOnRay * closestApproachOnRay;
      //   //auto closestApproachToCenter = std::sqrt(org2cdSqrd * org2cdSqrd - closestApproachOnRay * closestApproachOnRay);

      //   auto radiusSqrd = discs[idx].radius * discs[idx].radius;
      //   if (radiusSqrd > closestApproachToCenterSqrd) {
      //   //if (discs[idx].radius > closestApproachToCenter) {
      //     result.push_back(idx);
      //   }
      // }
