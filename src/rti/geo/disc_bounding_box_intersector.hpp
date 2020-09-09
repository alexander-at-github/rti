#pragma once

#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  class disc_bounding_box_intersector {
    using numeric_type = float;
  private:
  public:
    disc_bounding_box_intersector
    (rti::util::pair<rti::util::pair<numeric_type> > boundingbox) :
      //boundingbox(boundingbox),
      bbox({{boundingbox[0][0], boundingbox[0][1]}, {boundingbox[1][0], boundingbox[1][1]}}) {
      if (bbox.low.xx > bbox.high.xx) {
        rti::util::swap(bbox.low.xx, bbox.high.xx);
      }
      if (bbox.low.yy > bbox.high.yy) {
        rti::util::swap(bbox.low.yy, bbox.high.yy);
      }
      assert_invariants();
    }

    numeric_type
    area_inside
    (rti::util::quadruple<numeric_type> disc, rti::util::triple<numeric_type> dnormal)
    {
      assert_invariants();
      auto xx = disc[0];
      auto yy = disc[1];
      auto zz = disc[2];
      auto radius = disc[3];

      auto fulldiscarea = radius * radius * (numeric_type) rti::util::pi();
      
      // If fully outside
      if (xx + radius <= bbox.low.xx || bbox.high.xx <= xx - radius) {
        return 0;
      }
      if (yy + radius <= bbox.low.yy || bbox.high.yy <= yy - radius) {
        return 0;
      }
      // If fully inside
      if ((bbox.low.xx <= xx - radius && xx + radius <= bbox.high.xx) &&
          (bbox.low.yy <= yy - radius && yy + radius <= bbox.high.yy)) {
        return fulldiscarea;
      }
      // If disc parallel to bounding box
      auto bbXplain1 = rti::util::triple<rti::util::triple<numeric_type> >
        {rti::util::triple<numeric_type> {bbox.low.xx, bbox.low.yy, 0},
         rti::util::triple<numeric_type> {bbox.low.xx, bbox.low.yy, 1},
         rti::util::triple<numeric_type> {bbox.low.xx, bbox.high.yy, 0}};
      // TODO: We could check (cross_product <= eps) instead
      if (rti::util::normal_perpenticular_to_plain(dnormal, bbXplain1)) {
        if (bbox.low.xx <= xx && xx <= bbox.high.xx) {
          // TODO: intersect with y-achsis plain
          return fulldiscarea;
        } else {
          return 0;
        }
      }
      // If disc intersects boundary
      auto bbXnormal = rti::util::compute_normal(bbXplain1);
      rti::util::normalize(bbXnormal);
      { // TODO: refactor into function
        auto idir = get_intersection_vector(dnormal, bbXnormal);
        auto dpoint = rti::util::triple<numeric_type> {disc[0], disc[1], disc[2]};
        // right hand plain of bounding box (X)
        auto bbpoint = rti::util::triple<numeric_type> {bbox.high.xx, bbox.high.yy, 0};
        auto ipoint = find_one_intersection_point
          (idir, dnormal, dpoint, bbXnormal, bbpoint);
        // idir and ipoint define the intersection line between the plain of the bounding box and
        // the plain of the disc.
        //std::cerr << "idir == " << idir[0] << " " << idir[1] << " " << idir[2] << std::endl;
        //std::cerr << "ipoint == " << ipoint[0] << " " << ipoint[1] << " " << ipoint[2] << std::endl;

        assert(rti::util::is_normalized(idir) && "Assumption");
        auto ipointToCenterOfDisc = rti::util::diff(dpoint, ipoint); 
        auto closestApproachAlongIdir =
          rti::util::dot_product(ipointToCenterOfDisc, idir);
        auto closestPointOnIlineToDiscCenter =
          rti::util::sum(ipoint, rti::util::scale(closestApproachAlongIdir, idir));
        auto vecCenterOfDiscToClosestIline =
          rti::util::diff(closestPointOnIlineToDiscCenter, dpoint);
        auto distCenterOfDiscToClosestIline =
          rti::util::length_of_vec(vecCenterOfDiscToClosestIline);
        //std::cerr << "distCenterOfDiscToClosestIline == " << distCenterOfDiscToClosestIline << std::endl;
        if (radius <= distCenterOfDiscToClosestIline) {
          return fulldiscarea;
        }

        // Area of circular segment
        // See: https://de.wikipedia.org/wiki/Kreissegment
        assert(distCenterOfDiscToClosestIline < radius && "Assumption");
        auto angle = 2 * std::acos((double) distCenterOfDiscToClosestIline / radius);
        auto areacircularsegment = (double) radius * radius / 2 * (angle - std::sin(angle));
        //std::cout << "areacircularsegment == " << areacircularsegment << std::endl;
        if (disc[0] <= bbox.high.xx) {
          return fulldiscarea - areacircularsegment;
        } else {
          return areacircularsegment;
        }
      }
      
      // Error; for now
      return -1;
    }

    rti::util::triple<numeric_type>
    get_intersection_vector
    (rti::util::triple<numeric_type>& n1,
     rti::util::triple<numeric_type>& n2)
    {
      auto result = rti::util::cross_product(n1, n2);
      rti::util::normalize(result);
      return result;
    }

    rti::util::triple<numeric_type>
    find_one_intersection_point
    (rti::util::triple<numeric_type> idir, // cross product of the two normals
     rti::util::triple<numeric_type> n1, rti::util::triple<numeric_type> p1,
     rti::util::triple<numeric_type> n2, rti::util::triple<numeric_type> p2)
    {
      // See: https://stackoverflow.com/questions/16025620/finding-the-line-along-the-intersection-of-two-planes
      // and  http://geomalgorithms.com/a05-_intersect-1.html

      // Each plain is given by a point p_i and a normal vector n_i.
      // Each plain has an implicit equation n_i \dot (x, y, z) + d_i = 0.
      // To compute d_i we substitute the point p_i for (x, y, z) in the equaion.
      auto ixabs = std::abs(idir[0]);
      auto iyabs = std::abs(idir[1]);
      auto izabs = std::abs(idir[2]);
      // Pick component with highest absolut value
      auto maxidx = 0;
      if (ixabs >= iyabs && ixabs >= izabs)
        maxidx = 0;
      else if (iyabs >= ixabs && iyabs >= izabs)
        maxidx = 1;
      else
        maxidx = 2;

      //std::cerr << "maxidx == " << maxidx << std::endl;

      // Use indices in [0..2] but the maxidx
      //auto useidcs = rti::util::pair<size_t> { (maxidx + 1) % 3, (maxidx + 2) % 3};
      auto i0 = (maxidx + 1) % 3;
      auto i1 = (maxidx + 2) % 3;

      auto dot1 = rti::util::dot_product(n1, p1);
      auto dot2 = rti::util::dot_product(n2, p2);
      
      auto a1 = n1[i0];
      auto b1 = n1[i1];
      auto d1 = - dot1;
      auto a2 = n2[i0];
      auto b2 = n2[i1];
      auto d2 = - dot2;

      auto denominator = (a1 * b2) - (a2 * b1);
      auto a0 = ((b1 * d2) - (b2 * d1)) / denominator;
      auto b0 = ((a2 * d1) - (a1 * d2)) / denominator;

      auto result = rti::util::triple<numeric_type> {0, 0, 0};
      result[i0] = a0;
      result[i1] = b0;
      result[maxidx] = 0;
      
      return result;
    }
      
  private:
    void assert_invariants()
    {
      assert(bbox.low.xx <= bbox.high.xx && "Assertion");
      assert(bbox.low.yy <= bbox.high.yy && "Assertion");
    }

  private:
    // rti::util::pair<rti::util::pair<numeric_type> > boundingbox;
    struct {
      struct {
        numeric_type xx, yy;
      } low, high;
    } bbox;
  };
}}
