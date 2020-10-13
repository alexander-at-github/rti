#pragma once

#include <map>

#include "../util/utils.hpp"

namespace rti::geo {
  class disc_bounding_box_intersector {
    
  private:
    
    using numeric_type = float;
    // In this code it is helpful that ntriple and nquadruple are arrays with condecutive values because
    // of Embree (c-style) compatibility reasons.
    using ntriple = std::array<numeric_type, 3>;
    using nquadruple = std::array<numeric_type, 4>;

    using bbtype = struct {
      struct {
        numeric_type xx, yy;
      } low, high;
    };

  public:

    disc_bounding_box_intersector
    (numeric_type xmin, numeric_type ymin, numeric_type xmax, numeric_type ymax) :
      bbox({{xmin, ymin}, {xmax, ymax}})
    {
      fill_bboxtransforms(bbox, bboxtransforms);
      assert_invariants();
    }
      
    disc_bounding_box_intersector
    (rti::util::pair<rti::util::pair<numeric_type> > boundingbox) :
      disc_bounding_box_intersector(boundingbox[0][0],
                                    boundingbox[0][1],
                                    boundingbox[1][0],
                                    boundingbox[1][1]) {}

    // This function is thread-safe
    numeric_type
    area_inside (nquadruple& disc, ntriple& dnormal_)
    {
      assert_invariants();
      auto& xx = disc[0];
      auto& yy = disc[1];
      auto& zz = disc[2];
      auto& radius = disc[3];
      // Important: the disc normal needs to be normalized
      auto dnormal = dnormal_;
      rti::util::normalize(dnormal);
      if (radius <= 0) {
        return 0;
      }
      auto fulldiscarea = radius * radius * (numeric_type) rti::util::pi();
      // If fully inside
      if ((bbox.low.xx <= xx - radius && xx + radius <= bbox.high.xx) &&
          (bbox.low.yy <= yy - radius && yy + radius <= bbox.high.yy)) {
        return fulldiscarea;
      }
      // If fully outside
      if ((xx + radius <= bbox.low.xx || bbox.high.xx <= xx - radius) ||
          (yy + radius <= bbox.low.yy || bbox.high.yy <= yy - radius)) {
        return 0;
      }
      auto distobjs = compute_closest_approach(disc, dnormal);
      for (size_t idx = 0; idx < distobjs.size(); ++idx) {
        if (distobjs[idx].approach < -radius) {
          // fully outside
          return 0;
        }
      }
      auto areaoutside = compute_area_outside_given_closest_approach(disc, dnormal, distobjs);
      return fulldiscarea - areaoutside;
    }
    
    void print_bboxtransforms_member()
    {
      for (auto const& swapXY : std::vector<bool> {false, true}) {
        for (auto const& reflectX : std::vector<bool> {false, true}) {
          auto& bb = bboxtransforms[swapXY][reflectX];
          std::cout << "swapXY == " << swapXY << " refelctX == " << reflectX; // << std::endl;
          std::cout << " bb.low.xx == " << bb.low.xx
                    << " bb.low.yy == " << bb.low.yy
                    << " bb.high.xx == " << bb.high.xx
                    << " bb.high.yy == " << bb.high.yy << std::endl;
        }
      }
    }

  private:

    struct transferobj_t {
      numeric_type approach;
      rti::util::pair<bool> bbaccess;
    };

    ntriple
    compute_intersection_point_from_plain_and_x_y
    (ntriple& ppoint, ntriple& pnormal, numeric_type& xx, numeric_type& yy)
    {
      auto& p = ppoint;
      auto& n = pnormal;
      return {xx, yy, (n[0]*p[0] + n[1]*p[1] + n[2]*p[2] - n[0]*xx - n[1]*yy)/n[2]};
    }

    numeric_type
    compute_area_outside_given_closest_approach
    (nquadruple& disc, ntriple& dnormal, rti::util::quadruple<transferobj_t> aobj)
    {
      assert(rti::util::is_normalized(dnormal) && "Precondition");
      auto& radius = disc[3];
      auto area = (numeric_type) 0;
      // Iterate over the directions (x+, y-, x-, y+)
      for (size_t idx = 0; idx < aobj.size(); ++idx) {
        // this is the distance on disc from center of disc to the closest point on the intersection line
        auto& distDCtoCIL = aobj[idx].approach;
        if (-radius < distDCtoCIL && distDCtoCIL < radius) {
          auto angle = 2 * std::acos((double) distDCtoCIL / radius);
          //angle = rti::util::pi() * 2 - angle;
          auto circsegmentarea = radius * radius / 2 * (angle - std::sin(angle));
          area += circsegmentarea;
        }
      }
      // Iterate over the possible overlaps
      for (size_t idx = 0; idx < aobj.size(); ++idx) {
        auto& a1 = aobj[idx];
        auto& a2 = aobj[(idx+1)%aobj.size()];
        auto& d1 = a1.approach;
        auto& d2 = a2.approach;
        auto& swapXY1 = a1.bbaccess[0];
        auto& reflectX1 = a1.bbaccess[1];
        auto& swapXY2 = a2.bbaccess[0];
        auto& reflectX2 = a2.bbaccess[1];
        auto& bbt1 = bboxtransforms[swapXY1][reflectX1];
        auto& bbt2 = bboxtransforms[swapXY2][reflectX2];
        if (-radius < d1 && d1 < radius &&
            -radius < d2 && d2 < radius) {
          // Overlap possible
          auto dpoint = ntriple {disc[0], disc[1], disc[2]};
          // dnormal has the desired form already
          auto bbp1point = ntriple {bbt1.high.xx, bbt1.high.yy, 0};
          auto bbp2point = ntriple {bbt2.high.xx, bbt2.high.yy, 0};
          auto bbp1 = rti::util::triple<rti::util::triple<numeric_type> >
            {rti::util::triple<numeric_type> {bbt1.high.xx, bbt1.high.yy, 1},
             rti::util::triple<numeric_type> {bbt1.high.xx, bbt1.high.yy, 0},
             rti::util::triple<numeric_type> {bbt1.high.xx, bbt1.low.yy , 0}};
          auto bbp2 = rti::util::triple<rti::util::triple<numeric_type> >
            {rti::util::triple<numeric_type> {bbt2.high.xx, bbt2.high.yy, 1},
             rti::util::triple<numeric_type> {bbt2.high.xx, bbt2.high.yy, 0},
             rti::util::triple<numeric_type> {bbt2.high.xx, bbt2.low.yy , 0}};
          auto bbp1normal = rti::util::compute_normal(bbp1);
          auto bbp2normal = rti::util::compute_normal(bbp2);
          rti::util::normalize(bbp1normal);
          rti::util::normalize(bbp2normal);
          if (reflectX1) {
            // Reflect Y axis fist; See comments in fill_bboxtransforms()
            bbp1point[1] *= -1;
            bbp1normal[1] *= -1;
            // refelct X
            bbp1point[0] *= -1;
            bbp1normal[0] *= -1;
          }
          if (reflectX2) {
            // Reflect Y axis first; See comments in fill_bboxtransforms()
            bbp2point[1] *= -1;
            bbp2normal[1] *= -1;
            // refelct X
            bbp2point[0] *= -1;
            bbp2normal[0] *= -1;
          }
          if (swapXY1) {
            // Reflect Y axis first; See comments in fill_bboxtransforms()
            bbp1point[1] *= -1;
            bbp1normal[1] *= -1;
            // swap X and Y
            rti::util::swap(bbp1point[0], bbp1point[1]);
            rti::util::swap(bbp1normal[0], bbp1normal[1]);
          }
          if (swapXY2) {
            // Reflect Y axis first; See comments in fill_bboxtransforms()
            bbp2point[1] *= -1;
            bbp2normal[1] *= -1;
            // swap X and Y
            rti::util::swap(bbp2point[0], bbp2point[1]);
            rti::util::swap(bbp2normal[0], bbp2normal[1]);
          }
          auto idir1 = get_intersection_vector(dnormal, bbp1normal);
          auto idir2 = get_intersection_vector(dnormal, bbp2normal);
          // The cormals bbp1normal and bbp2normal are facing inwards.
          // We want the direction vectors idir1 and idir2 of the intersection to face outwards
          // (with respect to the bounding box).
          // We do that by comparing idir1 with bbp2normal and idir2 with bbp1normal.
          if (vec_same_direction(idir1, bbp2normal)) {
            idir1 = rti::util::inv(idir1);
          }
          if (vec_same_direction(idir2, bbp1normal)) {
            idir2 = rti::util::inv(idir2);
          }
          // bbp2point contains a point on the corner of the bounding box which (possibly)
          // is located in the disc.
          // That is, bbp2point contains the x and y coordinate of the intersection point.
          // When we compute the z axis we will have the point where the two plains of the bounding box
          // and the plain of the disc intersect.
          auto& xinter = bbp2point[0];
          auto& yinter = bbp2point[1];
          auto intersectionpoint = compute_intersection_point_from_plain_and_x_y(dpoint, dnormal, xinter, yinter);
          auto& ipoint1 = intersectionpoint;
          auto& ipoint2 = intersectionpoint;
          if (rti::util::distance(dpoint, intersectionpoint) >= radius) {
            // No overlap
            continue;
          }
          // Definitely an overlap
          assert(rti::util::is_normalized(idir1) && "Assumption");
          auto ipoint1ToCenterOfDisc = rti::util::diff(dpoint, ipoint1); 
          auto closestApproachAlongIdir1 =
            rti::util::dot_product(ipoint1ToCenterOfDisc, idir1);
          auto idir1copy = idir1;
          auto closestPointOnIline1ToDiscCenter =
            rti::util::sum(ipoint1, rti::util::scale(closestApproachAlongIdir1, idir1copy));
          auto thc1 = std::sqrt(radius * radius - d1 * d1);
          auto dirvar1 = idir1; // prepare
          rti::util::scale(thc1, dirvar1);
          // Since idir1 is facing outward (of the bounding box) we know that the correct intersection point
          // is (closestPointOnIline1ToDiscCenter + dirvar1) and not (closestPointOnIline1ToDiscCenter - dirvar1).
          auto q1 = rti::util::sum (closestPointOnIline1ToDiscCenter, dirvar1);
          // q1 is the point on the circumference of the disc which intersects the line defined by the point ipoint1
          // and the direction vector idir1.
          assert(rti::util::is_normalized(idir2) && "Assumption");
          auto ipoint2ToCenterOfDisc = rti::util::diff(dpoint, ipoint2); 
          auto closestApproachAlongIdir2 =
            rti::util::dot_product(ipoint2ToCenterOfDisc, idir2);
          auto idir2copy = idir2;
          auto closestPointOnIline2ToDiscCenter =
            rti::util::sum(ipoint2, rti::util::scale(closestApproachAlongIdir2, idir2copy));
          auto thc2 = std::sqrt(radius * radius - d2 * d2);
          assert(thc2 >= 0 && "Correctness Assertion");
          auto dirvar2 = idir2; // prepare
          rti::util::scale(thc2, dirvar2);
          // Since idir2 is facing outward (of the bounding box) we know that the correct intersection point
          // is (closestPointOnIline2ToDiscCenter + dirvar2) and not (closestPointOnIline2ToDiscCenter - dirvar2).
          auto q2 = rti::util::sum (closestPointOnIline2ToDiscCenter, dirvar2);
          // q2 is the point on the circumference of the disc which intersects the line defined by the point ipoint2
          // and the direction vector idir2.
          auto disccenterTOq1 = rti::util::diff(q1, dpoint);
          auto disccenterTOq2 = rti::util::diff(q2, dpoint);
          auto angle = // angle between disccenterTOq1 and disccenterTOq2
            std::acos(
            rti::util::dot_product(disccenterTOq1, disccenterTOq2)
            / rti::util::length_of_vec(disccenterTOq1)
            / rti::util::length_of_vec(disccenterTOq2));
          auto overlapcircsegmentarea = radius * radius / 2 * (angle - std::sin(angle));
          auto cornerTOq1 = rti::util::diff(q1, intersectionpoint);
          auto cornerTOq2 = rti::util::diff(q2, intersectionpoint);
          auto overlaptrianglearea =
            (numeric_type) 0.5 * rti::util::length_of_vec(rti::util::cross_product(cornerTOq1, cornerTOq2));
          area -= overlapcircsegmentarea + overlaptrianglearea;
        }
      }
      return area;
    }
    
    static void fill_bboxtransforms
    (bbtype& bbox, std::map<bool, std::map<bool, bbtype> >& bboxtransforms)
    {
      bboxtransforms[false] = std::map<bool, bbtype> {};
      bboxtransforms[true]  = std::map<bool, bbtype> {};

      auto swapXY = false;
      auto reflectX = false;
      // When assigning a value to a new map-key the following happens:
      // (1) The map calls the default constructor of the map-value type.
      // (2) The map.operator[] returns a reference to this new object
      // (3) the assignment operator of the map-value type is called.
      bboxtransforms[swapXY][reflectX] = bbox;
      
      swapXY = true;
      reflectX = false;
      bboxtransforms[swapXY][reflectX] = bbox;
      auto* currentP = &bboxtransforms[swapXY][reflectX];
      rti::util::swap(currentP->low.xx , currentP->low.yy);
      rti::util::swap(currentP->high.xx, currentP->high.yy);
      // We also reflect along the Y axis. This makes the top right corner of the
      // transformed bounding box correspond to the corners (in the original
      // bounding box) in clock wise order.
      currentP->low.yy *= -1;
      currentP->high.yy *= -1;

      swapXY = false;
      reflectX = true;
      bboxtransforms[swapXY][reflectX] = bbox;
      currentP = &bboxtransforms[swapXY][reflectX];
      currentP->low.xx  *= -1;
      currentP->high.xx *= -1;
      // We also reflect along the Y axis. This makes the top right corner of the
      // transformed bounding box correspond to the corners (in the original
      // bounding box) in clock wise order.
      currentP->low.yy *= -1;
      currentP->high.yy *= -1;

      swapXY = true;
      reflectX = true;
      bboxtransforms[swapXY][reflectX] = bbox;
      currentP = &bboxtransforms[swapXY][reflectX];
      // First swap, then reflect X values along origin
      rti::util::swap(currentP->low.xx , currentP->low.yy);
      rti::util::swap(currentP->high.xx, currentP->high.yy);
      currentP->low.xx  *= -1;
      currentP->high.xx *= -1;
      // Here we do not have to reflect the Y axis (like in the two cases above), cause
      // we would have to do it twice which cancles it out.
      
      fix_bboxtransforms(bbox, bboxtransforms);
    }

    static void fix_bboxtransforms
    (bbtype& bbox, std::map<bool, std::map<bool, bbtype> >& bboxtransforms)
    {
      // Fix the original bounding box
      if (bbox.low.xx > bbox.high.xx) {
        rti::util::swap(bbox.low.xx, bbox.high.xx);
      }
      if (bbox.low.yy > bbox.high.yy) {
        rti::util::swap(bbox.low.yy, bbox.high.yy);
      }
      // Fix transformed bounding boxes
      for (auto const& v1 : std::vector<bool> {false, true}) {
        for (auto const& v2 : std::vector<bool> {false, true}) {
          auto& currentRef = bboxtransforms[v1][v2];
          if (currentRef.low.xx > currentRef.high.xx) {
            rti::util::swap(currentRef.low.xx, currentRef.high.xx);
          }
          if (currentRef.low.yy > currentRef.high.yy) {
            rti::util::swap(currentRef.low.yy, currentRef.high.yy);
          }
        }
      }
    }

    rti::util::quadruple<transferobj_t>
    compute_closest_approach
    (nquadruple& disc, ntriple& dnormal)
    {
      auto result = rti::util::quadruple<transferobj_t> {};
      auto& radius = disc[3];
      
      // The ordering of the values in the result array is: right, bottom, left, top
      // Note that other parts of the algorithm (the overlap calculation) do _not_
      // work with any ordering.
      auto tuples = rti::util::quadruple<std::pair<size_t, rti::util::pair<bool> > >
        {std::pair<size_t, rti::util::pair<bool> > {0, {false, false}},
         std::pair<size_t, rti::util::pair<bool> > {1, {true , true }},
         std::pair<size_t, rti::util::pair<bool> > {2, {false, true }},
         std::pair<size_t, rti::util::pair<bool> > {3, {true , false}}};
      for (auto& tuple : tuples) {
        auto& idx = tuple.first;
        auto& swapXY = tuple.second[0];
        auto& reflectXcoords = tuple.second[1];
        result[idx].approach = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
        result[idx].bbaccess[0] = swapXY;
        result[idx].bbaccess[1] = reflectXcoords;
        if (result[idx].approach < -radius) {
          // disc is fully outside the bounding box
          return result;
        }
      }
      return result;
    }

    numeric_type
    compute_closest_approach
    (nquadruple& disc_,
     ntriple& dnormal_,
     bool swapXY,
     bool reflectXcoords)
    {
      assert_invariants();
      assert(rti::util::is_normalized(dnormal_) && "Precondition");
      // Prepare: swap or reflect coordinates, if necessary
      auto xidx = 0;
      auto yidx = 1;
      auto zidx = 2;
      if (swapXY) {
        rti::util::swap(xidx, yidx);
      }
      auto xx = disc_[xidx];
      auto yy = disc_[yidx];
      auto zz = disc_[zidx];
      auto radius = disc_[3];
      auto nx = dnormal_[xidx];
      auto ny = dnormal_[yidx];
      auto nz = dnormal_[zidx];
      // First swap then reflect X values along origin
      if (reflectXcoords) {
        xx = -xx;
        nx = -nx;
      }
      // Here reflecint Y values is not necessary. It does not make a difference.
      auto const& bb = bboxtransforms[swapXY][reflectXcoords];

      assert(radius > 0);
      auto xterm =  radius * std::sqrt(nz * nz + ny * ny);
      assert(xterm >= 0);
      auto discXlimithigh = xx + xterm;
      if (discXlimithigh <= bb.high.xx) {
        // disc fully inside
        return std::numeric_limits<numeric_type>::max();
      }
      auto discXlimitlow = xx - xterm;
      if (discXlimitlow >= bb.high.xx) {
        // disc fully outside
        return std::numeric_limits<numeric_type>::lowest();
      }
      
      auto bbXplain = rti::util::triple<rti::util::triple<numeric_type> >
        {rti::util::triple<numeric_type> {bb.high.xx, bb.high.yy, 0},
         rti::util::triple<numeric_type> {bb.high.xx, bb.high.yy, 1},
         rti::util::triple<numeric_type> {bb.high.xx, bb.low.yy , 0}};
      auto dnormForAssertion = ntriple {nx, ny, nz};
      // assert( (xterm <= 1e-9) == (rti::util::normal_perpenticular_to_plain(dnormForAssertion, bbXplain))
      //        && "Correctness Assumption");
      //
      // Candidates to check parallelness of disc to plain:
      // (1) rti::util::normal_perpenticular_to_plain({nx, ny, nz}, bbXplain)
      // (2) normalize {nx, ny, nz} and bbXnormal, compute cross-product, check if length < eps
      // (3) xterm <= eps
      // If disc is parallel to the bounding plain
      if (xterm <= 1e-9) {
        // Disc is parallel to the bounding plain
        return std::numeric_limits<numeric_type>::max();
      }
      // Compute closest approach
      // If (bb.high.xx -xx) is positive the intersection is to the right of the center of the disc,
      // otherwise to the left.
      // This is the closest approach on the disc (not along the x-achsis, which is xterm).
      auto closestapproach = (bb.high.xx - xx) * radius / xterm;
      assert (std::abs(closestapproach) <= radius && "Correctness Assumption");
      return closestapproach;
    }

    ntriple
    get_intersection_vector
    (ntriple& n1, ntriple& n2)
    {
      auto result = rti::util::cross_product(n1, n2);
      rti::util::normalize(result);
      return result;
    }

    bool
    vec_same_direction(ntriple& v1, ntriple& v2)
    {
      return rti::util::dot_product(v1, v2) >= 0;
    }

    void assert_invariants()
    {
      for (auto const& v1 : std::vector<bool> {false, true}) {
        for (auto const& v2 : std::vector<bool> {false, true}) {
          assert(bboxtransforms[v1][v2].low.xx <= bboxtransforms[v1][v2].high.xx && "Assertion");
          assert(bboxtransforms[v1][v2].low.yy <= bboxtransforms[v1][v2].high.yy && "Assertion");
        }
      }
    }

  private:
    
    bbtype bbox;
    std::map<bool, std::map<bool, bbtype> > bboxtransforms;
  };
}
