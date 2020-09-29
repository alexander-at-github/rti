#pragma once

#include "rti/util/utils.hpp"

namespace rti::geo {
  class disc_bounding_box_intersector {
    
  private:
    
    using numeric_type = float;
    using pair = rti::util::pair<numeric_type>;
    // Using local definitions of triple and quadruple in order to be independent of changes in the utils
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
    (rti::util::pair<rti::util::pair<numeric_type> > boundingbox) :
      bbox ({{boundingbox[0][0], boundingbox[0][1]}, {boundingbox[1][0], boundingbox[1][1]}})
    {
      fill_bboxtransforms(bbox, bboxtransforms);
      assert_invariants();
    }

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
      auto fulldiscarea = radius * radius * (numeric_type) rti::util::pi();
      
      if (radius <= 0) {
        return 0;
      }

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
      //dev_assertion(distobjs, radius);

      auto areaoutside = compute_area_outside_given_closest_approach(disc, dnormal, distobjs);
      return fulldiscarea - areaoutside;
      //std::cerr << "areaoutside == " << areaoutside << std::endl;
      
      // // For now this can handle at most one intersection per disc
      // auto area = fulldiscarea;
      // for (size_t idx = 0; idx < distobjs.size(); ++idx) {
      //   // this is the distance on disc from center of disc to the closest point on the intersection line
      //   auto& distDCtoCIL = distobjs[idx].approach;
      //   if (-radius < distDCtoCIL && distDCtoCIL < radius) {
      //     // auto angle = 2 * std::acos((double) distDCtoCIL / radius);
      //     // auto areacircularsegment = (double) radius * radius / 2 * (angle - std::sin(angle));
      //     // //std::cout << "areacircularsegment == " << areacircularsegment << std::endl;
      //     // return fulldiscarea - areacircularsegment;
      //     auto angle = 2 * std::acos((double) distDCtoCIL / radius);
      //     //angle = rti::util::pi() * 2 - angle;
      //     auto circsegmentarea = (double) radius * radius / 2 * (angle - std::sin(angle));
      //     area -= circsegmentarea;
      //     std::cerr << "circsegmentarea == " << circsegmentarea << std::endl;
      //   }
      // }
      // return area;
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

    numeric_type
    compute_area_outside_given_closest_approach
    (nquadruple& disc, ntriple& dnormal, rti::util::quadruple<transferobj_t> aobj)
    {
      assert(rti::util::is_normalized(dnormal) && "Precondition");
      auto& radius = disc[3];
      std::cout << "radius == " << radius << std::endl;
      auto area = (numeric_type) 0;
      std::cout << "area == " << area << std::endl;
      // Iterate over the directions (x+, y-, x-, y+)
      for (size_t idx = 0; idx < aobj.size(); ++idx) {
        // this is the distance on disc from center of disc to the closest point on the intersection line
        auto& distDCtoCIL = aobj[idx].approach;
        if (-radius < distDCtoCIL && distDCtoCIL < radius) {
          // auto angle = 2 * std::acos((double) distDCtoCIL / radius);
          // auto areacircularsegment = (double) radius * radius / 2 * (angle - std::sin(angle));
          // //std::cout << "areacircularsegment == " << areacircularsegment << std::endl;
          // return fulldiscarea - areacircularsegment;
          auto angle = 2 * std::acos((double) distDCtoCIL / radius);
          //angle = rti::util::pi() * 2 - angle;
          auto circsegmentarea = (numeric_type) radius * radius / 2 * (numeric_type) (angle - std::sin(angle));
          area += circsegmentarea;
        }
      }
      std::cout << "area == " << area << std::endl;
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
        std::cerr << "d1 == " << d1 << " d2 == " << d2 << std::endl;
        std::cerr << "swapXY1 == " << (swapXY1 ? "true" : "false") << " "
                  << "reflectX1 == " << (reflectX1 ? "true" : "false") << std::endl;
        std::cerr << "swapXY2 == " << (swapXY2 ? "true" : "false") << " "
                  << "reflectX2 == " << (reflectX2 ? "true" : "false") << std::endl;
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

          // auto bbp1_ = bbp1;
          // if (reflectX1) {
          //   //p1point[0] *= -1;
          //   bbp1_[0][0] *= -1;
          //   bbp1_[1][0] *= -1;
          //   bbp1_[2][0] *= -1;
          // }
          // auto bbp2_ = bbp2;
          // if (reflectX2) {
          //   //p2point[0] *= -1;
          //   bbp2_[0][0] *= -1;
          //   bbp2_[1][0] *= -1;
          //   bbp2_[2][0] *= -1;
          // }
          // auto bbp1normal_ = rti::util::compute_normal(bbp1);
          // auto bbp2normal_ = rti::util::compute_normal(bbp2);
          // rti::util::normalize(bbp1normal_);
          // rti::util::normalize(bbp2normal_);
          // std::cerr << "bbp1normal_ == " << bbp1normal_[0] << " " << bbp1normal_[1] << " " << bbp1normal_[2] << std::endl;
          // std::cerr << "bbp2normal_ == " << bbp2normal_[0] << " " << bbp2normal_[1] << " " << bbp2normal_[2] << std::endl;

            
          
          // Fix bbp1point, bbp2point, bbp1normal, and bbp2normal according to the values
          // of the swapXY and reflectX variables.
          // This is necessary in order to work within the same basis.
          // Because we are unding the swap and reflect from the transforms of the
          // bounding box, we have to first reflect along the X axis and then swap the
          // X and Y values.
          if (reflectX1) {
            bbp1point[0] *= -1;
            bbp1normal[0] *= -1;
          }
          if (reflectX2) {
            bbp2point[0] *= -1;
            bbp2normal[0] *= -1;
          }
          if (swapXY1) {
            rti::util::swap(bbp1point[0], bbp1point[1]);
            rti::util::swap(bbp1normal[0], bbp1normal[1]);
          }
          if (swapXY2) {
            rti::util::swap(bbp2point[0], bbp2point[1]);
            rti::util::swap(bbp2normal[0], bbp2normal[1]);
          }

          std::cerr << "dnormal == " << dnormal[0] << " " << dnormal[1] << " " << dnormal[2] << std::endl;
          std::cerr << "bbp1normal == " << bbp1normal[0] << " " << bbp1normal[1] << " " << bbp1normal[2] << std::endl;
          std::cerr << "bbp2normal == " << bbp2normal[0] << " " << bbp2normal[1] << " " << bbp2normal[2] << std::endl;
          std::cerr << "bbp1point == " << bbp1point[0] << " " << bbp1point[1] << " " << bbp1point[2] << std::endl;
          std::cerr << "bbp2point == " << bbp2point[0] << " " << bbp2point[1] << " " << bbp2point[2] << std::endl;
          
          auto idir1 = get_intersection_vector(dnormal, bbp1normal);
          auto idir2 = get_intersection_vector(dnormal, bbp2normal);
          // The cormals bbp1normal and bbp2normal are facing inwards.
          // We want the direction vectors idir1 and idir2 of the intersection to face outwards (with respect to
          // the bounding box).
          // We do that by comparing idir1 with bbp2normal and idir2 with bbp1normal.
          if (vec_same_direction(idir1, bbp2normal)) {
            idir1 = rti::util::inv(idir1);
          }
          if (vec_same_direction(idir2, bbp1normal)) {
            idir2 = rti::util::inv(idir2);
          }
          
          auto ipoint1 = find_one_intersection_point(idir1, dnormal, dpoint, bbp1normal, bbp1point);
          auto ipoint2 = find_one_intersection_point(idir2, dnormal, dpoint, bbp2normal, bbp2point);
          std::cerr << "idir1 == " << idir1[0] << " " << idir1[1] << " " << idir1[2] << std::endl;
          std::cerr << "idir2 == " << idir2[0] << " " << idir2[1] << " " << idir2[2] << std::endl;
          // We actually know the x and y coordinates of the intersection point. It is saved
          // in bbt1.high.xx and bbt1.high.yy. We just would need to comput the z-value.
          // TODO: improve
          std::cout << "TODO: improve the following call." << std::endl;
          auto intersectionpoint = compute_intersection_point(ipoint1, idir1, ipoint2, idir2);
          //
          
          if (rti::util::distance(dpoint, intersectionpoint) >= radius) {
            // No overlap
            continue;
          }
          // Definitely an overlap
          std::cerr
            << "intersectionpoint == " << intersectionpoint[0] << " "
            << intersectionpoint[1] << " " << intersectionpoint[2] << std::endl;

          assert(rti::util::is_normalized(idir1) && "Assumption");
          auto ipoint1ToCenterOfDisc = rti::util::diff(dpoint, ipoint1); 
          auto closestApproachAlongIdir1 =
            rti::util::dot_product(ipoint1ToCenterOfDisc, idir1);
          std::cerr << "idir1 == " << idir1[0] << " " << idir1[1] << " " << idir1[2] << std::endl;
          auto idir1copy = idir1;
          auto closestPointOnIline1ToDiscCenter =
            rti::util::sum(ipoint1, rti::util::scale(closestApproachAlongIdir1, idir1copy));
          std::cout << "closestPointOnIline1ToDiscCenter == "
                    << closestPointOnIline1ToDiscCenter[0] << " "
                    << closestPointOnIline1ToDiscCenter[1] << " "
                    << closestPointOnIline1ToDiscCenter[2] << std::endl;
          auto thc1 = std::sqrt(radius * radius - d1 * d1);
          std::cout << "thc1 == " << thc1 << std::endl;
          auto dirvar1 = idir1; // prepare
          std::cerr << "dirvar1 == " << dirvar1[0] << " " << dirvar1[1] << " " << dirvar1[2] << std::endl;
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

          // // Choose q1 and q2
          // auto q1 = ntriple {}; // initialze to any value
          // auto q2 = ntriple {};
          // // TODO: move that decision to the idir1 and idir2 variables!
          // q1 = ! vec_same_direction(idir1, bbp2normal) ? q11 : q12;
          // q2 = ! vec_same_direction(idir2, bbp1normal) ? q21 : q22;
          //
          // // In case the disc is centered exactly on the corner of the bounding box and the
          // // disc is perpenticular to both plains of the bounding box, then this method
          // // may chose wrong points q1 and q2. Anyway, this is not a problem, since then
          // // the angle is simply 90 degrees and any pair of points fullfilling this
          // // property will produce a correct area value.
          // auto candidates = rti::util::quadruple<rti::util::pair<decltype(q1)> >
          //   {rti::util::pair<decltype(q1) > { q11, q21 },
          //    rti::util::pair<decltype(q1) > { q11, q22 },
          //    rti::util::pair<decltype(q1) > { q12, q21 },
          //    rti::util::pair<decltype(q1) > { q12, q22 }};
          // auto mindist = std::numeric_limits<numeric_type>::max();
          // for (auto const& pair : candidates) {
          //   auto currdist = rti::util::distance(pair[0], pair[1]);
          //   if (currdist < mindist) {
          //     std::cout << "setting q1 and q2" << std::endl;
          //     mindist = currdist;
          //     q1 = pair[0];
          //     q2 = pair[1];
          //   }
          // }
          //
          // if (rti::util::distance(intersectionpoint, q11) < rti::util::distance(intersectionpoint, q12)) {
          //   q1 = q11;
          // } else {
          //   q1 = q12;
          // }
          // if (rti::util::distance(intersectionpoint, q21) < rti::util::distance(intersectionpoint, q22)) {
          //   q2 = q21;
          // } else {
          //   q1 = q22;
          // }
          std::cout << "q1 == " << q1[0] << " " << q1[1] << " " << q1[2] << std::endl;
          std::cout << "q2 == " << q2[0] << " " << q2[1] << " " << q2[2] << std::endl;
          auto disccenterTOq1 = rti::util::diff(q1, dpoint);
          auto disccenterTOq2 = rti::util::diff(q2, dpoint);
          auto angle = // angle between disccenterTOq1 and disccenterTOq2
            std::acos(
            rti::util::dot_product(disccenterTOq1, disccenterTOq2)
            / rti::util::length_of_vec(disccenterTOq1)
            / rti::util::length_of_vec(disccenterTOq2));
          std::cout << "angle == " << angle << std::endl;
            
          auto overlapcircsegmentarea = radius * radius / 2 * (angle - std::sin(angle));
          std::cerr << "overlapcircsegmentarea == " << overlapcircsegmentarea << std::endl;
          auto cornerTOq1 = rti::util::diff(q1, intersectionpoint);
          auto cornerTOq2 = rti::util::diff(q2, intersectionpoint);
          auto overlaptrianglearea =
            (float) 0.5 * rti::util::length_of_vec(rti::util::cross_product(cornerTOq1, cornerTOq2));
          std::cerr << "overlaptrianglearea == " << overlaptrianglearea << std::endl;

          area -= overlapcircsegmentarea + overlaptrianglearea;
        }
      }
      std::cout << "area == " << area << std::endl;
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
      // (3) with <=> the assignment operator of the map-value type is called.
      bboxtransforms[swapXY][reflectX] = bbox;
      
      swapXY = true;
      reflectX = false;
      bboxtransforms[swapXY][reflectX] = bbox;
      auto* currentP = &bboxtransforms[swapXY][reflectX];
      rti::util::swap(currentP->low.xx , currentP->low.yy);
      rti::util::swap(currentP->high.xx, currentP->high.yy);

      swapXY = false;
      reflectX = true;
      bboxtransforms[swapXY][reflectX] = bbox;
      currentP = &bboxtransforms[swapXY][reflectX];
      currentP->low.xx  = -currentP->low.xx;
      currentP->high.xx = -currentP->high.xx;

      swapXY = true;
      reflectX = true;
      bboxtransforms[swapXY][reflectX] = bbox;
      currentP = &bboxtransforms[swapXY][reflectX];
      // First swap, then reflect X values along origin
      rti::util::swap(currentP->low.xx , currentP->low.yy);
      rti::util::swap(currentP->high.xx, currentP->high.yy);
      currentP->low.xx  = -currentP->low.xx;
      currentP->high.xx = -currentP->high.xx;

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
      // auto swapXY = false;
      // auto reflectXcoords = false;
      // swapXY = false;
      // reflectXcoords = false;
      // result[0].approach = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      // result[0].bbaccess[0] = swapXY;
      // result[0].bbaccess[1] = reflectXcoords;
      // if (result[0] < -radius) return result;
      // swapXY = true;
      // reflectXcoords = true;
      // result[1].approach = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      // result[1].bbaccess[0] = swapXY;
      // result[1].bbaccess[1] = reflectXcoords;
      // if (result[1] < -radius) return result;
      // swapXY = false;
      // reflectXcoords = true;
      // result[2].approach = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      // result[2].bbaccess[0] = swapXY;
      // result[2].bbaccess[1] = reflectXcoords;
      // if (result[2] < -radius) return result;
      // swapXY = true;
      // reflectXcoords = false;
      // result[3].approach = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      // result[3].bbaccess[0] = swapXY;
      // result[3].bbaccess[1] = reflectXcoords;
      // return result;
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
      auto const& bb = bboxtransforms[swapXY][reflectXcoords];

      assert(radius > 0);
      auto xterm =  radius * std::sqrt(nz * nz + ny * ny);
      // std::cerr << "swapXY == " << (swapXY ? "T" : "F")
      //           << " reflectXcoords == " << (reflectXcoords ? "T" : "F") << std::endl;
      // std::cerr << "nx == " << nx << " ny == " << ny << " nz == " << nz << std::endl;
      // std::cerr << "xterm == " << " " << xterm << std::endl;
      // std::cerr << (xterm <= 1e-9 ? "xterm <= 1e-9" : "not (xterm <= 1e-9)") << std::endl;
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
      assert( (xterm <= 1e-9) == (rti::util::normal_perpenticular_to_plain(dnormForAssertion, bbXplain))
              && "Correctness Assumption");
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
    compute_intersection_point
    (ntriple& ipoint1, ntriple& idir1, ntriple& ipoint2, ntriple& idir2)
    {
      std::cerr << " ipoint1 == " << ipoint1[0] << " " << ipoint1[1] << " " << ipoint1[1]
                << " idir1 == "   << idir1[0] << " " << idir1[1] << " " << idir1[2]
                << " ipoint2 == " << ipoint2[0] << " " << ipoint2[1] << " " << ipoint2[1]
                << " idir2 == "   << idir2[0] << " " << idir2[1] << " " << idir2[2] << std::endl;
      
      // https://math.stackexchange.com/questions/270767/find-intersection-of-two-3d-lines/271366
      if (ipoint1[0] == ipoint2[0] && ipoint1[1] == ipoint2[1] && ipoint1[2] == ipoint2[2]) {
        // ipoint1 (or ipoint2) is already the intersection point
        return ipoint1;
      }
      auto eps = 1e-6;
      auto vecp1p2 = rti::util::diff(ipoint2, ipoint1);
      auto id2cpvecp1p2 = rti::util::cross_product(idir2, vecp1p2);
      assert(rti::util::length_of_vec(id2cpvecp1p2) > eps && "Correctness Assertion");
      auto id2cpid1 = rti::util::cross_product(idir2, idir1);
      assert(rti::util::length_of_vec(id2cpid1) > eps && "Correctness Assertion");
      std::cerr << "rti::util::length_of_vec(id2cpvecp1p2) == " << rti::util::length_of_vec(id2cpvecp1p2) << std::endl;
      auto samedir = vec_same_direction(id2cpvecp1p2, id2cpid1);
      auto scalevalue = rti::util::length_of_vec(id2cpvecp1p2) / rti::util::length_of_vec(id2cpid1);
      auto offsetvec = idir1; // prepare
      rti::util::scale(scalevalue, offsetvec); // scale
      if ( ! samedir ) {
        offsetvec = rti::util::inv(offsetvec);
      }
      auto result = rti::util::sum(ipoint1, offsetvec);
      return result;
    }

    template<typename tt> bool same_sgn(tt& v1, tt& v2)
    {
      return (v1 < 0) == (v2 < 0);
    }
    
    bool
    vec_same_direction(ntriple& v1, ntriple& v2)
    {
      //return same_sgn(v1[0], v2[0]) && same_sgn(v1[1], v2[1]) && same_sgn(v1[2], v2[2]);
      return rti::util::dot_product(v1, v2) >= 0;
    }

    void dev_assertion(nquadruple& distobjs, numeric_type& radius)
    {
      auto cnt = 0u;
      for (size_t idx = 0; idx < distobjs.size(); ++idx) {
        assert(distobjs[idx] >= -radius && "Correctness Assertion");
        if (-radius <= distobjs[idx] && distobjs[idx] <= radius) {
          // std::cerr << "idx == " << idx << std::endl;
          cnt += 1;
        }
      }
      assert(cnt <= 1 && "Development Assumption");
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



    
//////////////////////////////////////////////////
// BELOW: Code which might be needed in the future.
//////////////////////////////////////////////////

  private:

    using transferobj = struct {
      numeric_type closestApproach;
      ntriple vecCntrOfDisc2ClosestIline;
      ntriple ilineDir;
    };
    

    transferobj
    compute_closest_approach_old
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
      auto bb = bboxtransforms[swapXY][reflectXcoords];

      // If fully inside the higher bounding plain towads the x-coordinate
      assert(radius > 0);
      auto xterm = std::abs(nz) * radius / std::sqrt(nz * nz + nx * nx);
      // std::cerr << "swapXY == " << (swapXY ? "T" : "F") << " reflectXcoords == " << (reflectXcoords ? "T" : "F") << std::endl;
      // std::cerr << "nx == " << nx << " ny == " << ny << " nz == " << nz << std::endl;
      // std::cerr << nz << " " << nx << " " << std::sqrt(nz * nz + nx * nx) << " " << xterm << std::endl;
      // std::cerr << (xterm <= 1e-9 ? "xterm <= 1e-9" : "not (xterm <= 1e-9)") << std::endl;
      assert(xterm >= 0);
      auto discXlimithigh = xx + xterm;
      if (discXlimithigh <= bb.high.xx) {
        // disc fully inside
        return {std::numeric_limits<numeric_type>::max(), {0, 0, 0}, {0, 0, 0}};
      }
      auto discXlimitlow = xx - xterm;
      if (discXlimitlow >= bb.high.xx) {
        // disc fully outside
        return {std::numeric_limits<numeric_type>::lowest(), {0, 0, 0}, {0, 0, 0}};
      }

      
      auto bbXplain = rti::util::triple<rti::util::triple<numeric_type> >
        {rti::util::triple<numeric_type> {bb.high.xx, bb.high.yy, 0},
         rti::util::triple<numeric_type> {bb.high.xx, bb.high.yy, 1},
         rti::util::triple<numeric_type> {bb.high.xx, bb.low.yy , 0}};

      auto dnormForAssertion = ntriple {nx, ny, nz};
      assert( (xterm <= 1e-9) == (rti::util::normal_perpenticular_to_plain(dnormForAssertion, bbXplain))
              && "Correctness Assumption");
      // Candidates to check parallelness of disc to plain:
      // (1) rti::util::normal_perpenticular_to_plain({nx, ny, nz}, bbXplain)
      // (2) normalize {nx, ny, nz} and bbXnormal, compute cross-product, check if length < eps
      // (3) xterm <= eps
      // If disc is parallel to the bounding plain
      if (xterm <= 1e-9) {
        // Disc is parallel to the bounding plain
        return {std::numeric_limits<numeric_type>::max(), {0, 0, 0}, {0, 0, 0}};
      }
      
      // Compute disc-boundary intersectsion
      auto dpoint = ntriple {xx, yy, zz};
      auto dnormal = ntriple {nx, ny, nz};
      // right hand plain of bounding box (X)
      auto bbpoint = ntriple {bb.high.xx, bb.high.yy, 0};
      auto bbXnormal = rti::util::compute_normal(bbXplain);
      rti::util::normalize(bbXnormal);
      auto idir = get_intersection_vector(dnormal, bbXnormal);
      assert(rti::util::is_normalized(idir) && "Assumption");

      auto ipoint = find_one_intersection_point(idir, dnormal, dpoint, bbXnormal, bbpoint);
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
      if (distCenterOfDiscToClosestIline >= radius) {
        return {std::numeric_limits<numeric_type>::max(), {0, 0, 0}, {0, 0, 0}};
      }
      assert(0 <= distCenterOfDiscToClosestIline && distCenterOfDiscToClosestIline < radius
             && "Correctness Assumption");
      return {distCenterOfDiscToClosestIline, vecCenterOfDiscToClosestIline, idir};
    }

    ntriple
    get_intersection_vector
    (ntriple& n1, ntriple& n2)
    {
      auto result = rti::util::cross_product(n1, n2);
      rti::util::normalize(result);
      return result;
    }

    ntriple
    find_one_intersection_point
    (ntriple& idir, // cross product of the two normals
     ntriple& n1, ntriple& p1,
     ntriple& n2, ntriple& p2)
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

      auto result = ntriple {0, 0, 0};
      result[i0] = a0;
      result[i1] = b0;
      result[maxidx] = 0;
      
      return result;
    }

    
  };
}
