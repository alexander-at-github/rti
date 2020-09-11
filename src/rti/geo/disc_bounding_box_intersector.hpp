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
        if (distobjs[idx] < -radius) {
          // fully outside
          return 0;
        }
      }
      dev_assertion(distobjs, radius);

      // For now this can handle at most one intersection per disc
      for (size_t idx = 0; idx < distobjs.size(); ++idx) {
        // this is the distance on disc from center of disc to the closest point on the intersection line
        auto& distDCtoCIL = distobjs[idx];
        if (-radius < distDCtoCIL && distDCtoCIL < radius) {
          // auto angle = 2 * std::acos((double) distDCtoCIL / radius);
          // auto areacircularsegment = (double) radius * radius / 2 * (angle - std::sin(angle));
          // //std::cout << "areacircularsegment == " << areacircularsegment << std::endl;
          // return fulldiscarea - areacircularsegment;
          auto angle = 2 * std::acos((double) distDCtoCIL / radius);
          angle = rti::util::pi() * 2 - angle;
          auto area = (double) radius * radius / 2 * (angle - std::sin(angle));
          return area;
        }
      }
      return fulldiscarea;
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

    nquadruple
    compute_closest_approach
    (nquadruple& disc, ntriple& dnormal)
    {
      auto result = nquadruple {};
      auto& radius = disc[3];
      
      auto swapXY = false;
      auto reflectXcoords = false;

      // The ordering of the values in the result array is: right, bottom, left, top
      swapXY = false;
      reflectXcoords = false;
      result[0] = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      if (result[0] < -radius) return result;
      swapXY = true;
      reflectXcoords = true;
      result[1] = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      if (result[1] < -radius) return result;
      swapXY = false;
      reflectXcoords = true;
      result[2] = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
      if (result[2] < -radius) return result;
      swapXY = true;
      reflectXcoords = false;
      result[3] = compute_closest_approach(disc, dnormal, swapXY, reflectXcoords);
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
