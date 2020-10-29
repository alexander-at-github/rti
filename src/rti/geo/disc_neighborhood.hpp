#pragma once

namespace rti { namespace geo {
  template<typename numeric_type>
  class disc_neighborhood {
    
  public:

    // disc_neighborhood() {}

    void setup_neighborhood_naive
    (std::vector<rti::util::quadruple<numeric_type> >& points)
    {
      nbhd.clear();
      nbhd.resize(points.size(), std::vector<size_t> {});
      construct_neighborhood_naive_2(points);
    }


    void setup_neighborhood
    (std::vector<rti::util::quadruple<numeric_type> >& points,
     rti::util::triple<numeric_type>& min,
     rti::util::triple<numeric_type>& max)
    {
      nbhd.clear();
      nbhd.resize(points.size(), std::vector<size_t> {});
      construct_neighborhood_2(points, min, max);
    }

  private:

    void construct_neighborhood_2
    (std::vector<rti::util::quadruple<numeric_type> >& points,
     rti::util::triple<numeric_type>& min,
     rti::util::triple<numeric_type>& max)
    {
      std::cout << "points.size() == " << points.size() << std::endl;
      std::cout << "Starting devide and conquer neirest neighbor" << std::endl;
      auto diridx = 0;
      auto pivot = (max[diridx] + min[diridx]) / 2;
      auto s1maxrad = (numeric_type) 0;
      auto s2maxrad = (numeric_type) 0;
      auto s1 = std::vector<size_t> {};
      auto s2 = std::vector<size_t> {};
      for (size_t idx = 0; idx < points.size(); ++idx) {
        auto const& pp = points[idx];
        auto const& radius = pp[3];
        if (pp[diridx] <= pivot) {
          s1.push_back(idx);
          if (radius > s1maxrad) {
            s1maxrad = radius;
          }
        } else {
          s2.push_back(idx);
          if (radius > s2maxrad) {
            s2maxrad = radius;
          }
        }
      }
      divide_and_conquer(points, s1, s2, s1maxrad, s2maxrad, min, max, diridx, pivot);
    }

    void divide_and_conquer
    (std::vector<rti::util::quadruple<numeric_type> > const& pointdata,
     std::vector<size_t> const& s1,
     std::vector<size_t> const& s2,
     numeric_type const& s1maxrad,
     numeric_type const& s2maxrad,
     rti::util::triple<numeric_type> const& min,
     rti::util::triple<numeric_type> const& max,
     int const& diridx,
     numeric_type const& pivot)
    {
      assert (0 <= diridx && diridx < 3 && "Assumption");
      if (s1.size() + s2.size() <= 1) {
        return;
      }
      // sets of candidates
      auto s1c = std::vector<size_t> {};
      auto s2c = std::vector<size_t> {};

      auto newdiridx = (diridx + 1) % 3;
      auto newpivot = (max[newdiridx] + min[newdiridx]) / 2;
      auto s1r1maxr = (numeric_type) 0;
      auto s1r2maxr = (numeric_type) 0;
      auto s2r1maxr = (numeric_type) 0;
      auto s2r2maxr = (numeric_type) 0;
      // recursion sets
      auto s1r1set = std::vector<size_t> {};
      auto s1r2set = std::vector<size_t> {};
      auto s2r1set = std::vector<size_t> {};
      auto s2r2set = std::vector<size_t> {};

      for (size_t idx = 0; idx < s1.size(); ++idx) {
        auto& ee = pointdata[s1[idx]];
        assert(ee[diridx] <= pivot && "Correctness Assertion");
        auto& radius = ee[3];
        if (ee[newdiridx] <= newpivot) {
          s1r1set.push_back(s1[idx]);
          if (radius > s1r1maxr) {
            s1r1maxr = radius;
          }
        } else {
          s1r2set.push_back(s1[idx]);
          if (radius > s1r2maxr) {
            s1r2maxr = radius;
          }
        }
        if (ee[diridx] + radius + s2maxrad < pivot) {
          continue;
        }
        s1c.push_back(s1[idx]);
      }
      for (size_t idx = 0; idx < s2.size(); ++idx) {
        auto& ee = pointdata[s2[idx]];
        assert(ee[diridx] > pivot && "Correctness Assertion");
        auto& radius = ee[3];
        if (ee[newdiridx] <= newpivot) {
          s2r1set.push_back(s2[idx]);
          if (radius > s2r1maxr) {
            s2r1maxr = radius;
          }
        } else {
          s2r2set.push_back(s2[idx]);
          if (radius > s2r2maxr) {
            s2r2maxr = radius;
          }
        }
        if (ee[diridx] - radius - s1maxrad > pivot) {
          continue;
        }
        s2c.push_back(s2[idx]);
      }
      // std::cout << "s1c.size() == " << s1c.size() << std::endl;
      // std::cout << "s2c.size() == " << s2c.size() << std::endl;
      
      // Iterate over pairs of candidates
      if (s1c.size() > 0 && s2c.size() > 0) {
        for (size_t ci1 = 0; ci1 < s1c.size(); ++ci1) {
          for (size_t ci2 = ci1; ci2 < s2c.size(); ++ci2) {
            assert(std::abs(pointdata[s1c[ci1]][diridx] - pointdata[s2c[ci2]][diridx]) <= (2*(s1maxrad + s2maxrad)) &&
                   "Correctness Assertion");
            if ( check_dist(pointdata, s1c[ci1], s2c[ci2], diridx) ) {
              nbhd[s1c[ci1]].push_back(s1c[ci1]);
              nbhd[s2c[ci2]].push_back(s2c[ci2]);
            }
          }
        }
      }
      // Recurse
      if (s1.size() > 1) {
        auto news1max = max;
        news1max[diridx] = pivot; // old diridx and old pivot!
        divide_and_conquer(pointdata, s1r1set, s1r2set, s1r1maxr, s1r2maxr, min, news1max, newdiridx, newpivot);
      }
      if (s2.size() > 1) {
        auto news2min = min;
        news2min[diridx] = pivot; // old diridx and old pivot!
        divide_and_conquer(pointdata, s2r1set, s2r2set, s2r1maxr, s2r2maxr, news2min, max, newdiridx, newpivot);
      }
    }

   bool check_dist
   (std::vector<rti::util::quadruple<numeric_type> > const& pointdata,
    size_t const& i1,
    size_t const& i2,
    int const& diridx)
    {
      auto& p1 = pointdata[i1];
      auto& r1 = p1[3];
      auto& p2 = pointdata[i2];
      auto& r2 = p2[3];
      auto sr = r1 + r2;
      assert(p2[diridx] - p1[diridx] >= 0 && "Assumption");
      if (std::abs(p2[diridx] - p1[diridx]) >= sr) {
        return false;
      }
      if (std::abs(p1[diridx + 1 % 3] - p2[diridx + 1 % 3]) >= sr) {
        return false;
      }
      if (std::abs(p1[diridx + 2 % 3] - p2[diridx + 2 % 3]) >= sr) {
        return false;
      }
      // Check distance
      auto distance = rti::util::distance<numeric_type>({p1[0], p1[1], p1[2]}, {p2[0], p2[1], p2[2]});
      if (distance < r1+r2) {
      // auto squrddist = rti::util::squrd_distance<numeric_type>(p1, p2);
      // if (squrddist < sr * sr) {
        return true;
      }
      return false;
    }
    
    void construct_neighborhood_naive_1(std::vector<rti::util::quadruple<numeric_type> > const& points)
    {
      std::cout << "points.size() == " << points.size() << std::endl;
      std::cout << "Starting the quadratic loop v1" << std::endl;
      //#pragma omp parallel for
      for (size_t idx1 = 0; idx1 < points.size(); ++idx1) {
        for (size_t idx2 = idx1; idx2 < points.size(); ++idx2) {
          auto& p1 = points[idx1];
          auto& r1 = p1[3];
          auto& p2 = points[idx2];
          auto& r2 = p2[3];
          auto sr = r1 + r2;
          // auto distance = rti::util::distance<numeric_type>({p1[0], p1[1], p1[2]}, {p2[0], p2[1], p2[2]});
          // if (distance < r1+r2) {
          auto squrddist = rti::util::squrd_distance(p1, p2);
          if (squrddist < sr * sr) {
            nbhd[idx1].push_back(idx2);
            nbhd[idx2].push_back(idx1);
          }
        }
      }
    }

    void construct_neighborhood_naive_2(std::vector<rti::util::quadruple<numeric_type> > const& points)
    {
      std::cout << "points.size() == " << points.size() << std::endl;
      std::cout << "Starting the quadratic loop v2" << std::endl;
      // #pragma omp parallel for
      for (size_t idx1 = 0; idx1 < points.size(); ++idx1) {
        for (size_t idx2 = idx1; idx2 < points.size(); ++idx2) {
          if ( check_dist(points, idx1, idx2, 0) ) {
            nbhd[idx1].push_back(idx2);
            nbhd[idx2].push_back(idx1);
          }
        }
      }
    }

  private:

    // std::vector<rti::util::quadruple<numeric_type> > points;
    std::vector<std::vector<size_t> > nbhd;
  };
}}
