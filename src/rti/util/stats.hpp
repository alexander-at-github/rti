#pragma once

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <vector>

#include "utils.hpp"
#include "../trace/i_hit_accumulator.hpp"

namespace rti::util {
  template<typename numeric_type>
  class stats {
    
  public:

    stats(size_t size1, size_t size2) :
      s1(size1),
      s2(size2) {}
    
    // void record(util::triple<numeric_type> point)
    // {
    //   points.push_back(point);
    // }

    // using internal_type = double;
    // void kolmogorov_smirnov_test_uniform()
    // {
    //   auto nn = points.size();
    //   auto dd = (internal_type) 0.0;
    //   for (size_t idx = 0; idx < nn; ++idx) {
    //     auto& pp = points[idx];
    //     auto& xx = pp[0];
    //     auto& yy = pp[1];
    //     auto dq = get_data_quadrants(xx, yy);
    //     auto mq = get_uniform_model_quadrants(xx, yy);
    //     // Corrections
    //     // if (dq[0] > mq[0]) { dq[0] += 1/nn; }
    //     // if (dq[1] > mq[1]) { dq[1] += 1/nn; }
    //     // if (dq[2] > mq[2]) { dq[2] += 1/nn; }
    //     // if (dq[3] > mq[3]) { dq[3] += 1/nn; }
    //     using std::abs;
    //     dd = std::max({dd, abs(dq[0]-mq[0]), abs(dq[1]-mq[1]), abs(dq[2]-mq[2]), abs(dq[3]-mq[3])});
    //   }
    //   std::cout << "K-S: dd == " << dd << std::endl;

    //   compute_data_quadrants();
    // }

    // // Given an origin (xx, yy) and the instance variable points this method counts how
    // // of themare in each quadrant aroung the origin, and returns the normalized fractions.
    // // Quadrants are labeled alphabetically, counterclockwise from the upper right.
    // std::array<internal_type, 4> get_data_quadrants(numeric_type& xx, numeric_type& yy)
    // {
    //   return {0, 0, 0, 0};
    //   auto na = 0u;
    //   auto nb = 0u;
    //   auto nc = 0u;
    //   auto nd = 0u;
    //   for (auto const & pp : points) {
    //     auto const& px = pp[0];
    //     auto const& py = pp[1];
    //     if (px == xx && py ==yy) { continue; }
    //     if (py > yy) { px > xx ? ++na : ++nb; }
    //     else { px > xx ? ++nd : ++nc; }
    //   }
    //   using t = internal_type;
    //   return {(t) na, (t) nb, (t) nc, (t) nd};
    // }

    // void compute_data_quadrants()
    // {
    //   using nt = numeric_type;
    //   auto pXid = std::vector<size_t> {};
    //   auto pYid = std::vector<size_t> {};
    //   //pXid.reserve(points.size());
    //   pXid.resize(points.size());
    //   std::iota(pXid.begin(), pXid.end(), 0u); // fill vector with 0, 1, ...
    //   // pYid.reserve(points.size());
    //   pYid.resize(points.size());
    //   std::iota(pYid.begin(), pYid.end(), 0u);
    //   // sort according to xx values
    //   auto& ps = points; // local variable referencing points for lambda capture.
    //   std::sort(std::execution::par_unseq, pXid.begin(), pXid.end(),
    //             [&ps](size_t const& i1, size_t const& i2) { return ps[i1][0] < ps[i2][0]; });
    //   // sort accoring to yy values
    //   std::sort(std::execution::par_unseq, pYid.begin(), pYid.end(),
    //             [&ps](size_t const& i1, size_t const& i2) { return ps[i1][1] < ps[i2][1]; });
    //   std::cout << "HERE" << std::endl;
    //   std::cout << "Xid[0] == " << pXid[0] << std::endl;
    //   std::cout << "points[pXid[0]][0] == " << points[pXid[0]][0] << " " << points[pXid[0]][1] << " " << points[pXid[0]][2] << std::endl;
    //   std::cout << "points[pXid[1]][1] == " << points[pXid[1]][0] << " " << points[pXid[1]][1] << " " << points[pXid[1]][2] << std::endl;

    //   auto pXidInv = pXid; // TODO !
    //   auto pYidInv = pYid; // TODO !
      
    //   auto result = std::vector<std::array<internal_type, 4> > {};
    //   auto na = std::vector<size_t> {};
    //   auto nb = std::vector<size_t> {};
    //   auto nc = std::vector<size_t> {};
    //   auto nd = std::vector<size_t> {};
    //   na.resize(points.size(), 0u);
    //   nb.resize(points.size(), 0u);
    //   nc.resize(points.size(), 0u);
    //   nd.resize(points.size(), 0u);
    //   result.reserve(points.size());
    //   auto nleft = 0u;
    //   for (size_t idx = 0; idx < points.size(); ++idx) {
    //     auto xordpos = pXidInv[idx];
    //     auto yordpos = pYidInv[idx];
    //     result[idx] = {0, 0, 0, 0};
    //   }
    //   auto yidacc = std::vector<size_t> {};
    //   for (auto const& idx : pXid) {
    //     auto pair = get_above_below(yidacc, pXid);
    //     nb[idx] = pair.first;
    //     nc[idx] = pair.second;
        
    //   }
    // }

    // std::array<internal_type, 4> get_uniform_model_quadrants(numeric_type& xx, numeric_type& yy)
    // {
    //   return {0, 0, 0, 0};
    // }

    void chi_squared_test_uniform(trace::i_hit_accumulator<numeric_type>& acc)
    {
      auto cnts = acc.get_cnts();

      if (cnts.size() != s1 * s2) {
        std::cerr << "stats::get_expected() Error: argument of wrong size" << std::endl;
      }
      constexpr auto numCornerNodes = 4u;
      auto numEdgeNodes = (s1 - 2)*2 + (s2 - 2)*2;
      auto numInnerNodes = cnts.size() - numCornerNodes - numEdgeNodes;
      if (numEdgeNodes != 84) {
        std::cerr << "ERROR" << std::endl;
      }
      sort(cnts.begin(), cnts.end());
      auto corners = std::vector<size_t> (cnts.begin(), cnts.begin() + numCornerNodes);
      // auto corners = std::array<size_t, numCornerNodes> {};
      // std::copy_n(cornersvec.begin(), numCornerNodes, corners.begin());
      auto edges = std::vector<size_t> (cnts.begin() + numCornerNodes, cnts.begin() + numCornerNodes + numEdgeNodes);
      auto inner = std::vector<size_t> (cnts.begin() + numCornerNodes + numEdgeNodes, cnts.end());
      std::cout << "Corner values in [" << corners.front() << ", " << corners.back() << "]" << std::endl;
      std::cout << "Edge values in [" << edges.front() << ", " << edges.back() << "]" << std::endl;
      std::cout << "Inner values in [" << inner.front() << ", " << inner.back() << "]" << std::endl;

      using vv = std::vector<size_t>;
      auto totalcnts = std::accumulate<vv::const_iterator, vv::value_type>(
          cnts.cbegin(), cnts.cend(), 0u);

      auto correctcornerval =
          (double)totalcnts /
          (4 * numInnerNodes + 2 * numEdgeNodes + numCornerNodes);
      auto correctedgeval = 2 * correctcornerval;
      auto correctinnerval = 4 * correctcornerval;
      std::cout << "correctcornerval == " << correctcornerval
                << " correctedgeval == " << correctedgeval
                << " correctinnerval == " << correctinnerval << std::endl;

      auto chsq = 0.0;
      auto degfreedom = cnts.size() - 1; // -1 ; see book Numerical Recipes page 732
      auto handlevalues = [&chsq, &degfreedom](auto const& actualVal, auto const& expectedVal) {
        if (expectedVal < 0 || (expectedVal == 0 && actualVal > 0)) {
          assert(false && "Error in Chi-Squared computation");
        }
        if (expectedVal == 0 && actualVal == 0) {
          degfreedom -= 1;
          return;
        }
        auto dd = (double) expectedVal - actualVal;
        auto contribution = dd*dd/expectedVal;
        // if (contribution > 3) {
        //   std::cout << "chi-squared contribution == " << contribution << std::endl;
        // }
        chsq += contribution;
      };

      for (auto const& pair : std::vector<std::pair<std::vector<size_t>, double> > {{corners, correctcornerval}, {edges, correctedgeval}, {inner, correctinnerval}}) {
        auto const& expectedVal = pair.second;
        for (auto const& val : pair.first) {
          std::cout << val << "\t\t\t" << expectedVal << std::endl;
          handlevalues(val, expectedVal);
        }
      }

      // std::cout << "==========" << std::endl;
      // for (auto const& pair : std::vector<std::pair<std::vector<size_t>, double> > {{corners, correctcornerval}, {edges, correctedgeval}, {inner, correctinnerval}}) {
      //   for (auto const& val : pair.first) {
      //     std::cout << val << std::endl;
      //   }
      // }
      // std::cout << "==========" << std::endl;
      // for (auto const& pair : std::vector<std::pair<std::vector<size_t>, double> > {{corners, correctcornerval}, {edges, correctedgeval}, {inner, correctinnerval}}) {
      //   for (auto const& val : pair.first) {
      //     std::cout << pair.second << std::endl;
      //   }
      // }
      // std::cout << "==========" << std::endl;
      
      std::cout << "Chi-Squared == " << chsq << std::endl;
      std::cout << "degrees of freedom == " << degfreedom << std::endl;
    }
    
  private:
    size_t s1;
    size_t s2;
    std::vector<util::triple<numeric_type> > points;
  };
}
