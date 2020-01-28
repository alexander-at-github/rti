#pragma once

#include <cmath>

#include "rti/ray/adaptive_origin.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class adaptive_rectangle_origin_z : public rti::ray::adaptive_origin<Ty> {
  private:
    void swap(Ty& e1, Ty& e2) const {
      Ty tmp = e2;
      e2 = e1;
      e1 = tmp;
    }

  public:
    adaptive_rectangle_origin_z(Ty pZval, rti::util::pair<Ty> pC1, rti::util::pair<Ty> pC2) :
      zval(pZval),
      mC1(pC1),
      mC2(pC2) {
      // Rearrange corners if needed such that each value in the pair mC1
      // is smaller or equal to the corresponding value in mC2.
      if (mC1[0] > mC2[0]) this->swap(mC1[0], mC2[0]);
      if (mC1[1] > mC2[1]) this->swap(mC1[1], mC2[1]);
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Condition on ordering of corner points");

      std::cerr << "WARNING: This class is not thread-safe." << std::endl;

      auto xdivisions = 10u;
      auto ydivisions = 10u;
      auto stratuminitweight = 1.0;
      stratumweights.reserve(xdivisions);
      for (size_t xidx = 0; xidx < xdivisions; ++xidx) {
        stratumweights[xidx].reserve(ydivisions);
        for (size_t yidx = 0; yidx < ydivisions; ++yidx) {
          stratumweights[xidx][yidx] = stratuminitweight;
        }
      }
    }

  private:
    size_t number_of_x_divisions() {
      return stratumweights.size();
    }

    size_t number_of_y_divisions() {
      return stratumweights[0].size();
    }

    Ty size_of_x_division() {
      assert(mC1[0] <= mC2[0] && "Condition on ordering of corner points");
      return (mC2[0] - mC1[0]) / number_of_y_divisions();
    }

    Ty size_of_y_division() {
      assert(mC1[1] <= mC2[1] && "Condition on ordering of corner points");
      return (mC2[1] - mC1[1]) / number_of_y_divisions();
    }

  public:
    rti::util::triple<Ty> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Invariant on ordering of corner points");
      auto r1 = pRng.get(pRngState);
      auto r2 = (Ty) pRng.get(pRngState);
      auto r3 = (Ty) pRng.get(pRngState);

      auto numstratums = number_of_x_divisions() * number_of_y_divisions();
      auto r1scaled = scale_to_sumstratumweights(r1, pRng.min(), pRng.max());
      auto xrsidx = 0u;
      auto yrsidx = 0u;
      auto weightsum = 0.0;
      while (weightsum < r1scaled) {
        weightsum += stratumweights[xrsidx][yrsidx];
        yrsidx += 1;
        if (yrsidx >= stratumweights[xrsidx].size()) {
          xrsidx += 1;
          yrsidx = 0;
        }
        assert(xrsidx < stratumweights.sizs() && yrsidx < stratumweights[xrsidx].size() && "Assumption");
      }
      assert(weightsum < sumstratumweights && "Assumption");
      assert(xrsidx < stratumweights.sizs() && yrsidx < stratumweights[xrsidx].size() && "Assumption");

      std::cout
        << "[DEBUG]"
        << " xrsidx == " << xrsidx
        << " yrsidx == " << yrsidx << std::endl;

      auto xStratumMin = xrsidx * size_of_x_division();
      auto xStratumMax = (xrsidx + 1) * size_of_x_division();
      auto yStratumMin = yrsidx * size_of_y_division();
      auto yStratumMax = (yrsidx + 1) * size_of_y_division();

      auto randRngMin = (Ty) pRng.min();
      auto randRngMax = (Ty) pRng.max();
      auto r2prct = (r2 - randRngMin) / (randRngMax - randRngMin);
      auto r3prct = (r3 - randRngMin) / (randRngMax - randRngMin);

      Ty xx = xStratumMin + (xStratumMax - xStratumMin) * r2prct;
      Ty yy = yStratumMin + (yStratumMax - yStratumMin) * r3prct;
      return {xx, yy, zval};
    }

  private:
    double scale_to_sumstratumweights(uint64_t random, uint64_t min, uint64_t max) {
      auto percent = (random - min) / (max - min);
      return percent * sumstratumweights;
    }

    rti::util::pair<size_t> get_stratum_of_coord(rti::util::triple<Ty> coord) {
      auto xx = coord[0];
      auto yy = coord[1];
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Condition on ordering of corner points");
      assert(xx <= mC2[0] && yy <= mC2[1] && "Assumption");
      auto xidx = std::floor(xx / size_of_x_division());
      auto yidx = std::floor(yy / size_of_y_division());
      return {(size_t) xidx, (size_t) yidx};
    }

  public:
    // When tracing a ray, one can pass the maximum of all relative errors
    void consider(rti::util::triple<Ty> xyz, double relativeerror) override final {
      assert(xyz[2] == zval && "Assumption");
      auto indices = get_stratum_of_coord(xyz);
      auto xi = indices[0];
      auto yi = indices[1];
      s1res[xi][yi] += relativeerror;
      s2res[xi][yi] += relativeerror * relativeerror;
      cnts[xi][yi] += 1;
    }

  private:
    double estimate_of_relative_error_on_surface_at_location(size_t xidx, size_t yidx) {
      return s1res[xidx][yidx] / cnts[xidx][yidx];
    }

    double relative_error_of_relative_error_at_location(size_t xidx, size_t yidx) {
      auto s1reSquared = s1res[xidx][yidx] * s1res[xidx][yidx];
      return std::sqrt(s2res[xidx][yidx] / s1reSquared - 1.0 / cnts[xidx][yidx]);
    }

    void update_resre() {
      sumresre = 0;
      for(size_t xidx = 0; xidx < s1res.size(); ++xidx) {
        for (size_t yidx = 0; yidx < s1res[xidx].size(); ++yidx) {
          resre[xidx][yidx] = relative_error_of_relative_error_at_location(xidx, yidx);
          sumresre += resre[xidx][yidx];
        }
      }
    }

    // One should call update_resre() before calling update_weights()
    void update_weights() {
      sumstratumweights = 0;
      for (size_t xidx = 0; xidx < stratumweights.size(); ++xidx) {
        for (size_t yidx = 0; yidx < stratumweights[xidx].size(); ++yidx) {
          stratumweights[xidx][yidx] =
            resre[xidx][yidx] <= relativeerrorthreshold ?
            estimate_of_relative_error_on_surface_at_location(xidx, yidx) :
            1.0; // full weight
          if (stratumweights[xidx][yidx] < weightminimum)
            stratumweights[xidx][yidx] = weightminimum;
          sumstratumweights += stratumweights[xidx][yidx];
        }
      }
    }
  public:
    void update_adaptive_sampling_state() override final {
      update_resre();
      update_weights();
    }

  private:
    Ty zval;
    rti::util::pair<Ty> mC1;
    rti::util::pair<Ty> mC2;

    using doublevector = std::vector<std::vector<double> >;

    // Thershold for the relative error of the relative error
    double relativeerrorthreshold = 0.02;
    // A minimum for the weight of the strata
    double weightminimum = 0.05;

    doublevector stratumweights;
    double sumstratumweights;

    doublevector s1res; // sum of relative errors from this stratum
    doublevector s2res; // sum of squared relative errors from this stratum
    doublevector cnts; // number of samples from this stratum

    doublevector resre; // relative errors of relative error
    double sumresre; // sum of relative errors of relative error
    // TODO: Do we use sumresre?

  };
}}
