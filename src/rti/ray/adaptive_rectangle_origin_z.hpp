#pragma once

#include <cmath>

#include "rti/ray/i_adaptive_origin.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class adaptive_rectangle_origin_z : public rti::ray::i_adaptive_origin<Ty> {
  private:
    using numeric_type = double;
    using weightvector = std::vector<numeric_type>;
    using doubleweightvector = std::vector<std::vector<double> >;

    void swap(Ty& e1, Ty& e2) {
      Ty tmp = e2;
      e2 = e1;
      e1 = tmp;
    }

    template<typename localType> void
    init_double_vector(std::vector<std::vector<localType> >& vec, size_t size1, size_t size2, localType value) {
      vec.clear();
      vec.resize(size1);
      for (size_t xidx = 0; xidx < size1; ++xidx) {
        vec[xidx].clear();
        vec[xidx].resize(size2, value);
      }
    }

    void update_sumstratumweights_variable() {
      sumstratumweights = 0.0;
      for(auto const& vec : stratumweights) {
        for (auto const& val : vec) {
          sumstratumweights += val;
        }
      }
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
      init_double_vector(stratumweights, xdivisions, ydivisions, stratuminitweight);
      update_sumstratumweights_variable();
      init_double_vector(s1res, xdivisions, ydivisions, 0.0);
      init_double_vector(s2res, xdivisions, ydivisions, 0.0);
      init_double_vector(cnts, xdivisions, ydivisions, (size_t) 0);
      init_double_vector(resre, xdivisions, ydivisions, 0.0);
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

    rti::util::pair<size_t>
    get_random_stratum_indices_according_to_weights(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState)
    {
      auto r1 = pRng.get(pRngState);
      auto r1scaled = scale_to_sumstratumweights(r1, pRng.min(), pRng.max());
      auto xrsidx = 0ul;
      auto yrsidx = 0ul;
      auto weightsumlow = 0.0;
      auto weightsumhigh = stratumweights[0][0];
      while ( ! (weightsumlow < r1scaled && r1scaled <= weightsumhigh) ) {
        weightsumlow = weightsumhigh;
        weightsumhigh += stratumweights[xrsidx][yrsidx];
        yrsidx += 1;
        if (yrsidx >= stratumweights[xrsidx].size()) {
          xrsidx += 1;
          yrsidx = 0;
        }
        assert(xrsidx < stratumweights.size() && yrsidx < stratumweights[xrsidx].size() && "Assumption");
      }
      assert(weightsumlow <= sumstratumweights && weightsumhigh <= sumstratumweights && "Assumption");
      assert(xrsidx < stratumweights.size() && yrsidx < stratumweights[xrsidx].size() && "Assumption");
      return {xrsidx, yrsidx};
    }

    rti::util::quadruple<Ty>
    get_stratum_mins_and_maxs(size_t xidx, size_t yidx) {
      auto xmin = mC1[0] + xidx * size_of_x_division();
      auto xmax = mC1[0] + (xidx + 1) * size_of_x_division();
      auto ymin = mC1[1] + yidx * size_of_y_division();
      auto ymax = mC1[1] + (yidx + 1) * size_of_y_division();
      return {xmin, xmax, ymin, ymax};
    }

  public:
    rti::util::triple<Ty> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Invariant on ordering of corner points");
      auto xIdxYidx = get_random_stratum_indices_according_to_weights(pRng, pRngState);
      auto xrsidx = xIdxYidx[0];
      auto yrsidx = xIdxYidx[1];

      auto r2 = (Ty) pRng.get(pRngState);
      auto r3 = (Ty) pRng.get(pRngState);

      // std::cout
      //   << "[DEBUG]"
      //   << " xrsidx == " << xrsidx
      //   << " yrsidx == " << yrsidx << std::endl;

      auto stratumMinsAndMaxs = get_stratum_mins_and_maxs(xrsidx, yrsidx);
      auto xStratumMin = stratumMinsAndMaxs[0];
      auto xStratumMax = stratumMinsAndMaxs[1];
      auto yStratumMin = stratumMinsAndMaxs[2];
      auto yStratumMax = stratumMinsAndMaxs[3];

      auto randRngMin = (Ty) pRng.min();
      auto randRngMax = (Ty) pRng.max();
      auto r2prct = (r2 - randRngMin) / (randRngMax - randRngMin);
      auto r3prct = (r3 - randRngMin) / (randRngMax - randRngMin);

      Ty xx = xStratumMin + (xStratumMax - xStratumMin) * r2prct;
      Ty yy = yStratumMin + (yStratumMax - yStratumMin) * r3prct;
      assert(mC1[0] <= xx && xx <= mC2[0] && "Assumption");
      assert(mC1[1] <= yy && yy <= mC2[1] && "Assumption");
      return {xx, yy, zval};
    }

  private:
    double scale_to_sumstratumweights(uint64_t random, uint64_t min, uint64_t max) {
      auto percent = (((double)random - min) / (max - min));
      return percent * sumstratumweights;
    }

    rti::util::pair<size_t> get_stratum_of_coord(rti::util::triple<Ty> coord) {
      auto xx = coord[0];
      auto yy = coord[1];
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Condition on ordering of corner points");
      assert(mC1[0] <= xx && xx <= mC2[0] && "Correctness Assertion");
      assert(mC1[1] <= yy && yy <= mC2[1] && "Correctness Assertion");
      auto xidx = std::floor((xx - mC1[0]) / size_of_x_division());
      auto yidx = std::floor((yy - mC1[1]) / size_of_y_division());
      assert(0 <= xidx && xidx < number_of_x_divisions() && "Correctness Assertion");
      assert(0 <= yidx && yidx < number_of_y_divisions() && "Correctness Assertion");
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
      std::cout << "[update_adaptive_sampling_state()]" << std::endl;
      update_resre();
      update_weights();
    }

  private:
    Ty zval;
    rti::util::pair<Ty> mC1;
    rti::util::pair<Ty> mC2;

    // Thershold for the relative error of the relative error
    double relativeerrorthreshold = 0.02;
    // A minimum for the weight of the strata
    double weightminimum = 0.05;

    doubleweightvector stratumweights;
    double sumstratumweights;

    doubleweightvector s1res; // sum of relative errors from this stratum
    doubleweightvector s2res; // sum of squared relative errors from this stratum
    std::vector<std::vector<size_t> > cnts; // number of samples from this stratum

    doubleweightvector resre; // relative errors of relative error
    double sumresre; // sum of relative errors of relative error
    // TODO: Do we use sumresre?

  };
}}
