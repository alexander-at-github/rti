#pragma once

#include "rti/trace/i_hit_accumulator.hpp"

namespace rti { namespace trace {
  template<typename Ty>
  class hit_accumulator_with_checks : public rti::trace::i_hit_accumulator<Ty> {

    using internal_numeric_type = double;

  public:

    hit_accumulator_with_checks(size_t pSize) :
      mAcc(pSize, 0), // pSize number of elements initialized to 0.
      mCnts(pSize, 0),
      mTotalCnts(0),
      exposedareas(pSize, 0),
      mS1s(pSize, 0),
      mS2s(pSize, 0),
      mS3s(pSize,0),
      mS4s(pSize,0) {}

    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const& pA) :
      mAcc(pA.mAcc), // copy construct the vector member
      mCnts(pA.mCnts),
      mTotalCnts(pA.mTotalCnts),
      exposedareas(pA.exposedareas),
      mS1s(pA.mS1s),
      mS2s(pA.mS2s),
      mS3s(pA.mS3s),
      mS4s(pA.mS4s) {}

    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const&& pA) :
      mAcc(std::move(pA.mAcc)), // move the vector member
      mCnts(std::move(pA.mCnts)),
      mTotalCnts(std::move(pA.mTotalCnts)),
      exposedareas(std::move(exposedareas)),
      mS1s(std::move(pA.mS1s)),
      mS2s(std::move(pA.mS2s)),
      mS3s(std::move(pA.mS3s)),
      mS4s(std::move(pA.mS4s)) {}

    // A copy constructor which can accumulate values from two instances
    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const& pA1,
                                hit_accumulator_with_checks<Ty> const& pA2) :
      // Precondition: the size of the accumulators are equal
      hit_accumulator_with_checks(pA1) { // copy construct from the first argument
      assert(pA1.mAcc.size() == pA2.mAcc.size() &&
             pA1.mCnts.size() == pA2.mCnts.size() &&
             pA1.mS1s.size() == pA2.mS1s.size() &&
             pA1.mS2s.size() == pA2.mS2s.size() &&
             pA1.mS3s.size() == pA2.mS3s.size() &&
             pA1.mS4s.size() == pA2.mS4s.size() &&
             "Error: size missmatch");
      for (size_t idx = 0; idx < mAcc.size(); ++idx) {
        mAcc[idx] += pA2.mAcc[idx];
        mCnts[idx] += pA2.mCnts[idx];
        mS1s[idx] += pA2.mS1s[idx];
        mS2s[idx] += pA2.mS2s[idx];
        mS3s[idx] += pA2.mS3s[idx];
        mS4s[idx] += pA2.mS4s[idx];
      }
      mTotalCnts = pA1.mTotalCnts + pA2.mTotalCnts;
      /* Assertions about the exposed areas saved in the input instances */
      assert(pA1.exposedareas.size() == pA2.exposedareas.size());
      assert(exposedareas.size() == pA1.exposedareas.size());
      for (size_t idx = 0; idx < pA1.exposedareas.size(); ++idx) {
        assert(
          ( pA1.exposedareas[idx] == 0 ||
            pA2.exposedareas[idx] == 0 ||
            pA1.exposedareas[idx] == pA2.exposedareas[idx] )
          && "Correctness Assumption");
        exposedareas[idx] =
          pA1.exposedareas[idx] > pA2.exposedareas[idx] ?
          pA1.exposedareas[idx] : pA2.exposedareas[idx];
      }
    }

    // Assignment operators corresponding to the constructors
    hit_accumulator_with_checks<Ty>& operator=(hit_accumulator_with_checks<Ty> const& pOther) {
      if (this != &pOther) {
        // copy from pOther to this
        mAcc.clear();
        mAcc = pOther.mAcc;
        mCnts.clear();
        mCnts = pOther.mCnts;
        mTotalCnts = pOther.mTotalCnts;
        exposedareas.clear();
        exposedareas = pOther.exposedareas;
        mS1s.clear();
        mS1s = pOther.mS1s;
        mS2s.clear();
        mS2s = pOther.mS2s;
        mS3s.clear();
        mS3s = pOther.mS3s;
        mS4s.clear();
        mS4s = pOther.mS4s;
      }
      return *this;
    }

    hit_accumulator_with_checks<Ty>& operator=(hit_accumulator_with_checks<Ty> const&& pOther) {
      if (this != &pOther) {
        // move from pOther to this
        mAcc.clear();
        mAcc = std::move(pOther.mAcc);
        mCnts.clear();
        mCnts = std::move(pOther.mCnts);
        mTotalCnts = pOther.mTotalCnts;
        exposedareas.clear();
        exposedareas = std::move(pOther.exposedareas);
        mS1s.clear();
        mS1s = std::move(pOther.mS1s);
        mS2s.clear();
        mS2s = std::move(pOther.mS2s);
        mS3s.clear();
        mS3s = std::move(pOther.mS3s);
        mS4s.clear();
        mS4s = std::move(pOther.mS4s);
      }
      return *this;
    }


    void use(unsigned int pPrimID, Ty value) override final {
      assert(pPrimID < mAcc.size() && "primitive ID is out of bounds");
      mAcc[pPrimID] += (internal_numeric_type) value;
      mCnts[pPrimID] += 1;
      mTotalCnts += 1;

      mS1s[pPrimID] += (internal_numeric_type) value;
      mS2s[pPrimID] += (internal_numeric_type) value * value;
      mS3s[pPrimID] += (internal_numeric_type) value * value * value;
      mS4s[pPrimID] += (internal_numeric_type) value * value * value * value;
    }

    std::vector<internal_numeric_type> get_values() override final {
      return mAcc;
    }

    std::vector<size_t> get_cnts() override final {
      return mCnts;
    }

    size_t get_cnts_sum() override final {
      return mTotalCnts;
    }

    std::vector<internal_numeric_type> get_relative_error() override final {
      std::cerr << "### Fix! relative error" << std::endl;
      auto result = std::vector<internal_numeric_type>(mS1s.size(), 0); // size, initial values
      for (size_t idx = 0; idx < result.size(); ++idx) {
        auto s1square = mS1s[idx] * mS1s[idx];
        if (s1square == 0) {
          result[idx] = std::numeric_limits<internal_numeric_type>::max();
          continue;
        }
        // We require at least 2 samples to compute the relative error.
        if (mCnts[idx] <= 1) {
          result[idx] = std::numeric_limits<internal_numeric_type>::max();
          continue;
        }
        assert(mCnts[idx] != 0 && "Assumption");
        // This is an approximation of the relative error assuming sqrt(N-1) =~ sqrt(N)
        // For details and an exact formula see the book Exploring Monte Carlo Methods by Dunn and Shultis
        // page 83 and 84.
        result[idx] = (internal_numeric_type) (std::sqrt(mS2s[idx] / s1square - 1 / mCnts[idx]));
        // Debug
        // if (result[idx] != std::numeric_limits<internal_numeric_type>::max()) {
        //   std::cerr << "mCnts[idx] == " << mCnts[idx] << std::endl;
        //   std::cerr << "s1square == " << s1square << std::endl;
        //   std::cerr << "mS2ss[idx] == " << mS2s[idx] << std::endl;
        //   std::cerr << "result[idx] == " << result[idx] << std::endl;
        // }
      }
      return result;
    }

    std::vector<internal_numeric_type> get_vov() override final { // variance of variance
      std::cerr << "### Fix! vov" << std::endl;
      auto result = std::vector<internal_numeric_type>(mS1s.size(), 0); // size, initial values
      for (size_t idx = 0; idx < result.size(); ++idx) {
        if (mCnts[idx] == 0) {
          result[idx] = std::numeric_limits<internal_numeric_type>::max();
          continue;
        }
        // We require a minimum number of samples.
        // This also prevents numericals errors in the computation of the variance of the variance.
        // magic number
        if (mCnts[idx] <= 8) {
          result[idx] = std::numeric_limits<internal_numeric_type>::max();
          continue;
        }
        auto s1 = mS1s[idx];
        auto s2 = mS2s[idx];
        auto s3 = mS3s[idx];
        auto s4 = mS4s[idx];
        auto s1square = s1 * s1;
        auto s1FourthPow = s1square * s1square;
        auto s2square = s2 * s2;
        auto nn = mCnts[idx];
        auto nnsquare = nn * nn;
        auto nncube = nn * nn * nn;
        // For details and explanation of the formula see the book Exploring Monte Carlo Methods
        // by Dunn and Shultis page 84 and 85.
        auto numerator = (s4) - (4 * s1 * s3 / nn) + (8 * s2 * s1square / nnsquare)
                         - (4 * s1FourthPow / nncube) - (s2square / nn);
        auto denomroot = s2 - s1square / nn;
        auto denominator = denomroot * denomroot;
        auto denomroot2 = nn * s2 - s1square;
        auto denominator2 = denomroot2 * denomroot2;
        if (denominator == 0) {
          result[idx] = std::numeric_limits<internal_numeric_type>::max();
          continue;
        }
        //result[idx] = (internal_numeric_type) (numerator / denominator);
        result[idx] = (internal_numeric_type) (numerator / denominator2 * nn * nn);
        // Debug
        // if (result [idx] >= 100000) {
        //   std::cerr << std::setprecision(10);
        //   std::cerr << "vov >= 100000 holds" << std::endl;
        //   std::cerr << "nn == " << nn << "    "
        //             << "s1 == " << s1 << "    "
        //             << "s1square == " << s1square << "    "
        //             << "s2 == " << s2 << "    "
        //             << "s3 == " << s3 << "    "
        //             << "s4 == " << s4 << std::endl;
        //   std::cerr << "numerator == " << numerator << "    " << "denominator == " << denominator << "    "
        //             << "denomroot == " << denomroot << std::endl;
        //   std::cerr << "denominator2 == " << denominator2 << "    "
        //             << "denomroot2 == " << denomroot2 << std::endl;
        //   std::cerr << "s1square / nn == " << s1square / nn << std::endl;
        //   std::cerr << "vov == " << result[idx] << "    "
        //             << "relative error == " << this->get_relative_error()[idx] << std::endl;
        // }
      }
      return result;
    }

    void set_exposed_areas(std::vector<internal_numeric_type>& areas) override final
    {
      exposedareas = areas;
    }

    std::vector<internal_numeric_type> get_exposed_areas() override final
    {
      return exposedareas;
    }

    void print(std::ostream& pOs) const override final {
      pOs << "(";
      auto const* separator = " ";
      auto const* sep       = "";
      for (auto& vv : mAcc) {
        pOs << sep << vv;
        sep = separator;
      }
      pOs << ")" << std::endl;
    }
  private:
    std::vector<internal_numeric_type> mAcc;
    std::vector<size_t> mCnts;
    size_t mTotalCnts;
    std::vector<internal_numeric_type> exposedareas;

    // actuall - for now - mAcc und mS1s do the same thing!
    // We might want to remove one of them later.

    // S1 denotes the sum of sample values
    std::vector<internal_numeric_type> mS1s;
    // S2 denotes the sum of squared sample values
    std::vector<internal_numeric_type> mS2s;
    // S3 denotes the sum of cubed sample values
    std::vector<internal_numeric_type> mS3s;
    // S4 denotes the sum of the 4th-power of the sample values
    std::vector<internal_numeric_type> mS4s;
  };
}}
