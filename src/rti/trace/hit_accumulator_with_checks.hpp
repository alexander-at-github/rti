#pragma once

#include "rti/trace/i_hit_accumulator.hpp"

namespace rti { namespace trace {
  template<typename Ty>
  class hit_accumulator_with_checks : public i_hit_accumulator<Ty> {
  public:
    // Constructors
    hit_accumulator_with_checks(size_t pSize) :
      mAcc(pSize, 0), // pSize number of elements initialized to 0.
      mCnts(pSize, 0),
      mS1s(pSize, 0),
      mS2s(pSize, 0),
      mS3s(pSize,0),
      mS4s(pSize,0) {}

    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const& pA) :
      mAcc(pA.mAcc), // copy construct the vector member
      mCnts(pA.mCnts),
      mS1s(pA.mS1s),
      mS2s(pA.mS2s),
      mS3s(pA.mS3s),
      mS4s(pA.mS4s) {}

    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const&& pA) :
      mAcc(std::move(pA.mAcc)), // move the vector member
      mCnts(std::move(pA.mCnts)),
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
    }

    // Assignment operators corresponding to the constructors
    hit_accumulator_with_checks<Ty>& operator=(hit_accumulator_with_checks<Ty> const& pOther) {
      if (this != &pOther) {
        // copy from pOther to this
        mAcc.clear();
        mAcc = pOther.mAcc;
        mCnts.clear();
        mCnts = pOther.mCnts;
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

    // Member Functions
    void use(unsigned int pPrimID, Ty value) override final {
      assert(pPrimID < mAcc.size() && "primitive ID is out of bounds");
      mAcc[pPrimID] += value;
      mCnts[pPrimID] += 1;

      mS1s[pPrimID] += value;
      mS2s[pPrimID] += value * value;
      mS3s[pPrimID] += value * value * value;
      mS4s[pPrimID] += value * value * value * value;
    }

    std::vector<Ty> get_values() override final {
      return mAcc;
    }

    std::vector<size_t> get_cnts() override final {
      return mCnts;
    }

    std::vector<Ty> get_relative_errors() override final {
      auto result = std::vector<Ty>(mS1s.size(), 0); // size of vector, all initial values are equal to zero
      for (size_t idx = 0; idx < result.size(); ++idx) {
        result[idx] = (Ty) get_relative_error_for_id(idx);
      }
      return result;
    }

  public:
    double get_relative_error_for_id(unsigned int pPrimID) override final {
      auto s1square = mS1s[pPrimID] * mS1s[pPrimID];
      if (s1square == 0) {
        return std::numeric_limits<double>::max();
      }
      // We require at least 2 samples to compute the relative error.
      if (mCnts[pPrimID] <= 1) {
        return std::numeric_limits<double>::max();
      }
      assert(mCnts[pPrimID] != 0 && "Assumption");
      // This is an approximation of the relative error assuming sqrt(N-1) =~ sqrt(N)
      // For details and an exact formula see the book Exploring Monte Carlo Methods by Dunn and Shultis
      // page 83 and 84.
      return (std::sqrt(mS2s[pPrimID] / s1square - 1 / mCnts[pPrimID]));
    }

    std::vector<Ty> get_vov() override final { // variance of variance
      auto result = std::vector<Ty>(mS1s.size(), 0); // size of vectore, all initial values are equal to zero
      for (size_t idx = 0; idx < result.size(); ++idx) {
        if (mCnts[idx] == 0) {
          result[idx] = std::numeric_limits<Ty>::max();
          continue;
        }
        // We require a minimum number of samples.
        // This also prevents numericals errors in the computation of the variance of the variance.
        // magic number
        if (mCnts[idx] <= 8) {
          result[idx] = std::numeric_limits<Ty>::max();
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
          result[idx] = std::numeric_limits<Ty>::max();
          continue;
        }
        //result[idx] = (Ty) (numerator / denominator);
        result[idx] = (Ty) (numerator / denominator2 * nn * nn);
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
        //             << "relative error == " << this->get_relative_errors()[idx] << std::endl;
        // }
      }
      return result;
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
    std::vector<Ty> mAcc;
    std::vector<size_t> mCnts;

    // actuall - for now - mAcc und mS1s do the same thing!
    // We might want to remove one of them later.

    // S1 denotes the sum of sample values
    std::vector<double> mS1s;
    // S2 denotes the sum of squared sample values
    std::vector<double> mS2s;
    // S3 denotes the sum of cubed sample values
    std::vector<double> mS3s;
    // S4 denotes the sum of the 4th-power of the sample values
    std::vector<double> mS4s;
  };
}}
