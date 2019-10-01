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
      mS2s(pSize, 0) {}

    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const& pA) :
      mAcc(pA.mAcc), // copy construct the vector member
      mCnts(pA.mCnts),
      mS1s(pA.mS1s),
      mS2s(pA.mS2s) {}

    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const&& pA) :
      mAcc(std::move(pA.mAcc)), // move the vector member
      mCnts(std::move(pA.mCnts)),
      mS1s(std::move(pA.mS1s)),
      mS2s(std::move(pA.mS2s)) {}

    // A copy constructor which can accumulate values from two instances
    hit_accumulator_with_checks(hit_accumulator_with_checks<Ty> const& pA1,
                                hit_accumulator_with_checks<Ty> const& pA2) :
      // Precondition: the size of the accumulators are equal
      hit_accumulator_with_checks(pA1) { // copy construct from the first argument
      assert(pA1.mAcc.size() == pA2.mAcc.size() &&
             pA1.mCnts.size() == pA2.mCnts.size() &&
             pA1.mS1s.size() == pA2.mS1s.size() &&
             pA1.mS2s.size() == pA2.mS2s.size() && "Error: size missmatch");
      for (size_t idx = 0; idx < mAcc.size(); ++idx) {
        mAcc[idx] += pA2.mAcc[idx];
        mCnts[idx] += pA2.mCnts[idx];
        mS1s[idx] += pA2.mS1s[idx];
        mS2s[idx] += pA2.mS2s[idx];
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
    }

    std::vector<Ty> get_values() override final {
      return mAcc;
    }

    std::vector<size_t> get_cnts() override final {
      return mCnts;
    }

    std::vector<Ty> get_relative_error() override final {
      auto result = std::vector<Ty>(mS1s.size(), 0); // size of vector, all initial values are equal to zero
      for (size_t idx = 0; idx < result.size(); ++idx) {
        auto s1squared = mS1s[idx] * mS1s[idx];
        if (s1squared == 0) {
          result[idx] = std::numeric_limits<Ty>::max();
          continue;
        }
        assert(mCnts[idx] != 0 && "Assumption");
        // This is an approximation of the relative error assuming sqrt(N-1) =~ sqrt(N)
        // For details and an exact formula see the book Exploring Monte Carlo Methods by Dunn and Shultis
        // page 83 and 84.
        result[idx] = std::sqrt(mS2s[idx] / s1squared - 1 / mCnts[idx]);
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
    std::vector<Ty> mS1s;
    // S2 denotes the sum of squared sample values
    std::vector<Ty> mS2s;
  };
}}
