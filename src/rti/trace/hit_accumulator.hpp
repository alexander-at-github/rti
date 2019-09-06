#pragma once

#include "rti/trace/i_hit_accumulator.hpp"

namespace rti { namespace trace {
  template<typename Ty>
  class hit_accumulator : public i_hit_accumulator<Ty> {
  public:
    // Constructors
    hit_accumulator(size_t pSize) :
      mAcc(pSize, 0) {} // pSize number of elements initialized to 0.

    hit_accumulator(hit_accumulator<Ty> const& pA) :
      mAcc(pA.mAcc) {} // copy construct the vector member

    hit_accumulator(hit_accumulator<Ty> const&& pA) :
      mAcc(std::move(pA.mAcc)) {} // move the vector member

    // A copy constructor which can accumulate values from two instances
    hit_accumulator(hit_accumulator<Ty> const& pA1, hit_accumulator<Ty> const& pA2) :
      // Precondition: the size of the accumulators are equal
      hit_accumulator(pA1) { // copy construct from the first argument
      assert(pA1.mAcc.size() == pA2.mAcc.size() && "Error: size missmatch");
      for (size_t idx = 0; idx < mAcc.size(); ++idx) {
        mAcc[idx] += pA2.mAcc[idx];
      }
    }

    // Assignment operators corresponding to the constructors
    hit_accumulator<Ty>& operator=(hit_accumulator<Ty> const& pOther) {
      if (this != &pOther) {
        // copy from pOther to this
        mAcc.clear();
        mAcc = pOther.mAcc;
      }
      return *this;
    }

    hit_accumulator<Ty>& operator=(hit_accumulator<Ty> const&& pOther) {
      if (this != &pOther) {
        // move from pOther to this
        mAcc.clear();
        mAcc = std::move(pOther.mAcc);
      }
      return *this;
    }

    // Member Functions
    void use(unsigned int pPrimID, Ty value) override final {
      assert(pPrimID < mAcc.size() && "primitive ID is out of bounds");
      mAcc[pPrimID] += value;
    }

    std::vector<Ty> get_counts() override final {
      return mAcc;
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
  };
}}
