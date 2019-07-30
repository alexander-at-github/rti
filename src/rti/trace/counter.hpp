#pragma once

#include "rti/trace/i_hit_counter.hpp"

namespace rti { namespace trace {
  class counter : public rti::trace::i_hit_counter {
  public:

    //counter() {}

    counter(size_t pCSize) :
      mCnts(pCSize, 0) {} // pCSize number of elements initialized to 0.

    counter(counter const& pDc) :
      mCnts(pDc.mCnts) {} // copy construct the vector member

    counter(counter const&& pDc) :
      mCnts(std::move(pDc.mCnts)) {} // move the vector member

    counter& operator=(counter const& pOther) {
      if (this != &pOther) {
        // copy stuff from pOther to this
        mCnts.clear();
        mCnts = pOther.mCnts;
      }
      return *this;
    }

    counter& operator=(counter const&& pOther) {
      if (this != &pOther) {
        // move stuff from pOther to this
        mCnts.clear();
        mCnts = std::move(pOther.mCnts);
      }
      return *this;
    }

    counter(counter const& pDc1, counter const& pDc2) :
      // Precondition: the sizes of the two parameters are equal
      counter(pDc1) {
      assert(pDc1.mCnts.size() == pDc2.mCnts.size() && "Size missmatch");
      for (size_t idx = 0; idx < mCnts.size(); ++idx) {
        mCnts[idx] += pDc2.mCnts[idx];
      }
    }

    void use(const RTCRayHit& pRayhit) override final {
      // Precondition: the primitive ID of the hit needs to be within the range of
      // the mCnts vector (which is initialized by the constructor).
      assert (pRayhit.hit.primID < mCnts.size() && "ray-hit primitive ID out of bounds");
      mCnts[pRayhit.hit.primID] += 1;
    }

    std::vector<size_t> get_counts() override final {
      return mCnts;
    }

    void print(std::ostream& pOs) const override final {
      pOs << "(";
      for (size_t idx = 0; idx < mCnts.size(); ++idx) {
        pOs << mCnts[idx];
        if (idx < mCnts.size() - 1) {
          pOs << " ";
        }
      }
      pOs << ")" << std::endl;
    }

  private:
    std::vector<size_t> mCnts;
  };
}} // namespace
