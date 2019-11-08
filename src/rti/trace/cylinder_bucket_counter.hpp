#pragma once

#include "rti/trace/i_hit_accumulator.hpp"
#include "rti/trace/hit_accumulator.hpp"

namespace rti { namespace trace {
  template<typename Ty>
  class cylinder_bucket_counter_z : public i_hit_accumulator<Ty> {
  public:
    // Constructors
    cylinder_bucket_counter_z(size_t pNumBuckets, int pZmin, int pZmax, size_t pHitAccumulatorSize) :
      mZmin(pZmin),
      mZmax(pZmax),
      mHitAcc(pHitAccumulatorSize) {
      mBuckets = std::vector(pNumBuckets, 0); // number of elements and initial values equal to zero
    }
    cylinder_bucket_counter_z(cylinder_bucket_counter_z const& pA) :
      mZmin(pZmin),  // copy everything
      mZmax(pZmax),
      mHitAcc(pA.mHitAcc),
      mBuckets(pA.mBuckets) {}
    cylinder_bucket_counter_z(cylinder_bucket_counter_z const&& pA) :
      mZmin(pZmin),  // move everything
      mZmax(pZmax),
      mHitAcc(std::move(pA.mHitAcc)),
      mBuckets(std::move(pA.mBuckets)) {}
    cylinder_bucket_counter_z(cylinder_bucket_counter_z const& pA1, cylinder_bucket_counter_z const& pA2) :
      // precondition: the size of the accumulators and the size of the vectors are equal
      cylinder_bucket_counter_z(pA1), // copy construct from first argument
      mHitAcc(pA1.mHitAcc, pA2.mHitAcc) { // cobine hit accumulator values by using the hit accumulator constructor
      assert(pA1.mBuckets.size() == pA2.mBuckets.size() && "Error: size missmatch");
      for (size_t idx = 0; idx < mBuckets.size(); ++idx) {
        mBuckets[idx] += pA2.mBuckets[idx];
      }
    }
    // Assignment Operators
    cylinder_bucket_counter_z<Ty>& operator= (cylinder_bucket_counter_z<Ty> const& pOther) {
      if (this != &pOther) {
        // copy
        mZmin = pOther.mZmin;
        mZmax = pOther.mZmax;
        mHitAcc = pOther.mHitAcc; // Calls the copy assignment operator from hit_accumulator
        mBuckets.clear();
        mBuckets = pOther.mBuckets;
      }
      return *this;
    }
    cylinder_bucket_counter_z<Ty>& operator= (cylinder_bucket_counter_z<Ty> const&& pOther) {
      if (this != &pOther) {
        // move
        mZmin = pOther.mZmin;
        mZmax = pOther.mZmax;
        mHitAcc = std::move(pOther.mHitAcc);
        mBuckets.clear();
        mBuckets = std::move(pOther.mBuckets);
      }
      return *this;
    }
    // Member Functions
    //void use(unsigned int pPrimID, Ty value) override final {} // TODO
  private:
    int mZmin = 0;
    int mZmax = 0;
    // holds a normal hit accumulator
    rti::trace::hit_accumulator mHitAcc;
    std::vector<Ty> mBuckets;
  };
}} // namespace
