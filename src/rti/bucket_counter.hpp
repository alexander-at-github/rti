#pragma once

#include "rti/i_intersection_handler.hpp"

namespace rti {
  class bucket_counter : public i_intersection_handler {
  public:

    bucket_counter(double pLength, size_t pNumBuckets) :
      mLength(pLength),
      mNumBuckets(pNumBuckets),
      mSlap(pLength / pNumBuckets),
      mCounts(pNumBuckets) {
      for (auto& count : mCounts) {
        count = 0;
      }
    }

    // Copy constructor (for tbb::enumerable_thread_specific)
    bucket_counter(const bucket_counter& pOther) :
      bucket_counter(pOther.mLength, pOther.mNumBuckets) {}

    bucket_counter(std::vector<bucket_counter> pOthers) :
      bucket_counter(pOthers[0].mLength, pOthers[0].mNumBuckets) {
      for (auto other : pOthers) {
        if (other.mLength != this->mLength || other.mNumBuckets != this->mNumBuckets ||
            other.mSlap != this->mSlap || other.mCounts.size() != this->mCounts.size()) {
          assert(false && "Violation of preconditions");
          return;
        }
      }
      for (auto other : pOthers) {
        for (size_t idx = 0; idx < mCounts.size(); ++idx) {
          mCounts[idx] += other.mCounts[idx];
        }
      }
    }

    void use(const RTCRayHit& pRayHit) override final {
      //double hit_x = pRayHit.ray.org_x + pRayHit.ray.dir_x * pRayHit.ray.tfar;
      double hit_offset_x = pRayHit.ray.dir_x * pRayHit.ray.tfar;
      //RLOG_DEBUG << "hit_offset_x == " << hit_offset_x << std::endl;
      assert(0 <= hit_offset_x && hit_offset_x <= mLength && "Hit-x-value not within expected range");

      double rr = hit_offset_x / mSlap;
      assert(0 <= rr && "Error in bucket calculation");
      size_t bucket = (size_t) rr; // floor(rr)
      assert(bucket < mCounts.size() && "Error in bucket calculation");
      mCounts[bucket] += 1;
    }

  private:
    double mLength;
    size_t mNumBuckets;
    double mSlap;
    std::vector<size_t> mCounts;
  };
} // namespace rti
