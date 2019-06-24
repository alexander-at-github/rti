#pragma once

#include <ostream>

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
      bucket_counter(pOther.mLength, pOther.mNumBuckets) {
      // Copy content to new object
      for (size_t idx = 0; idx < mCounts.size(); ++idx) {
        mCounts[idx] = pOther.mCounts[idx];
      }
    }

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
      //std::cout << "[bucket_counter::use()] mCounts[bucket] ==" << mCounts[bucket] << std::endl;
    }

    void print(std::ostream& pOs) const {
      print_normalized(pOs);
      //print_percentage(pOs);
    }

    void print_percentage(std::ostream& pOs) const {
      size_t totalCount = 0;
      for (auto const& count : mCounts) {
        totalCount += count;
      }
      for (size_t idx = 0; idx < mCounts.size(); ++idx) {
        // One needs to make sure that there is no division by 0.
        double percentageCount =
          totalCount == 0 ?
          0 :
          ((double) mCounts[idx]) / totalCount;
        pOs << (((double) idx) / (mCounts.size() - 1)) << " " << percentageCount << std::endl;
      }
    }

    void print_normalized(std::ostream& pOs) const {
      size_t maxCount = *std::max_element(mCounts.begin(), mCounts.end());
      for (size_t idx = 0; idx < mCounts.size(); ++idx) {
        // One needs to make sure that there is no division by 0.
        double normalizedCount =
          maxCount == 0 ?
          0 :
          ((double) mCounts[idx]) / maxCount;
        pOs << (((double) idx) / (mCounts.size() - 1)) << " " << normalizedCount << std::endl;
      }
    }

    static bucket_counter combine(bucket_counter pBc1, bucket_counter pBc2) {
      // Copy constructor copies entries to new object
      bucket_counter newBc(pBc1);
      for(size_t idx = 0; idx < newBc.mCounts.size(); ++idx) {
        newBc.mCounts[idx] += pBc2.mCounts[idx];
      }
      return newBc;
    }

  private:
    double mLength;
    size_t mNumBuckets;
    double mSlap;
    std::vector<size_t> mCounts;
  };
} // namespace rti
