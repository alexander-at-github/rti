#pragma once

#include <chrono>

namespace rti { namespace util {
  class timer {
  public:
    timer() :
      mStartTime(time_stamp_now())
    {}

    void restart() {
      mStartTime = time_stamp_now();
    }

    double elapsed_secounds() const {
      return double(time_stamp_now() - mStartTime) * 1e-9;
    }

    std::uint64_t elapsed_nanoseconds() const {
      return time_stamp_now() - mStartTime;
    }

  private:
    static std::uint64_t time_stamp_now() {
      return std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::steady_clock::now().time_since_epoch()).count();
    }

    std::uint64_t mStartTime;

  };
}} // namespace
