#pragma once

#include <embree3/rtcore.h>

namespace rti { namespace trace {
  // wraping the Embree rtc context and the rti context into a C-style struct
  template<typename Ty>
  struct context_c_wrapper {
    // See https://www.embree.org/api.html#rtcinitintersectcontext
    // Data layout:
    // The first member HAS TO BE the RTC context, that is, context from the Embree library.
    RTCIntersectContext mRtcContext; // not a reference or pointer but full data stored here
    rti::trace::absc_context<Ty>& mAbscContext;
  };
}}
