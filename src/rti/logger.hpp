#pragma once

#include <iostream>

// An effective logger

namespace rti {
// Set the log level by commenting the following lines out
//
//#define RLOG_LEVEL_TRACE
//#define RLOG_LEVEL_DEBUG
#define RLOG_LEVEL_INFO
//#define RLOG_LEVEL_WARNING
#define RLOG_LEVEL_ERROR

  // Set the logging destination here
  static std::ostream& mLStream = std::cerr;

#ifdef RLOG_LEVEL_TRACE
#define RLOG_TRACE rti::mLStream
#else
#define RLOG_TRACE if (true) {} else rti::mLStream
#endif

#ifdef RLOG_LEVEL_DEBUG
#define RLOG_DEBUG rti::mLStream
#else
#define RLOG_DEBUG if (true) {} else rti::mLStream
#endif

#ifdef RLOG_LEVEL_INFO
#define RLOG_INFO rti::mLStream
#else
#define RLOG_INFO if (true) {} else rti::mLStream
#endif

#ifdef RLOG_LEVEL_WARNING
#define RLOG_WARNING rti::mLStream
#else
#define RLOG_WARNING if (true) {} else rti::mLStream
#endif

#ifdef RLOG_LEVEL_ERROR
#define RLOG_ERROR rti::mLStream
#else
#define RLOG_ERROR if (true) {} else rti::mLStream
#endif

} // namespace rti
