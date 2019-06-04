#pragma once

//#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
//#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>

namespace blt = boost::log::trivial;
namespace bls = boost::log::sources;

namespace rti {
  // Create global logger
  bls::severity_logger<blt::severity_level> mRLogger;
  namespace logger {
    void init() {
      // Log to the console
      boost::log::add_console_log(std::cout);
      // optional: set boost log format

      // Set logging severity level
      boost::log::core::get()->set_filter(blt::severity >= blt::warning);
    }
  }
}


// CONTINUE: WRITE NEW LOGGER WITH COMPILE TIME FEATURE

// #include <iostream>

// namespace rti {
//   class logger_t {
//   public:
//     logger_t(std::ostream& pOStream) : mOutStream(pOStream) {};
//     enum class log_level_t {trace = 0, debug, warning, error, fatal, none};
//     std::ostream& get_outStream() {
//       return mOutStream;
//     }
//     void log(log_level_t pLevel) {
//       if (this->sLogLevel <= pLevel) {
//         return this->get_outStream();
//       }
//       return ... // TODO
//     }
//   private:
//     std::ostream& mOutStream;
//     static constexpr log_level_t sLogLevel = log_level_t::warning;
//   };
//   logger_t mLogger(std::cout);
// }


