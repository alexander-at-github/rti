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
      boost::log::core::get()->set_filter(blt::severity >= blt::trace);
    }
  }
}
