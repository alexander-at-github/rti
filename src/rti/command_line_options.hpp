#pragma once

//
// THIS FILE IS DEPRECATED AND WILL BE REMOVED EVENTUALLY
//

#include <unordered_map>

#include "rti/enum_class_hash_function.hpp"
#include "rti/logger.hpp"

namespace rti {
  class command_line_options {
    // This is a singleton class
  public:
    static void init(int argc, char* argv []) {
      init_or_get(argc, argv);
    }
    // One needs to call init() before get_instance(), otherwise the command_line_options
    // instance does not contain any entries.
    static command_line_options& get_instance() {
      return init_or_get(0, nullptr);
    }
    // Delete default constructor, copy constructor, and copy assignment operator
    command_line_options() = delete;
    command_line_options(command_line_options const&) = delete;
    void operator=(command_line_options const&) = delete;

    enum class option_type {INFILE_NAME, MAX_THREADS};

    std::string get_option_value(rti::command_line_options::option_type pType) {
      for (auto & oo : options) {
        if (oo.first /* first is the key */ == pType) {
          return oo.second.value; // second is the value of the map
        }
      }
      return "";
    }
  private:
    struct option_spec {
      std::string name;
      std::string optStr;
      std::string value;
    };
    std::unordered_map<option_type, option_spec, rti::enum_class_hash_function> options {
      {option_type::INFILE_NAME, option_spec{"input file option", "--infile", ""}},
      {option_type::MAX_THREADS, option_spec {"maximum number of threads option", "--max-threads", ""}}
    };
    static command_line_options& init_or_get(int argc, char* argv []) {
      static command_line_options instance(argc, argv);
      return instance;
    }
    // Private constructor
    command_line_options(int argc, char *argv[]) {
      RLOG_DEBUG << "Reading command line" << std::endl;
      for (int idx = 0; idx < argc; ++idx) {
        RLOG_TRACE << "argv[" << idx << "] == " << argv[idx] << std::endl;
      }
      for (int idx = 0; idx < argc; ++idx) {
        for (auto & oo : options) {
          if (oo.second.optStr.compare(argv[idx]) == 0 && idx < (argc-1)) {
            oo.second.value = argv[idx+1];
            RLOG_DEBUG
              << "command line option " << oo.second.optStr << " " << oo.second.value << " found" << std::endl;
            // Break the inner loop
            break;
          }
        }
      }
    }
  };
} // namespace rti
