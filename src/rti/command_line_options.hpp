#pragma once

#include <unordered_map>

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

    enum class option_type {MESH_FILE, MAX_THREADS};

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
    std::unordered_map<option_type, option_spec> options {
      {option_type::MESH_FILE, option_spec{"mesh file option", "--msh-file", ""}},
      {option_type::MAX_THREADS, option_spec {"maximum number of threads option", "--max-threads", ""}}
    };
    static command_line_options& init_or_get(int argc, char* argv []) {
      static command_line_options instance(argc, argv);
      return instance;
    }
    // Private constructor
    command_line_options(int argc, char *argv[]) {
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "Reading command line";
      for (int idx = 0; idx < argc; ++idx) {
        BOOST_LOG_SEV(rti::mRLogger, blt::trace) << "argv[" << idx << "] == " << argv[idx];
      }
      for (int idx = 0; idx < argc; ++idx) {
        for (auto & oo : options) {
          if (oo.second.optStr.compare(argv[idx]) == 0 && idx < (argc-1)) {
            oo.second.value = argv[idx+1];
            BOOST_LOG_SEV(rti::mRLogger, blt::debug)
              << "command line option " << oo.second.optStr << " " << oo.second.value << " found";
            // Break the inner loop
            break;
          }
        }
      }
    }
  };
} // namespace rti
