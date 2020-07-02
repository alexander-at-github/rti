#pragma once

#include <assert.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>

namespace rti { namespace util { namespace clo {

  // A command line option base class
  class i_option {
  public:
    i_option(std::string pIdStr, std::vector<std::string> pStrings, std::string pHelpStr) :
      mIdStr(pIdStr),
      mOptStrs(pStrings),
      mHelpStr(pHelpStr) {}
    virtual ~i_option() {};
    //protected:
    std::string mIdStr;
    std::vector<std::string> mOptStrs;
    bool mNecessary;
    std::string mHelpStr;
  };

  // An option of the form "-o" or "--foo"
  class bool_option : public i_option {
  public:
    // Inheriting the constructor from i_option
    using i_option::i_option;
    //private:
    bool value = false;
  };

  // An option of the form "-g <string>" or "--g-man <string>"
  class string_option : public i_option {
  public:
    using i_option::i_option;
    string_option(std::string pIdStr,
                  std::vector<std::string> pStrings,
                  std::string pHelpStr,
                  bool pNecessary) :
      i_option(pIdStr, pStrings, pHelpStr),
      mNecessary(pNecessary) {}
    //private:
    bool mNecessary;
    std::string value;
  };

  // Option Manager
  class manager {
  public:
    manager() {}

    manager* addCmlParam(bool_option pBoolOpt) {
      mBoolOpts.insert({pBoolOpt.mIdStr, pBoolOpt});
      return this;
    }

    manager* addCmlParam(string_option pStrOpt) {
      mStrOpts.insert({pStrOpt.mIdStr, pStrOpt});
      return this;
    }

    bool parse_args(int argc, char** argv) {
      mArgc = argc;
      mArgv = argv;
      // starting from 1 because argv[0] is equal to the name of the executable
      for (int idx = 1; idx < argc; /*empty*/) {
        bool argUsed = false;
        for (auto& som : mStrOpts) { // Checking string options
          auto& so = som.second; // second equals the value in the map entry
          for (auto& str : so.mOptStrs) {
            if (str == argv[idx]) {
              if (idx >= argc - 1) {
                // Missing argument for this command line options
                return false;
              } else {
                so.value = std::string(argv[idx+1]);
                argUsed = true;
              }
            }
          }
        }
        if (argUsed) {
          idx += 2; // option string and its argument are parsed
          // It is important to increase the idx variable not earlier than that.
          // Otherwise it is really difficult to tell how often it will be increased.
          continue;
        }
        for (auto& bom : mBoolOpts) { // Checking boolean options
          auto& bo = bom.second; // second equals the value in the map entry
          for (auto& str : bo.mOptStrs) {
            if (str == argv[idx]) {
              bo.value = true;
              argUsed = true;
            }
          }
        }
        if (argUsed) {
          idx += 1;
          // It is important to increase the idx variable not earlier than that.
          // Otherwise it is really difficult to tell how often it will be increased.
          continue;
        }
        // else
        // We do not know the meaning of the given argument
        return false;
      }
      // Check if all necessary arguments have been passed
      for (auto& som : mStrOpts) {
        auto so = som.second; // second equals the value in the map entry
        if (so.mNecessary && so.value == "") // Necessary argument was not given.
          return false;
      }
      return true;
    }

    bool get_bool_option_value(std::string pOptIdStr) {
      for (auto& bome : mBoolOpts)
        if (bome.first == pOptIdStr)
          return bome.second.value;
      std::cerr << "Error in retrieving command line option value" << std::endl;
      assert(false && "Error in clo");
      return false;
    }

    std::string get_string_option_value(std::string pOptIdStr) {
      for (auto& some : mStrOpts)
        if (some.first == pOptIdStr)
          return some.second.value;
      std::cerr << "Error in retrieving command line option value" << std::endl;
      assert(false && "Error in clo");
      return "";
    }

    std::string get_usage_msg() {
      std::stringstream msg;
      msg << "Usage: " << mArgv[0] << " [options]";
      // Only string options can be necessary
      for (auto& som : mStrOpts) {
        auto& so = som.second;
        if (so.mNecessary)
          msg << " " << so.mOptStrs[0] << " <value>"; // Here we use simply the first option string
      }
      msg << std::endl << std::endl;
      msg << "Options:" << std::endl;
      std::string _or = "  or  ";
      std::string spaces = "    ";
      for (auto& bom : mBoolOpts) {
        auto& bo = bom.second;
        std::stringstream tmp;
        tmp << spaces;
        for (auto& os : bo.mOptStrs) {
          // write the options string and space with a word denoting disjunction (or).
          tmp << os << _or;
        }
        std::string tmpS = tmp.str();
        // Remove the last word
        tmpS.erase(tmpS.end() - _or.size(), tmpS.end());
        msg << tmpS << std::endl;
        // Add the help message as specified by the user of this class
        msg << spaces << "   " << bo.mHelpStr << std::endl; // add two spaces of indent of helper string
      }
      for (auto& som : mStrOpts) {
        auto& so = som.second;
        std::stringstream tmp;
        tmp << spaces;
        for (auto& os : so.mOptStrs) {
          tmp << os << " <value>" << _or;
        }
        std::string tmpS = tmp.str();
        // Remove the last word
        tmpS.erase(tmpS.end() - _or.size(), tmpS.end());
        msg << tmpS << std::endl;
        // Add the help message as specified by the user of this class
        msg << spaces << "   " << so.mHelpStr << std::endl;
      }
      return msg.str();
    }

  private:
    std::unordered_map<std::string, bool_option> mBoolOpts;
    std::unordered_map<std::string, string_option> mStrOpts;
    int mArgc = 0;
    char **mArgv = nullptr;
  };
}}}
