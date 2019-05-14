#pragma once

namespace rti {
  class test_result {
  public:
    std::string to_string() {
      std::stringstream strstream;
      strstream << "(:class test_result :not-implemented)";
      return strstream.str();
    }
  };
} // namespace rti
