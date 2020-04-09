#pragma once

// In the C++11 (gcc 4.9.2) standard there is a problem in hashing enum classes
// which causes compile time errors. This behaviour is considered a defect and
// was fixed in C++14.
//
// This class provides a hash function for an enum class for C++11.

namespace rti { namespace util {
class enum_class_hash_function {
public:
  template<typename Ty>
    std::size_t operator()(Ty t) const {
      return static_cast<std::size_t>(t);
    }
};
}}
