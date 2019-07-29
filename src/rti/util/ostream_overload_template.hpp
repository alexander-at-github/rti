#pragma once

namespace rti { namespace util {
  // This template allows you to print any object with the '<<' operator; if
  // the object has a 'T::print(std::ostream&) const' member function
  // (irrespective of the return type). (Note the const specifier, otherwise it
  // does not work.)
  template<typename T>
  auto operator<<(std::ostream& os, const T& t) -> decltype(t.print(os), os) {
    t.print(os);
    return os;
  }
}} // namespace
