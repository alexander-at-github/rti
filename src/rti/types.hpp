#pragma once

namespace rti {

  // As an alternative to that one can define template classes for pairs and triples.
  // The advantage of classes is, that one gets compile time bounds checking.
  // Yes we should do that. TODO
  template<typename T>
  using pair_t = std::array<T, 2>;

  template<typename T>
  using triple_t = std::array<T, 3>;
} // namespace rti
