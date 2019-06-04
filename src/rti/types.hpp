#pragma once

namespace rti {
  template<typename T>
  class pair {
  public:
    T frst, scnd;
  };

  template<typename T>
  class triple {
  public:
    T frst, scnd, thrd;
    // A hack to not need to implement the iterator, which is really
    // complicated in C++.
    // This iterator cannot modify the content.
    std::vector<T> get_iterable() {
      return {frst, scnd, thrd};
    }
  };
} // namespace rti
