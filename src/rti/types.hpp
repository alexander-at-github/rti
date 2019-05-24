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
    std::vector<T> getIterable() {
      // TODO: test whether you can modify the content with that.
      return {frst, scnd, thrd};
    }
  };

  // template<typename T>
  // class triple_perfect {
  // public:
  //   T frst, scnd, thrd;
  //   Iterator iter() {
  //     return {frst, scnd, thrd}.iter();
  //   }
  // }

  // template<typename T>
  // class triple_e : public std::tuple<T, T, T> {
  // public:
  //   T gFrst() {
  //     std::get<0>(*this);
  //   }
  //   void sFrst(T v) {
      
  //   }
  // };

  // template<typename T>
  // class triple_ee : public std::array<T, 3> {
  //   z
  // }
} // namespace rti
