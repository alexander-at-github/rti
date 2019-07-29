#pragma once

#include <array>

namespace rti {

  template<typename Ty>
  using pair = std::array<Ty, 2>;
  // class pair : public std::array<Ty, 2> {
  //   template<size_t idx>
  //   Ty aa() {
  //     return std::get<idx>(this);
  //   }
  // };

  template<typename Ty>
  using triple = std::array<Ty, 3>;
  // class triple : public std::array<Ty, 3> {
  //   template<size_t idx>
  //   constexpr Ty const& mat() {
  //     return std::get<idx> (this);
  //   }
  // };

  template<typename Ty>
  using quadruple = std::array<Ty, 4>;
  // class quadruple : public std::array<Ty, 4> {
  // public:
  //   template<size_t idx>
  //   constexpr Ty const& sg() {
  //     return std::get<idx> (*this);
  //   }
  // };

  // template<typename T>
  // class pair {
  // public:
  //   T frst, scnd;
  // };

  // template<typename Ty>
  // class quadruple {
  // public:
  //   Ty frst, scnd, thrd, frth;
  // }

  // template<typename T>
  // class triple {
  // public:
  //   T frst, scnd, thrd;

  //   // A hack to not need to implement the c++ STL iterator.
  //   // This iterator cannot modify the content.
  //   std::vector<T> get_iterable() const {
  //     return {frst, scnd, thrd};
  //   }

  //   void print(std::ostream& pOs) const {
  //     pOs << "(" << frst << " " << scnd << " " << thrd << ")";
  //   }

  //   bool operator==(const triple<T>& pO) const {
  //     return frst == pO.frst && scnd == pO.scnd && thrd == pO.thrd;
  //   }

  //   bool operator!=(const triple<T>& pO) const {
  //     return !(*this == pO);
  //   }
  // };

  // This function modifies the arguments when called
  template<typename T>
  triple<T> scale(T pF, triple<T>& pT) {
    pT[0] *= pF;
    pT[1] *= pF;
    pT[2] *= pF;
    return pT;
  }

  template<typename T>
  T dot_product(const triple<T>& pF, const triple<T>& pS) {
    return pF[0] * pS[0] + pF[1] * pS[1] + pF[2] * pS[2];
  }

  template<typename T>
  triple<T> cross_product(const triple<T>& pF, const triple<T>& pS) {
    triple<T> rr;
    rr[0] = pF[1] * pS[2] - pF[2] * pS[1];
    rr[1] = pF[2] * pS[0] - pF[0] * pS[2];
    rr[2] = pF[0] * pS[1] - pF[1] * pS[0];
    return rr;
  }

  template<typename T>
  triple<T> sum(const triple<T>& pF, const triple<T>& pS) {
    return {pF[0] + pS[0], pF[1] + pS[1], pF[2] + pS[2]};
  }

  template<typename T>
  triple<T> sum(const triple<T>& pF, const triple<T>& pS, const triple<T>& pT) {
    return {pF[0] + pS[0] + pT[0], pF[1] + pS[1] + pT[1], pF[2] + pS[2] + pT[2]};
  }

  template<typename T>
  triple<T> inv(const triple<T>& pT) {
    return {-pT[0], -pT[1], -pT[2]};
  }

  template<typename T>
  triple<T> diff(const triple<T>& pF, const triple<T>& pS) {
    return sum(pF, inv(pS));
  }

  // This function modifies its argument when called
  template<typename T>
  void normalize(triple<T>& pV) {
    T thrdNorm = std::sqrt(pV[0] * pV[0] + pV[1] * pV[1] + pV[2] * pV[2]);
    pV[0] /= thrdNorm;
    pV[1] /= thrdNorm;
    pV[2] /= thrdNorm;
  }

  template<typename Ty>
  bool is_normalized(triple<Ty> const& pV) {
    auto epsilon = 1e-6f;
    Ty length = std::sqrt(pV[0] * pV[0] + pV[1] * pV[1] + pV[2] * pV[2]);
    return 1-epsilon <= length && length <= 1+epsilon;
  }

  // Compute normal of a triangle
  template<typename Ty>
  rti::triple<Ty> compute_normal(rti::triple<rti::triple<Ty> >& pTri) {
    auto uu = rti::diff(pTri[1], pTri[0]);
    auto vv = rti::diff(pTri[2], pTri[0]);
    return rti::cross_product(uu, vv);
  }
} // namespace rti
