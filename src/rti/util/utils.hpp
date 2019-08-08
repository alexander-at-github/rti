#pragma once

#include <array>
#include <cmath>

namespace rti { namespace util {

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

  // template<typename Ty>
  // class pair {
  // public:
  //   Ty frst, scnd;
  // };

  // template<typename Ty>
  // class quadruple {
  // public:
  //   Ty frst, scnd, thrd, frth;
  // }

  // template<typename Ty>
  // class triple {
  // public:
  //   Ty frst, scnd, thrd;

  //   // A hack to not need to implement the c++ STL iterator.
  //   // This iterator cannot modify the content.
  //   std::vector<Ty> get_iterable() const {
  //     return {frst, scnd, thrd};
  //   }

  //   void print(std::ostream& pOs) const {
  //     pOs << "(" << frst << " " << scnd << " " << thrd << ")";
  //   }

  //   bool operator==(const triple<Ty>& pO) const {
  //     return frst == pO.frst && scnd == pO.scnd && thrd == pO.thrd;
  //   }

  //   bool operator!=(const triple<Ty>& pO) const {
  //     return !(*this == pO);
  //   }
  // };

  // This function modifies the arguments when called
  template<typename Ty>
  triple<Ty> scale(Ty pF, triple<Ty>& pT) {
    pT[0] *= pF;
    pT[1] *= pF;
    pT[2] *= pF;
    return pT;
  }

  template<typename Ty>
  Ty dot_product(const triple<Ty>& pF, const triple<Ty>& pS) {
    return pF[0] * pS[0] + pF[1] * pS[1] + pF[2] * pS[2];
  }

  template<typename Ty>
  triple<Ty> cross_product(const triple<Ty>& pF, const triple<Ty>& pS) {
    triple<Ty> rr;
    rr[0] = pF[1] * pS[2] - pF[2] * pS[1];
    rr[1] = pF[2] * pS[0] - pF[0] * pS[2];
    rr[2] = pF[0] * pS[1] - pF[1] * pS[0];
    return rr;
  }

  template<typename Ty>
  triple<Ty> sum(const triple<Ty>& pF, const triple<Ty>& pS) {
    return {pF[0] + pS[0], pF[1] + pS[1], pF[2] + pS[2]};
  }

  template<typename Ty>
  triple<Ty> sum(const triple<Ty>& pF, const triple<Ty>& pS, const triple<Ty>& pT) {
    return {pF[0] + pS[0] + pT[0], pF[1] + pS[1] + pT[1], pF[2] + pS[2] + pT[2]};
  }

  template<typename Ty>
  triple<Ty> inv(const triple<Ty>& pT) {
    return {-pT[0], -pT[1], -pT[2]};
  }

  template<typename Ty>
  triple<Ty> diff(const triple<Ty>& pF, const triple<Ty>& pS) {
    return sum(pF, inv(pS));
  }

  // This function modifies its argument when called
  template<typename Ty>
  void normalize(triple<Ty>& pV) {
    Ty thrdNorm = std::sqrt(pV[0] * pV[0] + pV[1] * pV[1] + pV[2] * pV[2]);
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
  rti::util::triple<Ty> compute_normal(rti::util::triple<rti::util::triple<Ty> >& pTri) {
    auto uu = rti::util::diff(pTri[1], pTri[0]);
    auto vv = rti::util::diff(pTri[2], pTri[0]);
    return rti::util::cross_product(uu, vv);
  }

  template<typename Ty>
  bool each_normalized(std::vector<rti::util::triple<Ty> > pP) {
    for (auto const& nn : pP)
      if ( ! rti::util::is_normalized(nn))
        return false;
    return true;
  }

}} // namespace
