#pragma once

#include <vector>

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

    // A hack to not need to implement the c++ STL iterator.
    // This iterator cannot modify the content.
    std::vector<T> get_iterable() {
      return {frst, scnd, thrd};
    }

    void print(std::ostream& pOs) const {
      pOs << "(" << frst << " " << scnd << " " << thrd << ")";
    }

    bool operator==(const triple<T>& pO) const {
      return frst == pO.frst && scnd == pO.scnd && thrd == pO.thrd;
    }

    bool operator!=(const triple<T>& pO) const {
      return !(*this == pO);
    }
  };

  // This function modifies the arguments when called
  template<typename T>
  void scale(triple<T>& pT, T pF) {
    pT.frst *= pF;
    pT.scnd *= pF;
    pT.thrd *= pF;
  }

  template<typename T>
  T dot_product(const triple<T>& pF, const triple<T>& pS) {
    return pF.frst * pS.frst + pF.scnd * pS.scnd + pF.thrd * pS.thrd;
  }

  template<typename T>
  triple<T> cross_product(const triple<T>& pF, const triple<T>& pS) {
    triple<T> rr;
    rr.frst = pF.scnd * pS.thrd - pF.thrd * pS.scnd;
    rr.scnd = pF.thrd * pS.frst - pF.frst * pS.thrd;
    rr.thrd = pF.frst * pS.scnd - pF.scnd * pS.frst;
    return rr;
  }

  template<typename T>
  triple<T> sum(const triple<T>& pF, const triple<T>& pS) {
    return {pF.frst + pS.frst, pF.scnd + pS.scnd, pF.thrd + pS.thrd};
  }

  template<typename T>
  triple<T> sum(const triple<T>& pF, const triple<T>& pS, const triple<T>& pT) {
    return {pF.frst + pS.frst + pT.frst, pF.scnd + pS.scnd + pT.scnd, pF.thrd + pS.thrd + pT.thrd};
  }

  template<typename T>
  triple<T> inv(const triple<T>& pT) {
    return {-pT.frst, -pT.scnd, -pT.thrd};
  }

  template<typename T>
  triple<T> diff(const triple<T>& pF, const triple<T>& pS) {
    return sum(pF, inv(pS));
  }

  // This function modifies its argument when called
  template<typename T>
  void normalize(triple<T>& pV) {
    T thrdNorm = std::sqrt(
      pV.frst * pV.frst + pV.scnd * pV.scnd + pV.thrd * pV.thrd);
    pV.frst /= thrdNorm;
    pV.scnd /= thrdNorm;
    pV.thrd /= thrdNorm;
  }
} // namespace rti
