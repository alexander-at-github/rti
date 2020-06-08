#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <iostream>

namespace rti { namespace util {

  //constexpr double pi() { return std::atan(1)*4; }
  constexpr auto pi = std::acos(-1);

  template<typename Ty>
  using pair = std::array<Ty, 2>;

  template<typename Ty>
  using triple = std::array<Ty, 3>;

  template<typename Ty>
  using quadruple = std::array<Ty, 4>;

  template<typename T0, typename T1, typename T2>
  class tripleN {
  private:
    T0 i0;
    T1 i1;
    T2 i2;
  public:
    tripleN() {}
    tripleN(T0 i0, T1 i1, T2 i2) : i0(i0), i1(i1), i2(i2) {}
    T0& get0() { return i0; }
    T1& get1() { return i1; }
    T2& get2() { return i2; }
  };

  // This function modifies the arguments when called
  template<typename Ty>
  triple<Ty> scale(Ty pF, triple<Ty>& pT)
  {
    pT[0] *= pF;
    pT[1] *= pF;
    pT[2] *= pF;
    return pT;
  }

  template<typename Ty>
  Ty dot_product(const triple<Ty>& pF, const triple<Ty>& pS)
  {
    return pF[0] * pS[0] + pF[1] * pS[1] + pF[2] * pS[2];
  }

  template<typename Ty>
  triple<Ty> cross_product(const triple<Ty>& pF, const triple<Ty>& pS)
  {
    triple<Ty> rr;
    rr[0] = pF[1] * pS[2] - pF[2] * pS[1];
    rr[1] = pF[2] * pS[0] - pF[0] * pS[2];
    rr[2] = pF[0] * pS[1] - pF[1] * pS[0];
    return rr;
  }

  template<typename Ty>
  triple<Ty> sum(const triple<Ty>& pF, const triple<Ty>& pS)
  {
    return {pF[0] + pS[0], pF[1] + pS[1], pF[2] + pS[2]};
  }

  template<typename Ty>
  triple<Ty> sum(const triple<Ty>& pF, const triple<Ty>& pS, const triple<Ty>& pT)
  {
    return {pF[0] + pS[0] + pT[0], pF[1] + pS[1] + pT[1], pF[2] + pS[2] + pT[2]};
  }

  template<typename Ty>
  triple<Ty> inv(const triple<Ty>& pT)
  {
    return {-pT[0], -pT[1], -pT[2]};
  }

  template<typename Ty>
  triple<Ty> diff(const triple<Ty>& pF, const triple<Ty>& pS)
  {
    return sum(pF, inv(pS));
  }

  template<typename Ty>
  bool is_normalized(triple<Ty> const& pV)
  {
    auto epsilon = 1e-6f;
    Ty length = std::sqrt(pV[0] * pV[0] + pV[1] * pV[1] + pV[2] * pV[2]);
    return 1-epsilon <= length && length <= 1+epsilon;
  }

    // This function modifies its argument when called
    template<typename Ty>
    void normalize(triple<Ty>& pV)
    {
      Ty thrdNorm = std::sqrt(pV[0] * pV[0] + pV[1] * pV[1] + pV[2] * pV[2]);
      pV[0] /= thrdNorm;
      pV[1] /= thrdNorm;
      pV[2] /= thrdNorm;
      if ( ! is_normalized(pV) )
        std::cerr << "Warning: Assertion error is about to happen. thrdNorm == " << thrdNorm << std::endl;
      assert( is_normalized(pV) && "Postcondition" );
    }

  // Compute normal of a triangle
  template<typename Ty>
  rti::util::triple<Ty> compute_normal(rti::util::triple<rti::util::triple<Ty> >& pTri)
  {
    auto uu = rti::util::diff(pTri[1], pTri[0]);
    auto vv = rti::util::diff(pTri[2], pTri[0]);
    return rti::util::cross_product(uu, vv);
  }

  template<typename Ty>
  bool each_normalized(std::vector<rti::util::triple<Ty> > pP)
  {
    for (auto const& nn : pP)
      if ( ! rti::util::is_normalized(nn))
        return false;
    return true;
  }

  template<typename Ty>
  bool contains(rti::util::triple<Ty>& pT, Ty pE)
  {
    return pT[0] == pE || pT[1] == pE || pT[2] == pE;
  }

  // Input: a triple of triples where each inner triple holds the x, y and z coordinates
  // of a vertex of a triangle.
  template<typename Ty>
  static rti::util::triple<Ty> centroid(rti::util::triple<rti::util::triple<Ty> > pTriangle)
  {
    rti::util::triple<Ty> result;
    result[0] = (pTriangle[0][0] + pTriangle[1][0] + pTriangle[2][0]) / 3;
    result[1] = (pTriangle[0][1] + pTriangle[1][1] + pTriangle[2][1]) / 3;
    result[2] = (pTriangle[0][2] + pTriangle[1][2] + pTriangle[2][2]) / 3;
    return result;
  }

  template<typename Ty>
  static Ty distance(rti::util::pair<rti::util::triple<Ty> > pPnts)
  {
    auto p1 = pPnts[0][0] - pPnts[1][0];
    auto p2 = pPnts[0][1] - pPnts[1][1];
    auto p3 = pPnts[0][2] - pPnts[1][2];
    return std::sqrt(p1*p1 + p2*p2 + p3*p3);
  }

  template<typename Ty>
  static Ty length_of_vec(rti::util::triple<Ty> vec)
  {
    return distance<Ty>({vec, {0,0,0}});
  }

  template<typename Ty>
  static Ty area_of_triangle(rti::util::triple<rti::util::triple<Ty> > tri)
  {
    auto ab = diff(tri[1], tri[0]);
    auto ac = diff(tri[2], tri[0]);
    auto crossp = cross_product(ab, ac);
    auto length = length_of_vec(crossp);
    return ((Ty) 0.5) * length;
  }

  // A c-style-array foldl as generic auxiliary implementation.
  // This implementation should be reasonably efficient. It tries to minimize copies
  // of data.
  template<typename T1, typename T2>
  static T1& foldl_aux(std::function<T1 (T1&, T2 const&)>& pF,
                       T1& pT1,
                       T2 const* pT2,
                       size_t& pT2Length)
  {
    if (pT2Length <= 0) return pT1;
    pT1 = pF(pT1, pT2[0]); // apply the fold-function // modify the content of pT1
    return foldl_aux<T1, T2>(pF, pT1, &pT2[1], --pT2Length /* modify pT2Length */);
  }

  // A foldl for std::vector
  template<typename T1, typename T2>
  static T1 foldl(std::function<T1 (T1&, T2 const&)> pF,
                  T1 pT1,
                  std::vector<T2> pT2)
  {
    auto pT2Size = pT2.size(); // copy the size once into a separate memory location.
    return foldl_aux<T1,T2>(pF, pT1, pT2.data(), pT2Size);
  }

  template<typename TT>
  static void swap(TT& e1, TT& e2)
  {
    TT tmp = e2;
    e2 = e1;
    e1 = tmp;
  }

  // Returns some orthonormal basis containing a the input vector pVector
  // (possibly scaled) as the first element of the return value.
  // This function is deterministic, i.e., for one input it will return always
  // the same result.
  template<typename Ty> static
  rti::util::triple<rti::util::triple<Ty> >
  get_orthonormal_basis(const rti::util::triple<Ty> pVector) {
    rti::util::triple<rti::util::triple<Ty> > rr;
    rr[0] = pVector;

    // Calculate a vector (rr[1]) which is perpendicular to rr[0]
    // https://math.stackexchange.com/questions/137362/how-to-find-perpendicular-vector-to-another-vector#answer-211195
    rti::util::triple<Ty> candidate0 {rr[0][2], rr[0][2], -(rr[0][0] + rr[0][1])};
    rti::util::triple<Ty> candidate1 {rr[0][1], -(rr[0][0] + rr[0][2]), rr[0][1]};
    rti::util::triple<Ty> candidate2 {-(rr[0][1] + rr[0][2]), rr[0][0], rr[0][0]};
    // We choose the candidate which maximizes the sum of its components, because we
    // want to avoid numeric errors and that the result is (0, 0, 0).
    std::array<rti::util::triple<Ty>, 3> cc = {candidate0, candidate1, candidate2};
    auto sumFun = [](rti::util::triple<Ty> oo){return oo[0] + oo[1] + oo[2];};
    int maxIdx = 0;
    for (size_t idx = 1; idx < cc.size(); ++idx) {
      if (sumFun(cc[idx]) > sumFun(cc[maxIdx])) {
        maxIdx = idx;
      }
    }
    assert (maxIdx < 3 && "Error in computation of perpenticular vector");
    rr[1] = cc[maxIdx];

    rr[2] = rti::util::cross_product(rr[0], rr[1]);
    rti::util::normalize(rr[0]);
    rti::util::normalize(rr[1]);
    rti::util::normalize(rr[2]);
    // Sanity check
    Ty epsilon = 1e-6;
    assert(std::abs(rti::util::dot_product(rr[0], rr[1])) < epsilon &&
           "Error in orthonormal basis computation");
    assert(std::abs(rti::util::dot_product(rr[1], rr[2])) < epsilon &&
           "Error in orthonormal basis computation");
    assert(std::abs(rti::util::dot_product(rr[2], rr[0])) < epsilon &&
           "Error in orthonormal basis computation");
    return rr;
    }
}}
