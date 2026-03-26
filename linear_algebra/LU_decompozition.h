#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include "Matrix.h"
#include "Sequence.h"
#include <stdexcept>

template <class T> class LU_Decomposition {
public:
  LU_Decomposition(const Matrix<T> *input);
  ~LU_Decomposition();

  Matrix<T> *GetL() const;
  Matrix<T> *GetU() const;
  Sequence<int> *GetP() const;

  T GetDet() const;

  Sequence<T> *Solve(const Sequence<T> *b) const;

private:
  int n;

  Matrix<T> *L;
  Matrix<T> *U;
  Sequence<int> *P;

  int swaps;

  int neededrow(int skip) const;
  void swaprows(int step);
  void subtraction(int current);
  void decompose();
};

#include "LU_Decomposition.tpp"

#endif // LU_DECOMPOSITION_H