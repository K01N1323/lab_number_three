#ifndef LU_DECOMPOZITION_H
#define LU_DECOMPOZITION_H

#include "../Sequences/Sequence.h"
#include "Matrix.h"
#include <stdexcept>

template <class T> class LUDecompozition {
public:

  LUDecompozition(const Matrix<T> *input);
  ~LUDecompozition();

  // Доступ к результатам разложения
  Matrix<T> *GetL() const;
  Matrix<T> *GetU() const;
  Sequence<int> *GetP() const;

  T GetDet() const;

  // Решение СЛАУ Ax = b методом LU-разложения
  Sequence<T> *Solve(const Sequence<T> *b) const;

private:

  int n;
  Matrix<T> *L;
  Matrix<T> *U;
  Sequence<int> *P;
  int swaps;

  // Вспомогательные методы разложения
  int neededrow(int skip) const;
  void swaprows(int step);
  void subtraction(int current);
  void decompose();
};

#include "LUDecompozition.tpp"

#endif // LU_DECOMPOZITION_H
