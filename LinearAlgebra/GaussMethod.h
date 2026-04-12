#ifndef GAUSS_METHOD_H
#define GAUSS_METHOD_H

#include "../Sequences/Sequence.h"
#include "Matrix.h"

template <class T> class GaussMethod {
public:

  GaussMethod(const Matrix<T> *input);
  ~GaussMethod();

  Matrix<T> *GetInverseMatrix() const;
  T GetDet() const;
  void TakeReverse();

  Sequence<T> *Solve(const Sequence<T> *b);
  void SolveForTests(const Sequence<T> *b);

private:

  int n;
  Matrix<T> *MatrixPtr;
  Matrix<T> *ConMatrixPtr;

  T det;
  int swaps;

  int NeededRow(int skip) const;
  void SwapRows(int);
  void DivisionRow(int);
  void subtraction(int current);
  void triangle(void);
  void obrat(void);
  void reverse(void);
  T GetElement(int, int) const;
};

#include "GaussMethod.tpp"

#endif // GAUSS_METHOD_H
