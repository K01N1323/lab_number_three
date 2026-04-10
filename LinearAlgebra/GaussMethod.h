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

  int neededrow(int skip) const;
  void swaprows(int);
  void divisionrow(int);
  void subtraction(int);
  void triangle();
  void obrat();
  void reverse();
  T GetElement(int, int) const;
};

#include "GaussMethod.tpp"

#endif // GAUSS_METHOD_H
