#ifndef GAUSS_METHOD_H
#define GAUSS_METHOD_H

#include "Sequence.h"
#include "matrix.h"

template <class T> class gauss_method {
public:
  gauss_method(const matrix<T> *input);
  ~gauss_method();

  matrix<T> *GetInverseMatrix() const;
  T GetDet() const;
  void TakeReverse();

  Sequence<T> *Solve(const Sequence<T> *b);
  void SolveForTests(const Sequence<T> *b);

private:
  int n;

  matrix<T> *matrix_ptr;
  matrix<T> *conmatrix_ptr;

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

#include "gauss_method.tpp"

#endif // GAUSS_METHOD_H