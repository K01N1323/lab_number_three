#include "../Sequences/MutableArraySequence.h"
#include "GaussMethod.h"
#include <cmath>

template <class T> T clone_mapper(const T &val) { return val; }
template <class T> T zero_mapper(const T &val) { return T(0); }

template <class T> GaussMethod<T>::GaussMethod(const Matrix<T> *input) {
  this->n = input->GetRows();
  this->det = T(1.0);
  this->swaps = 0;

  this->MatrixPtr = input->Map(clone_mapper<T>);
  this->ConMatrixPtr = input->Map(zero_mapper<T>);

  for (int row = 0; row < n; row++) {
    this->ConMatrixPtr->SetIJ(row, row, T(1.0));
  }
}

template <class T> GaussMethod<T>::~GaussMethod() {
  delete MatrixPtr;
  delete ConMatrixPtr;
}

template <class T> T GaussMethod<T>::GetDet() const {
  if (swaps % 2 == 0)
    return det;
  else
    return T(-1.0) * det;
}

template <class T> int GaussMethod<T>::neededrow(int skip) const {
  T mx = std::abs(MatrixPtr->GetIJ(skip, skip));
  int numrow = skip;

  for (int row = skip + 1; row < n; row++) {
    T loc = std::abs(MatrixPtr->GetIJ(row, skip));
    if (loc > mx) {
      mx = loc;
      numrow = row;
    }
  }

  return numrow;
}

template <class T> void GaussMethod<T>::swaprows(int we_swap1) {
  int we_swap = neededrow(we_swap1);

  if (we_swap1 == we_swap)
    return;

  swaps++;

  for (int col = 0; col < n; col++) {
    T loc = MatrixPtr->GetIJ(we_swap1, col);
    T loccon = ConMatrixPtr->GetIJ(we_swap1, col);

    MatrixPtr->SetIJ(we_swap1, col, MatrixPtr->GetIJ(we_swap, col));
    MatrixPtr->SetIJ(we_swap, col, loc);

    ConMatrixPtr->SetIJ(we_swap1, col, ConMatrixPtr->GetIJ(we_swap, col));
    ConMatrixPtr->SetIJ(we_swap, col, loccon);
  }
}

template <class T> void GaussMethod<T>::divisionrow(int num) {
  T el = MatrixPtr->GetIJ(num, num);

  if (el == T(0)) {
    det = T(0);
    return;
  }

  det *= el;

  for (int col = 0; col < n; col++) {
    MatrixPtr->SetIJ(num, col, MatrixPtr->GetIJ(num, col) / el);
    ConMatrixPtr->SetIJ(num, col, ConMatrixPtr->GetIJ(num, col) / el);
  }
}

template <class T> void GaussMethod<T>::subtraction(int current) {
  for (int row = current + 1; row < n; row++) {
    T mnozh = MatrixPtr->GetIJ(row, current);
    for (int col = 0; col < n; col++) {
      MatrixPtr->SetIJ(row, col,
                        MatrixPtr->GetIJ(row, col) -
                            mnozh * MatrixPtr->GetIJ(current, col));
      ConMatrixPtr->SetIJ(row, col,
                           ConMatrixPtr->GetIJ(row, col) -
                               mnozh * ConMatrixPtr->GetIJ(current, col));
    }
  }
}

template <class T> void GaussMethod<T>::triangle() {
  for (int row = 0; row < n; row++) {
    swaprows(row);
    divisionrow(row);
    subtraction(row);
  }
}

template <class T> void GaussMethod<T>::obrat() {
  for (int pivotRow = n - 1; pivotRow >= 0; pivotRow--) {
    for (int targetRow = pivotRow - 1; targetRow >= 0; targetRow--) {
      T koof = MatrixPtr->GetIJ(targetRow, pivotRow);
      for (int col = 0; col < n; col++) {
        ConMatrixPtr->SetIJ(targetRow, col,
                             ConMatrixPtr->GetIJ(targetRow, col) -
                                 koof * ConMatrixPtr->GetIJ(pivotRow, col));
        MatrixPtr->SetIJ(targetRow, col,
                          MatrixPtr->GetIJ(targetRow, col) -
                              koof * MatrixPtr->GetIJ(pivotRow, col));
      }
    }
  }
}

template <class T> void GaussMethod<T>::reverse() {
  triangle();
  obrat();
}

template <class T> Sequence<T> *GaussMethod<T>::Solve(const Sequence<T> *b) {
  for (int row = 0; row < n; row++) {
    ConMatrixPtr->SetIJ(row, 0, b->Get(row));
  }

  triangle();
  obrat();

  Sequence<T> *result = new MutableArraySequence<T>();
  for (int row = 0; row < n; row++) {
    result->Append(ConMatrixPtr->GetIJ(row, 0));
  }

  return result;
}

template <class T> void GaussMethod<T>::SolveForTests(const Sequence<T> *b) {
  for (int row = 0; row < n; row++) {
    ConMatrixPtr->SetIJ(row, 0, b->Get(row));
  }
  triangle();
  obrat();
}

template <class T> void GaussMethod<T>::TakeReverse() { reverse(); }

template <class T> Matrix<T> *GaussMethod<T>::GetInverseMatrix() const {
  return ConMatrixPtr->Map(clone_mapper<T>);
}

template <class T> T GaussMethod<T>::GetElement(int row, int col) const {
  return MatrixPtr->GetIJ(row, col);
}
