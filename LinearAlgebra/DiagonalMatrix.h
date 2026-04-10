#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"

template <class T> class DiagonalMatrix : public BaseMatrix<T> {
public:
  DiagonalMatrix(int size) : BaseMatrix<T>(size, size) {
    this->data = new MutableArraySequence<T>();
    for (int index = 0; index < size; index++) {
      this->data->Append(T(0));
    }
  }

  const T &GetIJ(int row, int col) const override {
    if (row != col) {
      return this->ZeroValue;
    }
    return this->data->Get(row);
  }

  void SetIJ(int row, int col, const T &item) override {
    if (row != col && item != T(0)) {
      throw std::invalid_argument(
          "Нельзя задать ненулевой элемент вне диагонали");
    } else if (row == col) {
      this->data->Set(row, item);
    }
  }

  T GetDet() const override {
    T det = T(1);
    for (int row = 0; row < this->rows; row++) {
        det *= this->GetIJ(row, row);
    }
    return det;
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    return new DiagonalMatrix<T>(rows);
  }
};

#endif // DIAGONAL_MATRIX_H