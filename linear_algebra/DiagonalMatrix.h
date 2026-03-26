#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "BaseMatrix.h"
#include "MutableArraySequence.h"

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
      return this->zero_value;
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

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    return new DiagonalMatrix<T>(rows);
  }
};

#endif // DIAGONAL_MATRIX_H