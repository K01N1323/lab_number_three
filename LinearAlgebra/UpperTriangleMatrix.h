#ifndef UPPER_TRIANGULAR_MATRIX_H
#define UPPER_TRIANGULAR_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"
#include <stdexcept>

template <class T> class UpperTriangularMatrix : public BaseMatrix<T> {
public:
  UpperTriangularMatrix(int size) : BaseMatrix<T>(size, size) {
    this->data = new MutableArraySequence<T>();
    for (int index = 0; index < size * (size + 1) / 2; index++) {
      this->data->Append(T(0));
    }
  }

  const T &GetIJ(int row, int col) const override {

    if (row > col) {
      return this->zero_value;
    }

    int n = this->rows;
    int index = (row * (2 * n - row + 1)) / 2 + (col - row);

    return this->data->Get(index);
  }

  void SetIJ(int row, int col, const T &item) override {
    if (row > col) {
      if (item == T(0))
        return;
      throw std::invalid_argument(
          "Индекс вне рабочей зоны верхнетреугольной матрицы");
    }

    int n = this->rows;
    int index = (row * (2 * n - row + 1)) / 2 + (col - row);

    this->data->Set(index, item);
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    if (rows != cols)
      throw std::invalid_argument("Треугольная матрица должна быть квадратной");
    return new UpperTriangularMatrix<T>(rows);
  }
};

#endif