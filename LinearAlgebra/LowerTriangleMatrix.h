#ifndef LOWER_TRIANGULAR_MATRIX_H
#define LOWER_TRIANGULAR_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"
#include <stdexcept>

// Нижнетреугольная матрица: хранит только элементы на главной диагонали и ниже
template <class T> class LowerTriangularMatrix : public BaseMatrix<T> {
public:
  LowerTriangularMatrix(int size) : BaseMatrix<T>(size, size) {
    this->data = new MutableArraySequence<T>();
    for (int index = 0; index < size * (size + 1) / 2; index++) {
      this->data->Append(T(0));
    }
  }

  const T &GetIJ(int row, int col) const override {

    if (row < col) {
      return this->zero_value;
    }

    int index = (row * (row + 1) / 2) + col;

    return this->data->Get(index);
  }

  void SetIJ(int row, int col, const T &item) override {
    if (row < col) {
      if (item == T(0))
        return;
      throw std::invalid_argument(
          "Индекс вне рабочей зоны нижнетреугольной матрицы");
    }

    int index = (row * (row + 1) / 2) + col;
    this->data->Set(index, item);
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    if (rows != cols)
      throw std::invalid_argument("Треугольная матрица должна быть квадратной");
    return new LowerTriangularMatrix<T>(rows);
  }
};

#endif // LOWER_TRIANGULAR_MATRIX_H