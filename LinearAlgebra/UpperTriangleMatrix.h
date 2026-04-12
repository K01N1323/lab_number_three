#ifndef UPPER_TRIANGULAR_MATRIX_H
#define UPPER_TRIANGULAR_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"
#include <stdexcept>

template <class T> class UpperTriangularMatrix : public BaseMatrix<T> {
private:
  static Sequence<T> *CreateUpperTriangleData(const T *data, int size) {
    Sequence<T> *result = new MutableArraySequence<T>();

    if (data == nullptr) {
      for (int index = 0; index < size * (size + 1) / 2; index++) {
        result->append(T(0));
      }
    } else {
      for (int count = 0; count < size * (size + 1) / 2; count++) {
        result->append(data[count]);
      }
    }

    return result;
  }

public:
  UpperTriangularMatrix(int size)
      : BaseMatrix<T>(CreateUpperTriangleData(nullptr, size), size, size) {}

  UpperTriangularMatrix(const T *data, int size)
      : BaseMatrix<T>(CreateUpperTriangleData(data, size), size, size) {}

  const T &GetIJ(int row, int col) const override {

    if (row > col) {
      return this->ZeroValue;
    }

    int n = this->rows;
    int index = (row * (2 * n - row + 1)) / 2 + (col - row);

    return this->data->get(index);
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

    this->data->set(index, item);
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    if (rows != cols)
      throw std::invalid_argument("Треугольная матрица должна быть квадратной");
    return new UpperTriangularMatrix<T>(rows);
  }
};

#endif