#ifndef LOWER_TRIANGULAR_MATRIX_H
#define LOWER_TRIANGULAR_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"
#include <stdexcept>

template <class T> class LowerTriangularMatrix : public BaseMatrix<T> {
private:
  static Sequence<T> *CreateLowerTriangleData(const T *data, int size) {
    Sequence<T> *result = new MutableArraySequence<T>();

    if (data == nullptr) {
      for (int count = 0; count < size * (size + 1) / 2; count++) {
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
  LowerTriangularMatrix(int size)
      : BaseMatrix<T>(CreateLowerTriangleData(nullptr, size), size, size) {}

  LowerTriangularMatrix(const T *data, int size)
      : BaseMatrix<T>(CreateLowerTriangleData(data, size), size, size) {}

  const T &GetIJ(int row, int col) const override {

    if (row < col) {
      return this->ZeroValue;
    }

    int index = (row * (row + 1) / 2) + col;

    return this->data->get(index);
  }

  void SetIJ(int row, int col, const T &item) override {
    if (row < col) {
      if (item == T(0))
        return;
      throw std::invalid_argument(
          "Индекс вне рабочей зоны нижнетреугольной матрицы");
    }

    int index = (row * (row + 1) / 2) + col;
    this->data->set(index, item);
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    if (rows != cols)
      throw std::invalid_argument("Треугольная матрица должна быть квадратной");
    return new LowerTriangularMatrix<T>(rows);
  }
};

#endif // LOWER_TRIANGULAR_MATRIX_H