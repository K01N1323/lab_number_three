#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"

template <class T> class SquareMatrix : public BaseMatrix<T> {
private:
  static Sequence<T> *CreateSquareData(const T *data, int size) {
    Sequence<T> *result = new MutableArraySequence<T>();

    if (data == nullptr) {
      for (int count = 0; count < size * size; count++) {
        result->append(T(0));
      }
    } else {
      for (int count = 0; count < size * size; count++) {
        result->append(data[count]);
      }
    }

    return result;
  }

public:
  SquareMatrix(int size)
      : BaseMatrix<T>(CreateSquareData(nullptr, size), size, size) {}

  SquareMatrix(const T *data, int size)
      : BaseMatrix<T>(CreateSquareData(data, size), size, size) {}

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    if (rows != cols)
      throw std::invalid_argument("Ожидается квадратная матрица");
    return new SquareMatrix<T>(rows);
  }
};

#endif // SQUARE_MATRIX_H