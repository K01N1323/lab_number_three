#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

#include "BaseMatrix.h"
#include "MutableArraySequence.h"

template <class T> class SquareMatrix : public BaseMatrix<T> {
public:
  SquareMatrix(int size) : BaseMatrix<T>(size, size) {
    this->data = new MutableArraySequence<T>();
    for (int index = 0; index < size * size; index++) {
      this->data->Append(T(0));
    }
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    if (rows != cols)
      throw std::invalid_argument("Ожидается квадратная матрица");
    return new SquareMatrix<T>(rows);
  }
};

#endif // SQUARE_MATRIX_H