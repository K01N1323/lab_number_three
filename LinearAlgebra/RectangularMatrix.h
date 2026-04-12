#ifndef RECTANGULAR_MATRIX_H
#define RECTANGULAR_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"

template <class T> class RectangularMatrix : public BaseMatrix<T> {
private:
  static Sequence<T> *CreateRectangularData(const T *data, int rows, int cols) {
    Sequence<T> *result = new MutableArraySequence<T>();

    if (data == nullptr) {
      for (int count = 0; count < rows * cols; count++) {
        result->append(T(0));
      }
    } else {
      for (int count = 0; count < rows * cols; count++) {
        result->append(data[count]);
      }
    }
    return result;
  }

public:
  RectangularMatrix(int rows, int cols)
      : BaseMatrix<T>(CreateRectangularData(nullptr, rows, cols), rows, cols) {}
  RectangularMatrix(const T *data, int rows, int cols)
      : BaseMatrix<T>(CreateRectangularData(data, rows, cols), rows, cols) {}

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    return new RectangularMatrix<T>(rows, cols);
  }
};

#endif // RECTANGULAR_MATRIX_H