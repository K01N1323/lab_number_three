#ifndef RECTANGULAR_MATRIX_H
#define RECTANGULAR_MATRIX_H

#include "BaseMatrix.h"
#include "MutableArraySequence.h"

template <class T> class RectangularMatrix : public BaseMatrix<T> {
public:
  RectangularMatrix(int rows, int cols) : BaseMatrix<T>(rows, cols) {
    this->data = new MutableArraySequence<T>();
    for (int index = 0; index < rows * cols; index++) {
      this->data->Append(T(0));
    }
  }

protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    return new RectangularMatrix<T>(rows, cols);
  }
};

#endif // RECTANGULAR_MATRIX_H