#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "../Sequences/MutableArraySequence.h"
#include "BaseMatrix.h"

template <class T> class DiagonalMatrix : public BaseMatrix<T> {
private:
  static Sequence<T> *CreateDiagonalData(const T *data, int size) {
    Sequence<T> *result = new MutableArraySequence<T>();

    if (data == nullptr) {
      for (int count = 0; count < size; count++) {
        result->append(T(0));
      }
    } else {
      for (int count = 0; count < size; count++) {
        result->append(data[count]);
      }
    }

    return result;
  }

public:
  DiagonalMatrix(int size)
      : BaseMatrix<T>(CreateDiagonalData(nullptr, size), size, size) {}

  DiagonalMatrix(const T *data, int size)
      : BaseMatrix<T>(CreateDiagonalData(data, size), size, size) {}

  const T &GetIJ(int row, int col) const override {
    if (row != col) {
      return this->ZeroValue;
    }
    return this->data->get(row);
  }

  void SetIJ(int row, int col, const T &item) override {
    if (row != col && item != T(0)) {
      throw std::invalid_argument(
          "Нельзя задать ненулевой элемент вне диагонали");
    } else if (row == col) {
      this->data->set(row, item);
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