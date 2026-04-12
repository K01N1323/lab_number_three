#ifndef IMMUTABLE_MATRIX_H
#define IMMUTABLE_MATRIX_H

#include "../Sequences/Sequence.h"
#include "Matrix.h"
#include <stdexcept>

template <class T> class ImmutableMatrix : public Matrix<T> {
private:

  Matrix<T> *CoreMatrix;

public:

  ImmutableMatrix(Matrix<T> *matrix) {
    if (matrix == nullptr) {
      throw std::invalid_argument("Передан нулевой указатель на матрицу");
    }
    CoreMatrix = matrix;
  }

  ~ImmutableMatrix() { delete CoreMatrix; }

  int GetRows() const override { return CoreMatrix->GetRows(); }
  int GetCols() const override { return CoreMatrix->GetCols(); }

  const T &GetIJ(int row, int col) const override {
    return CoreMatrix->GetIJ(row, col);
  }

  double GetFrobeniusNorm() const override {
    return CoreMatrix->GetFrobeniusNorm();
  }

  Sequence<T> *SolveSlauGauss(const Sequence<T> *b) const override {
    return CoreMatrix->SolveSlauGauss(b);
  }

  Sequence<T> *SolveSlauLU(const Sequence<T> *b) const override {
    return CoreMatrix->SolveSlauLU(b);
  }

  Matrix<T> *GetL() const override { return CoreMatrix->GetL(); }
  Matrix<T> *GetU() const override { return CoreMatrix->GetU(); }

  void SetIJ(int row, int col, const T &item) override {
    throw std::logic_error("Попытка изменить неизменяемую матрицу");
  }

  void set(int index, const T &item) override {
    throw std::logic_error("Попытка изменить неизменяемую матрицу");
  }

  T GetDet() const override { return CoreMatrix->GetDet(); }
  Matrix<T> *GetInverseMatrix() const override { return CoreMatrix->GetInverseMatrix(); }
  
  Matrix<T> *MakeGilbert() override {
    throw std::logic_error("Попытка изменить неизменяемую матрицу");
  }
  
  Matrix<T> *MakeOnes() override {
    throw std::logic_error("Попытка изменить неизменяемую матрицу");
  }

  Matrix<T> *operator+(const Matrix<T> &rv) const override { return *CoreMatrix + rv; }
  Matrix<T> *operator-(const Matrix<T> &rv) const override { return *CoreMatrix - rv; }
  Matrix<T> *operator*(const Matrix<T> &rv) const override { return *CoreMatrix * rv; }
  Matrix<T> *operator*(const T scalar) const override { return *CoreMatrix * scalar; }

  Matrix<T> *map(T (*mapper)(const T &)) const override { return CoreMatrix->map(mapper); }
  
  T reduce(T (*ReduceFunc)(const T &, const T &), const T &StartValue) const override {
    return CoreMatrix->reduce(ReduceFunc, StartValue);
  }
  
  IEnumerator<T> *GetEnumerator() const override { return CoreMatrix->GetEnumerator(); }

  protected:
  Matrix<T> *CreateEmpty(int rows, int cols) const override {
    throw std::logic_error("CreateEmpty не доступен для неизменяемой матрицы");
  }
};

#endif // IMMUTABLE_MATRIX_H