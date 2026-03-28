#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include "../Sequences/Sequence.h"
#include "Matrix.h"

#include <cmath>
#include <stdexcept>

template <class T>
static double sum_of_squares_func(const double &acc, const T &val) {
  double v = static_cast<double>(std::abs(val));
  return acc + v * v;
}

template <class T> class BaseMatrix : public Matrix<T> {
protected:
  Sequence<T> *data;
  int rows;
  int cols;

  T zero_value = T(0);

public:
  BaseMatrix(int rows, int cols) : BaseMatrix(nullptr, rows, cols) {}

  BaseMatrix(Sequence<T> *data, int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    this->data = data;
  }

  virtual ~BaseMatrix() override { delete this->data; }

  virtual const T &GetIJ(int row, int col) const override {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
      throw std::out_of_range("Индекс вне границ матрицы");
    }
    return data->Get(row * cols + col);
  }

  virtual void SetIJ(int row, int col, const T &item) override {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
      throw std::out_of_range("Индекс вне границ матрицы");
    }
    data->Set(row * cols + col, item);
  }

  virtual void Set(int index, const T &item) override {
    if (index >= rows * cols) {
      throw std::out_of_range("Индекс вне массива");
    }

    data->Set(index, item);
  }

  virtual int GetRows() const override { return rows; }
  virtual int GetCols() const override { return cols; }

  virtual T GetDet() const override;
  virtual Matrix<T> *GetInverseMatrix() const override;
  virtual Sequence<T> *SolveSlauGauss(const Sequence<T> *b) const override;
  virtual Sequence<T> *SolveSlauLU(const Sequence<T> *b) const override;
  virtual Matrix<T> *GetL() const override;
  virtual Matrix<T> *GetU() const override;

  virtual double GetFrobeniusNorm() const override {

    double sum = this->data->Reduce(sum_of_squares_func<T>, 0.0);

    return std::sqrt(sum);
  }

  virtual Matrix<T> *MakeGilbert() override {
    if (rows != cols) {
      throw std::invalid_argument("Матрица Гильберта должна быть квадратной");
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        this->SetIJ(i, j, T(1.0) / T(i + j + 1));
      }
    }
    return this;
  }

  virtual Matrix<T> *MakeOnes() override {
    if (rows != cols) {
      throw std::invalid_argument("Единичная матрица обязана быть квадратной");
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        if (i == j)
          this->SetIJ(i, j, T(1.0));
        else
          this->SetIJ(i, j, T(0.0));
      }
    }
    return this;
  }

  virtual Matrix<T> *operator+(const Matrix<T> &rv) const override {
    if (rv.GetCols() != cols || rv.GetRows() != rows) {
      throw std::invalid_argument("Размерности матриц не совпадают");
    }
    Matrix<T> *result = this->CreateEmpty(rows, cols);
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        result->SetIJ(row, col, this->GetIJ(row, col) + rv.GetIJ(row, col));
      }
    }
    return result;
  }

  virtual Matrix<T> *operator-(const Matrix<T> &rv) const override {
    if (rv.GetCols() != cols || rv.GetRows() != rows) {
      throw std::invalid_argument("Размерности матриц не совпадают");
    }
    Matrix<T> *result = this->CreateEmpty(rows, cols);
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        result->SetIJ(row, col, this->GetIJ(row, col) - rv.GetIJ(row, col));
      }
    }
    return result;
  }

  virtual Matrix<T> *operator*(const Matrix<T> &rv) const override {
    if (cols != rv.GetRows()) {
      throw std::invalid_argument("Матрицы несовместимы для умножения");
    }
    int M = this->GetRows();
    int N = this->GetCols();
    int P = rv.GetCols();

    Matrix<T> *result = this->CreateEmpty(M, P);

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < P; j++) {
        T sum = T(0);
        for (int k = 0; k < N; k++) {
          sum += this->GetIJ(i, k) * rv.GetIJ(k, j);
        }
        result->SetIJ(i, j, sum);
      }
    }
    return result;
  }

  virtual Matrix<T> *operator*(const T scalar) const override {
    Matrix<T> *result = this->CreateEmpty(rows, cols);
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        result->SetIJ(row, col, this->GetIJ(row, col) * scalar);
      }
    }
    return result;
  }

  virtual Matrix<T> *Map(T (*mapper)(const T &)) const override {
    Sequence<T> *MappedData = this->data->Map(mapper);
    Matrix<T> *res = this->CreateEmpty(rows, cols);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        int flatIndex = i * cols + j;
        res->SetIJ(i, j, MappedData->Get(flatIndex));
      }
    }
    delete MappedData;
    return res;
  }

  virtual IEnumerator<T> *GetEnumerator() const override {
    return this->data->GetEnumerator();
  }

  virtual T Reduce(T (*reduce_func)(const T &, const T &),
                   const T &start_value) const override {
    return this->data->Reduce(reduce_func, start_value);
  }
};

#include "GaussMethod.h"
#include "LUDecompozition.h"

template <class T>
T BaseMatrix<T>::GetDet() const {
  if (rows != cols) {
    throw std::invalid_argument("Для данной матрицы детерминант не определен");
  }
  GaussMethod<T> mat(this);
  mat.TakeReverse();
  return mat.GetDet();
}

template <class T>
Matrix<T> *BaseMatrix<T>::GetInverseMatrix() const {
  if (rows != cols) {
    throw std::invalid_argument("Для данной матрицы обратная не определена");
  }
  GaussMethod<T> mat(this);
  mat.TakeReverse();
  return mat.GetInverseMatrix();
}

template <class T>
Sequence<T> *BaseMatrix<T>::SolveSlauGauss(const Sequence<T> *b) const {
  if (this->rows != this->cols) {
    throw std::invalid_argument("Матрица системы должна быть квадратной");
  }
  if (b->GetLength() != this->rows) {
    throw std::invalid_argument("Размер вектора не совпадает с матрицей");
  }
  GaussMethod<T> solver(this);
  return solver.Solve(b);
}

template <class T>
Sequence<T> *BaseMatrix<T>::SolveSlauLU(const Sequence<T> *b) const {
  if (this->rows != this->cols) {
    throw std::invalid_argument("Матрица системы должна быть квадратной");
  }
  if (b->GetLength() != this->rows) {
    throw std::invalid_argument("Размер вектора не совпадает с матрицей");
  }
  LUDecompozition<T> solver(this);
  return solver.Solve(b);
}

template <class T>
Matrix<T> *BaseMatrix<T>::GetL() const {
  if (this->rows != this->cols) {
    throw std::invalid_argument("LU-разложение применимо только к квадратным матрицам");
  }
  LUDecompozition<T> lu(this);
  Matrix<T> *L_matrix = lu.GetL();
  Matrix<T> *result = this->CreateEmpty(this->rows, this->cols);
  for (int i = 0; i < this->rows; i++) {
    for (int j = 0; j < this->cols; j++) {
      result->SetIJ(i, j, L_matrix->GetIJ(i, j));
    }
  }
  return result;
}

template <class T>
Matrix<T> *BaseMatrix<T>::GetU() const {
  if (this->rows != this->cols) {
    throw std::invalid_argument("LU-разложение применимо только к квадратным матрицам");
  }
  LUDecompozition<T> lu(this);
  Matrix<T> *U_matrix = lu.GetU();
  Matrix<T> *result = this->CreateEmpty(this->rows, this->cols);
  for (int i = 0; i < this->rows; i++) {
    for (int j = 0; j < this->cols; j++) {
      result->SetIJ(i, j, U_matrix->GetIJ(i, j));
    }
  }
  return result;
}

#endif // BASE_MATRIX_H