#ifndef MATRIX_H
#define MATRIX_H

#include "../Sequences/IEnumerator.h"
#include "../Sequences/Sequence.h"

// Чисто виртуальный базовый класс для всех видов матриц
template <class T> class Matrix {
public:

  virtual ~Matrix() {}

  virtual const T &GetIJ(int row, int col) const = 0;
  virtual void SetIJ(int row, int col, const T &item) = 0;
  virtual void Set(int index, const T &item) = 0;

  virtual int GetRows() const = 0;
  virtual int GetCols() const = 0;

  virtual T GetDet() const = 0;
  virtual Matrix<T> *GetInverseMatrix() const = 0;

  virtual Sequence<T> *SolveSlauGauss(const Sequence<T> *b) const = 0;
  virtual Sequence<T> *SolveSlauLU(const Sequence<T> *b) const = 0;

  virtual double GetFrobeniusNorm() const = 0;

  virtual Matrix<T> *GetL() const = 0;
  virtual Matrix<T> *GetU() const = 0;

  virtual Matrix<T> *MakeGilbert() = 0;
  virtual Matrix<T> *MakeOnes() = 0;

  virtual Matrix<T> *operator+(const Matrix<T> &rv) const = 0;
  virtual Matrix<T> *operator-(const Matrix<T> &rv) const = 0;
  virtual Matrix<T> *operator*(const Matrix<T> &rv) const = 0;
  virtual Matrix<T> *operator*(const T scalar) const = 0;

  virtual Matrix<T> *Map(T (*mapper)(const T &)) const = 0;

  virtual T Reduce(T (*reduce_func)(const T &, const T &),
                   const T &start_value) const = 0;

  virtual IEnumerator<T> *GetEnumerator() const = 0;

protected:

  virtual Matrix<T> *CreateEmpty(int rows, int cols) const = 0;
};

#endif // MATRIX_H