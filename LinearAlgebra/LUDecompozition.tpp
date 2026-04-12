#include "../Sequences/MutableArraySequence.h"
#include "SquareMatrix.h"
#include <cmath>

template <class T> LUDecompozition<T>::LUDecompozition(const Matrix<T> *input) {

  this->n = input->GetRows();

  if (n != input->GetCols()) {
    throw std::invalid_argument("LU-разложение требует квадратную матрицу");
  }

  this->swaps = 0;

  // не получилось импользовать треугольные тк есть свапы, по замерам квадратная
  // со свапами лучше
  this->L = new SquareMatrix<T>(n);
  this->U = new SquareMatrix<T>(n);
  this->P = new MutableArraySequence<int>();

  for (int row = 0; row < n; row++) {
    this->P->append(row);

    for (int col = 0; col < n; col++) {
      this->U->SetIJ(row, col, input->GetIJ(row, col));
      this->L->SetIJ(row, col, (row == col) ? T(1.0) : T(0.0));
    }
  }

  decompose();
}

template <class T> LUDecompozition<T>::~LUDecompozition() {
  delete L;
  delete U;
  delete P;
}

template <class T> int LUDecompozition<T>::NeededRow(int skip) const {
  T mx = std::abs(U->GetIJ(skip, skip));
  int numrow = skip;

  for (int row = skip + 1; row < n; row++) {
    T loc = std::abs(U->GetIJ(row, skip));
    if (loc > mx) {
      mx = loc;
      numrow = row;
    }
  }

  return numrow;
}

template <class T> void LUDecompozition<T>::SwapRows(int step) {
  int we_swap = NeededRow(step);

  if (step == we_swap)
    return;

  swaps++;

  int tempP = P->get(step);
  P->set(step, P->get(we_swap));
  P->set(we_swap, tempP);

  for (int col = step; col < n; col++) {
    T tempU = U->GetIJ(step, col);
    U->SetIJ(step, col, U->GetIJ(we_swap, col));
    U->SetIJ(we_swap, col, tempU);
  }

  for (int col = 0; col < step; col++) {
    T tempL = L->GetIJ(step, col);
    L->SetIJ(step, col, L->GetIJ(we_swap, col));
    L->SetIJ(we_swap, col, tempL);
  }
}

template <class T> void LUDecompozition<T>::subtraction(int current) {
  T pivot = U->GetIJ(current, current);

  if (std::abs(pivot) < 1e-9)
    return;

  for (int row = current + 1; row < n; row++) {
    T mnozh = U->GetIJ(row, current) / pivot;
    L->SetIJ(row, current, mnozh);

    for (int col = current; col < n; col++) {
      U->SetIJ(row, col, U->GetIJ(row, col) - mnozh * U->GetIJ(current, col));
    }
  }
}

template <class T> void LUDecompozition<T>::decompose() {
  for (int row = 0; row < n - 1; row++) {
    SwapRows(row);
    subtraction(row);
  }
}

template <class T> Matrix<T> *LUDecompozition<T>::GetL() const { return L; }

template <class T> Matrix<T> *LUDecompozition<T>::GetU() const { return U; }

template <class T> Sequence<int> *LUDecompozition<T>::GetP() const { return P; }

template <class T> T LUDecompozition<T>::GetDet() const {
  T det = T(1.0);

  for (int row = 0; row < n; row++) {
    det *= U->GetIJ(row, row);
  }

  if (swaps % 2 != 0)
    det = -det;

  return det;
}

template <class T>
Sequence<T> *LUDecompozition<T>::Solve(const Sequence<T> *b) const {

  if (b->GetLength() != n) {
    throw std::invalid_argument(
        "Размер вектора свободных членов не совпадает с матрицей");
  }

  // Прямой ход
  Sequence<T> *y = new MutableArraySequence<T>();
  for (int row = 0; row < n; row++)
    y->append(T(0));

  for (int row = 0; row < n; row++) {
    T sum = b->get(P->get(row));
    for (int col = 0; col < row; col++) {
      sum -= L->GetIJ(row, col) * y->get(col);
    }
    y->set(row, sum);
  }

  // Обратный ход
  Sequence<T> *x = new MutableArraySequence<T>();
  for (int row = 0; row < n; row++)
    x->append(T(0));

  for (int row = n - 1; row >= 0; row--) {
    T sum = y->get(row);
    for (int col = row + 1; col < n; col++) {
      sum -= U->GetIJ(row, col) * x->get(col);
    }
    x->set(row, sum / U->GetIJ(row, row));
  }

  delete y;
  return x;
}
