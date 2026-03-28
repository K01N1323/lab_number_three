#include "LUDecompozition.h"
#include "../Sequences/MutableArraySequence.h"
#include "SquareMatrix.h"
#include <cmath>

template <class T>
LUDecompozition<T>::LUDecompozition(const Matrix<T> *input) {

  this->n = input->GetRows();

  if (n != input->GetCols()) {
    throw std::invalid_argument("LU-разложение требует квадратную матрицу");
  }

  this->swaps = 0;

  // Используем SquareMatrix: треугольная структура нарушается при перестановках
  this->L = new SquareMatrix<T>(n);
  this->U = new SquareMatrix<T>(n);
  this->P = new MutableArraySequence<int>();

  // Копируем матрицу в U, L инициализируем единичной, P — тождественной
  for (int i = 0; i < n; i++) {
    this->P->Append(i);

    for (int j = 0; j < n; j++) {
      this->U->SetIJ(i, j, input->GetIJ(i, j));
      this->L->SetIJ(i, j, (i == j) ? T(1.0) : T(0.0));
    }
  }

  decompose();
}

template <class T> LUDecompozition<T>::~LUDecompozition() {
  delete L;
  delete U;
  delete P;
}

// Возвращает индекс строки с максимальным по модулю элементом в столбце skip
template <class T> int LUDecompozition<T>::neededrow(int skip) const {
  T mx = std::abs(U->GetIJ(skip, skip));
  int numrow = skip;

  for (int i = skip + 1; i < n; i++) {
    T loc = std::abs(U->GetIJ(i, skip));
    if (loc > mx) {
      mx = loc;
      numrow = i;
    }
  }

  return numrow;
}

// Меняет строки в U и L местами (частичный выбор ведущего элемента)
template <class T> void LUDecompozition<T>::swaprows(int step) {
  int we_swap = neededrow(step);

  if (step == we_swap)
    return;

  swaps++;

  int tempP = P->Get(step);
  P->Set(step, P->Get(we_swap));
  P->Set(we_swap, tempP);

  for (int j = step; j < n; j++) {
    T tempU = U->GetIJ(step, j);
    U->SetIJ(step, j, U->GetIJ(we_swap, j));
    U->SetIJ(we_swap, j, tempU);
  }

  for (int j = 0; j < step; j++) {
    T tempL = L->GetIJ(step, j);
    L->SetIJ(step, j, L->GetIJ(we_swap, j));
    L->SetIJ(we_swap, j, tempL);
  }
}

// Вычитает текущую строку из нижележащих, заполняя L и обнуляя столбец в U
template <class T> void LUDecompozition<T>::subtraction(int current) {
  T pivot = U->GetIJ(current, current);

  if (std::abs(pivot) < 1e-9)
    return;

  for (int k = current + 1; k < n; k++) {
    T mnozh = U->GetIJ(k, current) / pivot;
    L->SetIJ(k, current, mnozh);

    for (int j = current; j < n; j++) {
      U->SetIJ(k, j, U->GetIJ(k, j) - mnozh * U->GetIJ(current, j));
    }
  }
}

// Основной цикл разложения PA = LU
template <class T> void LUDecompozition<T>::decompose() {
  for (int i = 0; i < n - 1; i++) {
    swaprows(i);
    subtraction(i);
  }
}

template <class T> Matrix<T> *LUDecompozition<T>::GetL() const { return L; }

template <class T> Matrix<T> *LUDecompozition<T>::GetU() const { return U; }

template <class T> Sequence<int> *LUDecompozition<T>::GetP() const { return P; }

// Произведение диагональных элементов U с учётом знака перестановок
template <class T> T LUDecompozition<T>::GetDet() const {
  T det = T(1.0);

  for (int i = 0; i < n; i++) {
    det *= U->GetIJ(i, i);
  }

  if (swaps % 2 != 0)
    det = -det;

  return det;
}

// Решает Ax = b через прямой ход Ly = Pb и обратный ход Ux = y
template <class T>
Sequence<T> *LUDecompozition<T>::Solve(const Sequence<T> *b) const {

  if (b->GetLength() != n) {
    throw std::invalid_argument(
        "Размер вектора свободных членов не совпадает с матрицей");
  }

  // Прямой ход
  Sequence<T> *y = new MutableArraySequence<T>();
  for (int i = 0; i < n; i++)
    y->Append(T(0));

  for (int i = 0; i < n; i++) {
    T sum = b->Get(P->Get(i));
    for (int j = 0; j < i; j++) {
      sum -= L->GetIJ(i, j) * y->Get(j);
    }
    y->Set(i, sum);
  }

  // Обратный ход
  Sequence<T> *x = new MutableArraySequence<T>();
  for (int i = 0; i < n; i++)
    x->Append(T(0));

  for (int i = n - 1; i >= 0; i--) {
    T sum = y->Get(i);
    for (int j = i + 1; j < n; j++) {
      sum -= U->GetIJ(i, j) * x->Get(j);
    }
    x->Set(i, sum / U->GetIJ(i, i));
  }

  delete y;
  return x;
}
