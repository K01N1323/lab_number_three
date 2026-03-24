#include "LU_decompozition.h"
#include <algorithm> // для std::swap
#include <cmath>

// Конструктор: выделяет память и выполняет декомпозицию
template <class T>
LU_Decomposition<T>::LU_Decomposition(T *matrix_data, int n) {
  this->n = n;
  this->swaps = 0;

  L = new T[n * n];
  U = new T[n * n];
  P = new int[n];
  x = new T[n];
  y = new T[n];

  for (int i = 0; i < n; i++) {
    P[i] = i;

    for (int j = 0; j < n; j++) {
      if (matrix_data != nullptr) {
        U[i * n + j] = matrix_data[i * n + j];
      } else {
        U[i * n + j] = T(0);
      }

      if (i == j) {
        L[i * n + j] = T(1.0);
      } else {
        L[i * n + j] = T(0.0);
      }
    }
  }

  decompose();
}

// Деструктор: освобождает выделенную память
template <class T>
LU_Decomposition<T>::~LU_Decomposition() {
  delete[] L;
  delete[] U;
  delete[] P;
  delete[] x;
  delete[] y;
}

// Ищет максимальный по модулю ведущий элемент в текущем столбце
template <class T>
int LU_Decomposition<T>::neededrow(int skip) const {
  T mx = std::abs(U[skip * n + skip]);
  int numrow = skip;

  for (int i = skip + 1; i < n; i++) {
    T loc = std::abs(U[i * n + skip]);

    if (loc > mx) {
      mx = loc;
      numrow = i;
    }
  }

  return numrow;
}

// Переставляет строки для устойчивости алгоритма
template <class T>
void LU_Decomposition<T>::swaprows(int step) {
  int we_swap = neededrow(step);

  if (step == we_swap) {
    return;
  }

  swaps++;

  std::swap(P[step], P[we_swap]);

  for (int j = step; j < n; j++) {
    std::swap(U[step * n + j], U[we_swap * n + j]);
  }

  for (int j = 0; j < step; j++) {
    std::swap(L[step * n + j], L[we_swap * n + j]);
  }
}

// Вычитает строки и формирует L и U
template <class T>
void LU_Decomposition<T>::subtraction(int current) {
  T pivot = U[current * n + current];

  if (std::abs(pivot) < 1e-9) { // 1e-9 might need cast to T though usually works
    return;
  }

  for (int k = current + 1; k < n; k++) {
    T mnozh = U[k * n + current] / pivot;
    L[k * n + current] = mnozh;

    for (int j = current; j < n; j++) {
      U[k * n + j] -= mnozh * U[current * n + j];
    }
  }
}

// Основной цикл LU-разложения
template <class T>
void LU_Decomposition<T>::decompose() {
  for (int i = 0; i < n - 1; i++) {
    swaprows(i);
    subtraction(i);
  }
}

// Возвращает матрицу L
template <class T>
T *LU_Decomposition<T>::GetL() const { return L; }

// Возвращает матрицу U
template <class T>
T *LU_Decomposition<T>::GetU() const { return U; }

// Возвращает массив перестановок
template <class T>
int *LU_Decomposition<T>::GetP() const { return P; }

// Вычисляет и возвращает определитель матрицы
template <class T>
T LU_Decomposition<T>::GetDet() const {
  T det = T(1.0);

  for (int i = 0; i < n; i++) {
    det *= U[i * n + i];
  }

  if (swaps % 2 != 0) {
    det = -det;
  }

  return det;
}

// Решает СЛАУ: сначала LY = Pb, затем UX = Y
template <class T>
T *LU_Decomposition<T>::Solve(T *b) const {
  for (int i = 0; i < n; i++) {
    y[i] = b[P[i]];

    for (int j = 0; j < i; j++) {
      y[i] -= L[i * n + j] * y[j];
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    x[i] = y[i];

    for (int j = i + 1; j < n; j++) {
      x[i] -= U[i * n + j] * x[j];
    }

    x[i] /= U[i * n + i];
  }

  return x;
}

// Решение СЛАУ без возврата ответа (для тестов скорости)
template <class T>
void LU_Decomposition<T>::SolveForTests(T *b) const {
  for (int i = 0; i < n; i++) {
    y[i] = b[P[i]];

    for (int j = 0; j < i; j++) {
      y[i] -= L[i * n + j] * y[j];
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    x[i] = y[i];

    for (int j = i + 1; j < n; j++) {
      x[i] -= U[i * n + j] * x[j];
    }

    x[i] /= U[i * n + i];
  }
}
