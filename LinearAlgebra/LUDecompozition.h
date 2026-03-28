#ifndef LU_DECOMPOZITION_H
#define LU_DECOMPOZITION_H

#include "../Sequences/Sequence.h"
#include "Matrix.h"
#include <stdexcept>

// Класс для LU-разложения квадратной матрицы с частичным выбором ведущего элемента
template <class T> class LUDecompozition {
public:

  // Конструктор принимает исходную матрицу и сразу строит разложение PA = LU
  LUDecompozition(const Matrix<T> *input);
  ~LUDecompozition();

  // Доступ к результатам разложения
  Matrix<T> *GetL() const;
  Matrix<T> *GetU() const;
  Sequence<int> *GetP() const;

  T GetDet() const;

  // Решение СЛАУ Ax = b методом LU-разложения
  Sequence<T> *Solve(const Sequence<T> *b) const;

private:

  int n;           // Размер матрицы
  Matrix<T> *L;    // Нижнетреугольная матрица
  Matrix<T> *U;    // Верхнетреугольная матрица
  Sequence<int> *P; // Вектор перестановок строк

  int swaps; // Количество выполненных перестановок (для знака определителя)

  // Вспомогательные методы разложения
  int neededrow(int skip) const;
  void swaprows(int step);
  void subtraction(int current);
  void decompose();
};

#include "LUDecompozition.tpp"

#endif // LU_DECOMPOZITION_H
