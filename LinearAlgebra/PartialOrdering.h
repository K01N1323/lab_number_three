#ifndef PARTIAL_ORDERING_H
#define PARTIAL_ORDERING_H

#include "../Sequences/MutableArraySequence.h"
#include "../Sequences/Sequence.h"
#include "SquareMatrix.h"

template <typename T> struct Pair {
  T First;
  T Second;
};

// Частичный порядок на конечном множестве элементов,
// заданный списком пар. Транзитивное замыкание строится при инициализации.
template <class T> class PartialOrdering {
private:

  Sequence<T> *Elements;
  SquareMatrix<int> *GraphMatrix;
  SquareMatrix<int> *ClosureMatrix;

  // Возвращает индекс элемента в списке Elements, или -1 если не найден
  int GetIndex(const T &item) const {
    for (int index = 0; index < Elements->GetLength(); index++) {
      if (Elements->Get(index) == item) {
        return index;
      }
    }
    return -1;
  }

  // Строит транзитивное замыкание отношения методом возведения матрицы в степень
  void BuildTransitiveClosure() {
    int n = Elements->GetLength();

    SquareMatrix<int> *E = new SquareMatrix<int>(n);
    E->MakeOnes();

    Matrix<int> *R = (*GraphMatrix) + (*E);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        ClosureMatrix->SetIJ(i, j, R->GetIJ(i, j));
      }
    }

    for (int count = 0; count < n - 1; count++) {
      Matrix<int> *LocMull = (*R) * (*ClosureMatrix);

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          ClosureMatrix->SetIJ(i, j, LocMull->GetIJ(i, j));
        }
      }

      delete LocMull;
    }

    // Нормируем: любое ненулевое значение заменяем на 1
    for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
        if (ClosureMatrix->GetIJ(row, col) > 0) {
          ClosureMatrix->SetIJ(row, col, 1);
        }
      }
    }

    delete E;
    delete R;
  }

public:

  PartialOrdering(Sequence<Pair<T>> *pairs) {
    Elements = new MutableArraySequence<T>();

    // Собираем уникальные элементы из всех пар
    for (int number = 0; number < pairs->GetLength(); number++) {
      T firstItem = pairs->Get(number).First;
      T secondItem = pairs->Get(number).Second;

      bool isFirstFound = false;
      bool isSecondFound = false;

      for (int index = 0; index < Elements->GetLength(); index++) {
        if (Elements->Get(index) == firstItem)  isFirstFound = true;
        if (Elements->Get(index) == secondItem) isSecondFound = true;
      }

      if (!isFirstFound)  Elements->Append(firstItem);
      if (!isSecondFound) Elements->Append(secondItem);
    }

    int n = Elements->GetLength();

    GraphMatrix = new SquareMatrix<int>(n);
    ClosureMatrix = new SquareMatrix<int>(n);

    // Заполняем матрицу смежности по переданным парам
    for (int number = 0; number < pairs->GetLength(); number++) {
      int row = GetIndex(pairs->Get(number).First);
      int col = GetIndex(pairs->Get(number).Second);

      if (row != -1 && col != -1) {
        GraphMatrix->SetIJ(row, col, 1);
      }
    }

    BuildTransitiveClosure();
  }

  // Возвращает true, если x1 <= x2 в данном частичном порядке
  bool IsLessOrEqual(const T &x1, const T &x2) const {
    int row = GetIndex(x1);
    int col = GetIndex(x2);

    if (row == -1 || col == -1) {
      throw std::invalid_argument(
          "Один или оба элемента отсутствуют в частичном порядке");
    }

    return ClosureMatrix->GetIJ(row, col) > 0;
  }

  // Возвращает все пары, для которых выполняется отношение порядка
  Sequence<Pair<T>> *GetMaterializedEdges() const {
    int n = Elements->GetLength();

    Sequence<Pair<T>> *result = new MutableArraySequence<Pair<T>>();

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (ClosureMatrix->GetIJ(i, j) > 0) {
          Pair<T> newPair;
          newPair.First = Elements->Get(i);
          newPair.Second = Elements->Get(j);
          result->Append(newPair);
        }
      }
    }

    return result;
  }

  ~PartialOrdering() {
    delete Elements;
    delete GraphMatrix;
    delete ClosureMatrix;
  }
};

#endif // PARTIAL_ORDERING_H