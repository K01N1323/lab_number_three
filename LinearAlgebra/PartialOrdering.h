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
// Транзитивное замыкание строим при инициализации
template <class T> class PartialOrdering {
private:
  Sequence<T> *Elements;
  SquareMatrix<int> *GraphMatrix;
  SquareMatrix<int> *ClosureMatrix;

  int GetIndex(const T &item) const {
    for (int index = 0; index < Elements->GetLength(); index++) {
      if (Elements->Get(index) == item) {
        return index;
      }
    }
    return -1;
  }

  // Строит транзитивное замыкание отношения алгоритмом Флойда-Уоршелла O(n^3)
  void BuildTransitiveClosure() {
    int n = Elements->GetLength();

    // Копируем граф в ClosureMatrix и добавляем рефлексивность (диагональ)
    for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
        ClosureMatrix->SetIJ(row, col, GraphMatrix->GetIJ(row, col));
      }
      ClosureMatrix->SetIJ(row, row, 1);
    }

    // Алгоритм Флойда-Уоршелла
    for (int k = 0; k < n; k++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (ClosureMatrix->GetIJ(i, k) > 0 && ClosureMatrix->GetIJ(k, j) > 0) {
            ClosureMatrix->SetIJ(i, j, 1);
          }
        }
      }
    }
  }

public:
  PartialOrdering(Sequence<Pair<T>> *pairs) {
    Elements = new MutableArraySequence<T>();

    // Собираем уникальные элементы из всех пар
    for (int number = 0; number < pairs->GetLength(); number++) {
      T FirstItem = pairs->Get(number).First;
      T SecondItem = pairs->Get(number).Second;

      bool IsFirstFound = false;
      bool IsSecondFound = false;

      for (int index = 0; index < Elements->GetLength(); index++) {
        if (Elements->Get(index) == FirstItem)
          IsFirstFound = true;
        if (Elements->Get(index) == SecondItem)
          IsSecondFound = true;
      }

      if (!IsFirstFound)
        Elements->Append(FirstItem);
      if (!IsSecondFound)
        Elements->Append(SecondItem);
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

  bool IsLessOrEqual(const T &x1, const T &x2) const {
    int row = GetIndex(x1);
    int col = GetIndex(x2);

    if (row == -1 || col == -1) {
      throw std::invalid_argument(
          "Один или оба элемента отсутствуют в частичном порядке");
    }

    return ClosureMatrix->GetIJ(row, col) > 0;
  }

  Sequence<Pair<T>> *GetMaterializedEdges() const {
    int n = Elements->GetLength();

    Sequence<Pair<T>> *result = new MutableArraySequence<Pair<T>>();

    for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
        if (ClosureMatrix->GetIJ(row, col) > 0) {
          Pair<T> NewPair;
          NewPair.First = Elements->Get(row);
          NewPair.Second = Elements->Get(col);
          result->Append(NewPair);
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