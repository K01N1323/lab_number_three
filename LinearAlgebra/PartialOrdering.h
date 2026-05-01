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

  // Строит транзитивное замыкание отношения методом возведения матрицы в
  // степень
  void BuildTransitiveClosure() {
    int n = Elements->GetLength();

    SquareMatrix<int> *E = new SquareMatrix<int>(n);
    E->MakeOnes();

    Matrix<int> *R = (*GraphMatrix) + (*E);

    for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
        ClosureMatrix->SetIJ(row, col, R->GetIJ(row, col));
      }
    }

    for (int count = 0; count < n - 1; count++) {
      Matrix<int> *LocMull = (*R) * (*ClosureMatrix);

      for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
          ClosureMatrix->SetIJ(row, col, LocMull->GetIJ(row, col));
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