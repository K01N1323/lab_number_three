#include "LinearAlgebraTests.h"

#include "../LinearAlgebra/DiagonalMatrix.h"
#include "../LinearAlgebra/GaussMethod.h"
#include "../LinearAlgebra/ImmutableMatrix.h"
#include "../LinearAlgebra/LUDecompozition.h"
#include "../LinearAlgebra/LowerTriangleMatrix.h"
#include "../LinearAlgebra/PartialOrdering.h"
#include "../LinearAlgebra/RectangularMatrix.h"
#include "../LinearAlgebra/SquareMatrix.h"
#include "../LinearAlgebra/UpperTriangleMatrix.h"
#include "../Sequences/MutableArraySequence.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

static int TestsPassed = 0;
static int TestsFailed = 0;

static void check(bool condition, const char *TestName) {
  if (condition) {
    cout << "  Успешно: " << TestName << "\n";
    TestsPassed++;
  } else {
    cout << "  Провалено: " << TestName << "\n";
    TestsFailed++;
  }
}

static bool NearEqual(double a, double b, double eps = 1e-9) {
  return fabs(a - b) < eps;
}

void TestSquareMatrix() {
  cout << "\n Тестирование SquareMatrix \n";

  {
    SquareMatrix<double> m(3);
    m.SetIJ(0, 0, 1.0);
    m.SetIJ(1, 2, 5.0);
    m.SetIJ(2, 1, -3.0);
    check(NearEqual(m.GetIJ(0, 0), 1.0), "Запись и чтение элемента [0,0]");
    check(NearEqual(m.GetIJ(1, 2), 5.0), "Запись и чтение элемента [1,2]");
    check(NearEqual(m.GetIJ(0, 1), 0.0), "Нулевой элемент по умолчанию");
    check(m.GetRows() == 3 && m.GetCols() == 3, "Размеры матрицы корректны");
  }

  {
    SquareMatrix<double> m(2);
    bool threw = false;
    try {
      m.GetIJ(2, 0);
    } catch (const out_of_range &) {
      threw = true;
    }
    check(threw, "Исключение при обращении вне границ");
  }

  {
    SquareMatrix<double> m(3);
    m.MakeOnes();
    bool ok = true;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (!NearEqual(m.GetIJ(i, j), (i == j) ? 1.0 : 0.0))
          ok = false;
    check(ok, "MakeOnes формирует единичную матрицу");
  }

  {
    // [[1,2],[3,4]] + [[5,6],[7,8]] = [[6,8],[10,12]]
    SquareMatrix<double> a(2), b(2);
    a.SetIJ(0, 0, 1);
    a.SetIJ(0, 1, 2);
    a.SetIJ(1, 0, 3);
    a.SetIJ(1, 1, 4);
    b.SetIJ(0, 0, 5);
    b.SetIJ(0, 1, 6);
    b.SetIJ(1, 0, 7);
    b.SetIJ(1, 1, 8);
    Matrix<double> *c = a + b;
    check(NearEqual(c->GetIJ(0, 0), 6.0) && NearEqual(c->GetIJ(1, 1), 12.0),
          "Сложение матриц");
    delete c;
  }

  {
    // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
    SquareMatrix<double> a(2), b(2);
    a.SetIJ(0, 0, 1);
    a.SetIJ(0, 1, 2);
    a.SetIJ(1, 0, 3);
    a.SetIJ(1, 1, 4);
    b.SetIJ(0, 0, 5);
    b.SetIJ(0, 1, 6);
    b.SetIJ(1, 0, 7);
    b.SetIJ(1, 1, 8);
    Matrix<double> *c = a * b;
    check(NearEqual(c->GetIJ(0, 0), 19) && NearEqual(c->GetIJ(0, 1), 22) &&
              NearEqual(c->GetIJ(1, 0), 43) && NearEqual(c->GetIJ(1, 1), 50),
          "Умножение матриц");
    delete c;
  }

  {
    // det([[1,2],[3,4]]) = -2
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 1);
    m.SetIJ(0, 1, 2);
    m.SetIJ(1, 0, 3);
    m.SetIJ(1, 1, 4);
    check(NearEqual(m.GetDet(), -2.0, 1e-9), "Определитель 2x2");
  }

  {
    // норма Фробениуса единичной матрицы 3x3 = sqrt(3)
    SquareMatrix<double> m(3);
    m.MakeOnes();
    check(NearEqual(m.GetFrobeniusNorm(), sqrt(3.0), 1e-9),
          "Норма Фробениуса единичной матрицы");
  }

  {
    // x + y = 3, 2x - y = 0  =>  x=1, y=2
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 1);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 2);
    m.SetIJ(1, 1, -1);
    MutableArraySequence<double> b;
    b.append(3.0);
    b.append(0.0);
    Sequence<double> *x = m.SolveSlauGauss(&b);
    check(NearEqual(x->get(0), 1.0) && NearEqual(x->get(1), 2.0),
          "Решение СЛАУ методом Гаусса");
    delete x;
  }

  {
    // то же уравнение, метод LU
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 1);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 2);
    m.SetIJ(1, 1, -1);
    MutableArraySequence<double> b;
    b.append(3.0);
    b.append(0.0);
    Sequence<double> *x = m.SolveSlauLU(&b);
    check(NearEqual(x->get(0), 1.0) && NearEqual(x->get(1), 2.0),
          "Решение СЛАУ методом LU");
    delete x;
  }

  {
    // inv([[1,2],[3,4]]) = [[-2,1],[1.5,-0.5]]
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 1);
    m.SetIJ(0, 1, 2);
    m.SetIJ(1, 0, 3);
    m.SetIJ(1, 1, 4);
    Matrix<double> *inv = m.GetInverseMatrix();
    check(NearEqual(inv->GetIJ(0, 0), -2.0) && NearEqual(inv->GetIJ(1, 0), 1.5),
          "Обратная матрица");
    delete inv;
  }
}

void TestRectangularMatrix() {
  cout << "\nТестирование RectangularMatrix \n";

  {
    RectangularMatrix<int> m(2, 3);
    m.SetIJ(0, 2, 7);
    m.SetIJ(1, 0, 4);
    check(m.GetIJ(0, 2) == 7 && m.GetIJ(1, 0) == 4,
          "Запись и чтение элементов");
    check(m.GetIJ(0, 0) == 0, "Нулевой элемент по умолчанию");
  }

  {
    RectangularMatrix<double> m(4, 7);
    check(m.GetRows() == 4 && m.GetCols() == 7, "Размеры матрицы корректны");
  }
}

void TestDiagonalMatrix() {
  cout << "\n Тестирование DiagonalMatrix \n";

  {
    DiagonalMatrix<double> m(3);
    m.SetIJ(0, 0, 2.0);
    m.SetIJ(1, 1, 5.0);
    m.SetIJ(2, 2, -1.0);
    check(NearEqual(m.GetIJ(0, 0), 2.0) && NearEqual(m.GetIJ(1, 1), 5.0),
          "Чтение диагональных элементов");
    check(NearEqual(m.GetIJ(0, 1), 0.0) && NearEqual(m.GetIJ(2, 0), 0.0),
          "Внедиагональные элементы равны нулю");
  }

  {
    DiagonalMatrix<double> m(3);
    bool threw = false;
    try {
      m.SetIJ(0, 1, 5.0);
    } catch (const invalid_argument &) {
      threw = true;
    }
    check(threw, "Исключение при записи вне диагонали");
  }

  {
    // det диагональной матрицы = произведение диагоналей = 2*3*4 = 24
    DiagonalMatrix<double> m(3);
    m.SetIJ(0, 0, 2);
    m.SetIJ(1, 1, 3);
    m.SetIJ(2, 2, 4);
    check(NearEqual(m.GetDet(), 24.0, 1e-7), "Определитель диагональной матрицы");
  }
}

void TestUpperTriangularMatrix() {
  cout << "\n Тестирование UpperTriangularMatrix \n";

  {
    UpperTriangularMatrix<double> m(3);
    m.SetIJ(0, 0, 1);
    m.SetIJ(0, 1, 2);
    m.SetIJ(0, 2, 3);
    m.SetIJ(1, 1, 4);
    m.SetIJ(1, 2, 5);
    m.SetIJ(2, 2, 6);
    check(NearEqual(m.GetIJ(1, 1), 4.0) && NearEqual(m.GetIJ(0, 2), 3.0),
          "Чтение надиагональных элементов");
    check(NearEqual(m.GetIJ(1, 0), 0.0) && NearEqual(m.GetIJ(2, 0), 0.0),
          "Поддиагональные элементы равны нулю");
  }

  {
    UpperTriangularMatrix<double> m(3);
    bool threw = false;
    try {
      m.SetIJ(2, 0, 1.0);
    } catch (const invalid_argument &) {
      threw = true;
    }
    check(threw, "Исключение при записи ниже диагонали");
  }
}

void TestLowerTriangularMatrix() {
  cout << "\n Тестирование LowerTriangularMatrix \n";

  {
    LowerTriangularMatrix<double> m(3);
    m.SetIJ(0, 0, 1);
    m.SetIJ(1, 0, 2);
    m.SetIJ(1, 1, 3);
    m.SetIJ(2, 2, 4);
    check(NearEqual(m.GetIJ(1, 0), 2.0) && NearEqual(m.GetIJ(2, 2), 4.0),
          "Чтение поддиагональных элементов");
    check(NearEqual(m.GetIJ(0, 1), 0.0) && NearEqual(m.GetIJ(0, 2), 0.0),
          "Надиагональные элементы равны нулю");
  }

  {
    LowerTriangularMatrix<double> m(3);
    bool threw = false;
    try {
      m.SetIJ(0, 2, 1.0);
    } catch (const invalid_argument &) {
      threw = true;
    }
    check(threw, "Исключение при записи выше диагонали");
  }
}

void TestImmutableMatrix() {
  cout << "\n Тестирование ImmutableMatrix \n";

  {
    SquareMatrix<double> *base = new SquareMatrix<double>(2);
    base->SetIJ(0, 0, 7.0);
    base->SetIJ(1, 1, 3.0);
    ImmutableMatrix<double> im(base);
    check(NearEqual(im.GetIJ(0, 0), 7.0) && NearEqual(im.GetIJ(1, 1), 3.0),
          "Чтение элементов через ImmutableMatrix");
  }

  {
    ImmutableMatrix<double> im(new SquareMatrix<double>(2));
    bool threw = false;
    try {
      im.SetIJ(0, 0, 1.0);
    } catch (const logic_error &) {
      threw = true;
    }
    check(threw, "Исключение при попытке изменить неизменяемую матрицу");
  }
}

void TestGaussMethod() {
  cout << "\n Тестирование GaussMethod\n";

  {
    // 3x + y = 9, x + 2y = 8  =>  x=2, y=3
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    MutableArraySequence<double> b;
    b.append(9.0);
    b.append(8.0);
    GaussMethod<double> g(&m);
    Sequence<double> *x = g.Solve(&b);
    check(NearEqual(x->get(0), 2.0) && NearEqual(x->get(1), 3.0), "Решение СЛАУ");
    delete x;
  }

  {
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    GaussMethod<double> g(&m);
    g.TakeReverse();
    check(NearEqual(g.GetDet(), 5.0, 1e-9), "Вычисление определителя");
  }

  {
    // inv([[2,1],[1,3]]) = [[3/5,-1/5],[-1/5,2/5]]
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 2);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 3);
    GaussMethod<double> g(&m);
    g.TakeReverse();
    Matrix<double> *inv = g.GetInverseMatrix();
    check(NearEqual(inv->GetIJ(0, 0), 3.0 / 5) &&
              NearEqual(inv->GetIJ(1, 0), -1.0 / 5),
          "Обратная матрица");
    delete inv;
  }
}

void TestLUDecompozition() {
  cout << "\n Тестирование LUDecompozition \n";

  {
    // 3x + y = 9, x + 2y = 8  =>  x=2, y=3
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    MutableArraySequence<double> b;
    b.append(9.0);
    b.append(8.0);
    LUDecompozition<double> lu(&m);
    Sequence<double> *x = lu.Solve(&b);
    check(NearEqual(x->get(0), 2.0) && NearEqual(x->get(1), 3.0), "Решение СЛАУ");
    delete x;
  }

  {
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    LUDecompozition<double> lu(&m);
    check(NearEqual(lu.GetDet(), 5.0, 1e-9), "Вычисление определителя");
  }

  {
    // Для единичной матрицы L и U должны быть единичными
    SquareMatrix<double> m(3);
    m.MakeOnes();
    LUDecompozition<double> lu(&m);
    Matrix<double> *L = lu.GetL();
    Matrix<double> *U = lu.GetU();
    bool lOk = true, uOk = true;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double e = (i == j) ? 1.0 : 0.0;
        if (!NearEqual(L->GetIJ(i, j), e))
          lOk = false;
        if (!NearEqual(U->GetIJ(i, j), e))
          uOk = false;
      }
    check(lOk, "Матрица L для единичной матрицы корректна");
    check(uOk, "Матрица U для единичной матрицы корректна");
  }

  {
    RectangularMatrix<double> m(2, 3);
    bool threw = false;
    try {
      LUDecompozition<double> lu(&m);
    } catch (const invalid_argument &) {
      threw = true;
    }
    check(threw, "Исключение для прямоугольной матрицы");
  }
}

void TestPartialOrdering() {
  cout << "\n Тестирование PartialOrdering \n";

  {
    // Отношение: 1 <= 2, 2 <= 3
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    pairs.append({2, 3});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 2) && po.IsLessOrEqual(2, 3),
          "Прямые дуги отношения");
  }

  {
    // Транзитивно: 1 <= 3 должно выполняться
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    pairs.append({2, 3});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 3), "Транзитивная дуга 1 <= 3");
  }

  {
    // Рефлексивность: каждый элемент <= себя
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 1) && po.IsLessOrEqual(2, 2),
          "Рефлексивность отношения");
  }

  {
    // Антисимметрия: 3 не должен быть <= 1
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    pairs.append({2, 3});
    PartialOrdering<int> po(&pairs);
    check(!po.IsLessOrEqual(3, 1) && !po.IsLessOrEqual(2, 1),
          "Отсутствие обратного отношения");
  }

  {
    // Цепь {1,2},{2,3} — транзитивное замыкание даёт 6 дуг (с рефлексивными)
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    pairs.append({2, 3});
    PartialOrdering<int> po(&pairs);
    Sequence<Pair<int>> *edges = po.GetMaterializedEdges();
    int count = edges->GetLength();
    delete edges;
    check(count == 6,
          "Количество материализованных дуг после транзитивного замыкания");
  }

  {
    // Длинная цепь: 1->2->3->4->5, транзитивно 1<=5
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    pairs.append({2, 3});
    pairs.append({3, 4});
    pairs.append({4, 5});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 5) && po.IsLessOrEqual(1, 4),
          "Транзитивность в длинной цепи");
    check(!po.IsLessOrEqual(5, 1),
          "Отсутствие обратного отношения в длинной цепи");
  }

  {
    MutableArraySequence<Pair<int>> pairs;
    pairs.append({1, 2});
    PartialOrdering<int> po(&pairs);
    bool threw = false;
    try {
      po.IsLessOrEqual(1, 99);
    } catch (const invalid_argument &) {
      threw = true;
    }
    check(threw, "Исключение при неизвестном элементе");
  }
}

void RunAllTests() {
  TestsPassed = 0;
  TestsFailed = 0;

  TestSquareMatrix();
  TestRectangularMatrix();
  TestDiagonalMatrix();
  TestUpperTriangularMatrix();
  TestLowerTriangularMatrix();
  TestImmutableMatrix();
  TestGaussMethod();
  TestLUDecompozition();
  TestPartialOrdering();

  cout << " Результаты тестов \n";
  cout << "Пройдено:  " << TestsPassed << "\n";
  cout << "Провалено: " << TestsFailed << "\n";

  if (TestsFailed == 0) {
    cout << "Все тесты успешно пройдены!\n";
  } else {
    cout << "Были обнаружены ошибки.\n";
  }
}
