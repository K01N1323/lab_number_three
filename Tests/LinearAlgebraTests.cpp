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

static int tests_passed = 0;
static int tests_failed = 0;

static void check(bool condition, const char *test_name) {
  if (condition) {
    cout << "  Успешно: " << test_name << "\n";
    tests_passed++;
  } else {
    cout << "  Провалено: " << test_name << "\n";
    tests_failed++;
  }
}

static bool near_eq(double a, double b, double eps = 1e-9) {
  return fabs(a - b) < eps;
}

void test_square_matrix() {
  cout << "\n Тестирование SquareMatrix \n";

  {
    SquareMatrix<double> m(3);
    m.SetIJ(0, 0, 1.0);
    m.SetIJ(1, 2, 5.0);
    m.SetIJ(2, 1, -3.0);
    check(near_eq(m.GetIJ(0, 0), 1.0), "Запись и чтение элемента [0,0]");
    check(near_eq(m.GetIJ(1, 2), 5.0), "Запись и чтение элемента [1,2]");
    check(near_eq(m.GetIJ(0, 1), 0.0), "Нулевой элемент по умолчанию");
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
        if (!near_eq(m.GetIJ(i, j), (i == j) ? 1.0 : 0.0))
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
    check(near_eq(c->GetIJ(0, 0), 6.0) && near_eq(c->GetIJ(1, 1), 12.0),
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
    check(near_eq(c->GetIJ(0, 0), 19) && near_eq(c->GetIJ(0, 1), 22) &&
              near_eq(c->GetIJ(1, 0), 43) && near_eq(c->GetIJ(1, 1), 50),
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
    check(near_eq(m.GetDet(), -2.0, 1e-9), "Определитель 2x2");
  }

  {
    // норма Фробениуса единичной матрицы 3x3 = sqrt(3)
    SquareMatrix<double> m(3);
    m.MakeOnes();
    check(near_eq(m.GetFrobeniusNorm(), sqrt(3.0), 1e-9),
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
    b.Append(3.0);
    b.Append(0.0);
    Sequence<double> *x = m.SolveSlauGauss(&b);
    check(near_eq(x->Get(0), 1.0) && near_eq(x->Get(1), 2.0),
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
    b.Append(3.0);
    b.Append(0.0);
    Sequence<double> *x = m.SolveSlauLU(&b);
    check(near_eq(x->Get(0), 1.0) && near_eq(x->Get(1), 2.0),
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
    check(near_eq(inv->GetIJ(0, 0), -2.0) && near_eq(inv->GetIJ(1, 0), 1.5),
          "Обратная матрица");
    delete inv;
  }
}

void test_rectangular_matrix() {
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

void test_diagonal_matrix() {
  cout << "\n Тестирование DiagonalMatrix \n";

  {
    DiagonalMatrix<double> m(3);
    m.SetIJ(0, 0, 2.0);
    m.SetIJ(1, 1, 5.0);
    m.SetIJ(2, 2, -1.0);
    check(near_eq(m.GetIJ(0, 0), 2.0) && near_eq(m.GetIJ(1, 1), 5.0),
          "Чтение диагональных элементов");
    check(near_eq(m.GetIJ(0, 1), 0.0) && near_eq(m.GetIJ(2, 0), 0.0),
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
    check(near_eq(m.GetDet(), 24.0, 1e-7), "Определитель диагональной матрицы");
  }
}

void test_upper_triangular_matrix() {
  cout << "\n Тестирование UpperTriangularMatrix \n";

  {
    UpperTriangularMatrix<double> m(3);
    m.SetIJ(0, 0, 1);
    m.SetIJ(0, 1, 2);
    m.SetIJ(0, 2, 3);
    m.SetIJ(1, 1, 4);
    m.SetIJ(1, 2, 5);
    m.SetIJ(2, 2, 6);
    check(near_eq(m.GetIJ(1, 1), 4.0) && near_eq(m.GetIJ(0, 2), 3.0),
          "Чтение надиагональных элементов");
    check(near_eq(m.GetIJ(1, 0), 0.0) && near_eq(m.GetIJ(2, 0), 0.0),
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

void test_lower_triangular_matrix() {
  cout << "\n Тестирование LowerTriangularMatrix \n";

  {
    LowerTriangularMatrix<double> m(3);
    m.SetIJ(0, 0, 1);
    m.SetIJ(1, 0, 2);
    m.SetIJ(1, 1, 3);
    m.SetIJ(2, 2, 4);
    check(near_eq(m.GetIJ(1, 0), 2.0) && near_eq(m.GetIJ(2, 2), 4.0),
          "Чтение поддиагональных элементов");
    check(near_eq(m.GetIJ(0, 1), 0.0) && near_eq(m.GetIJ(0, 2), 0.0),
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

void test_immutable_matrix() {
  cout << "\n Тестирование ImmutableMatrix \n";

  {
    SquareMatrix<double> *base = new SquareMatrix<double>(2);
    base->SetIJ(0, 0, 7.0);
    base->SetIJ(1, 1, 3.0);
    ImmutableMatrix<double> im(base);
    check(near_eq(im.GetIJ(0, 0), 7.0) && near_eq(im.GetIJ(1, 1), 3.0),
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

void test_gauss_method() {
  cout << "\n Тестирование GaussMethod\n";

  {
    // 3x + y = 9, x + 2y = 8  =>  x=2, y=3
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    MutableArraySequence<double> b;
    b.Append(9.0);
    b.Append(8.0);
    GaussMethod<double> g(&m);
    Sequence<double> *x = g.Solve(&b);
    check(near_eq(x->Get(0), 2.0) && near_eq(x->Get(1), 3.0), "Решение СЛАУ");
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
    check(near_eq(g.GetDet(), 5.0, 1e-9), "Вычисление определителя");
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
    check(near_eq(inv->GetIJ(0, 0), 3.0 / 5) &&
              near_eq(inv->GetIJ(1, 0), -1.0 / 5),
          "Обратная матрица");
    delete inv;
  }
}

void test_lu_decompozition() {
  cout << "\n Тестирование LUDecompozition \n";

  {
    // 3x + y = 9, x + 2y = 8  =>  x=2, y=3
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    MutableArraySequence<double> b;
    b.Append(9.0);
    b.Append(8.0);
    LUDecompozition<double> lu(&m);
    Sequence<double> *x = lu.Solve(&b);
    check(near_eq(x->Get(0), 2.0) && near_eq(x->Get(1), 3.0), "Решение СЛАУ");
    delete x;
  }

  {
    SquareMatrix<double> m(2);
    m.SetIJ(0, 0, 3);
    m.SetIJ(0, 1, 1);
    m.SetIJ(1, 0, 1);
    m.SetIJ(1, 1, 2);
    LUDecompozition<double> lu(&m);
    check(near_eq(lu.GetDet(), 5.0, 1e-9), "Вычисление определителя");
  }

  {
    // Для единичной матрицы L и U должны быть единичными
    SquareMatrix<double> m(3);
    m.MakeOnes();
    LUDecompozition<double> lu(&m);
    Matrix<double> *L = lu.GetL();
    Matrix<double> *U = lu.GetU();
    bool l_ok = true, u_ok = true;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double e = (i == j) ? 1.0 : 0.0;
        if (!near_eq(L->GetIJ(i, j), e))
          l_ok = false;
        if (!near_eq(U->GetIJ(i, j), e))
          u_ok = false;
      }
    check(l_ok, "Матрица L для единичной матрицы корректна");
    check(u_ok, "Матрица U для единичной матрицы корректна");
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

void test_partial_ordering() {
  cout << "\n Тестирование PartialOrdering \n";

  {
    // Отношение: 1 <= 2, 2 <= 3
    MutableArraySequence<Pair<int>> pairs;
    pairs.Append({1, 2});
    pairs.Append({2, 3});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 2) && po.IsLessOrEqual(2, 3),
          "Прямые дуги отношения");
  }

  {
    // Транзитивно: 1 <= 3 должно выполняться
    MutableArraySequence<Pair<int>> pairs;
    pairs.Append({1, 2});
    pairs.Append({2, 3});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 3), "Транзитивная дуга 1 <= 3");
  }

  {
    // Рефлексивность: каждый элемент <= себя
    MutableArraySequence<Pair<int>> pairs;
    pairs.Append({1, 2});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 1) && po.IsLessOrEqual(2, 2),
          "Рефлексивность отношения");
  }

  {
    // Антисимметрия: 3 не должен быть <= 1
    MutableArraySequence<Pair<int>> pairs;
    pairs.Append({1, 2});
    pairs.Append({2, 3});
    PartialOrdering<int> po(&pairs);
    check(!po.IsLessOrEqual(3, 1) && !po.IsLessOrEqual(2, 1),
          "Отсутствие обратного отношения");
  }

  {
    // Цепь {1,2},{2,3} — транзитивное замыкание даёт 6 дуг (с рефлексивными)
    MutableArraySequence<Pair<int>> pairs;
    pairs.Append({1, 2});
    pairs.Append({2, 3});
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
    pairs.Append({1, 2});
    pairs.Append({2, 3});
    pairs.Append({3, 4});
    pairs.Append({4, 5});
    PartialOrdering<int> po(&pairs);
    check(po.IsLessOrEqual(1, 5) && po.IsLessOrEqual(1, 4),
          "Транзитивность в длинной цепи");
    check(!po.IsLessOrEqual(5, 1),
          "Отсутствие обратного отношения в длинной цепи");
  }

  {
    MutableArraySequence<Pair<int>> pairs;
    pairs.Append({1, 2});
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

void run_all_tests() {
  tests_passed = 0;
  tests_failed = 0;

  test_square_matrix();
  test_rectangular_matrix();
  test_diagonal_matrix();
  test_upper_triangular_matrix();
  test_lower_triangular_matrix();
  test_immutable_matrix();
  test_gauss_method();
  test_lu_decompozition();
  test_partial_ordering();

  cout << " Результаты тестов \n";
  cout << "Пройдено:  " << tests_passed << "\n";
  cout << "Провалено: " << tests_failed << "\n";

  if (tests_failed == 0) {
    cout << "Все тесты успешно пройдены!\n";
  } else {
    cout << "Были обнаружены ошибки.\n";
  }
}
