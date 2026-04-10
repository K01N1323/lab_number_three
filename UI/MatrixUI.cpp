#include "MatrixUI.h"
#include "../LinearAlgebra/DiagonalMatrix.h"
#include "../LinearAlgebra/GaussMethod.h"
#include "../LinearAlgebra/LUDecompozition.h"
#include "../LinearAlgebra/LowerTriangleMatrix.h"
#include "../LinearAlgebra/RectangularMatrix.h"
#include "../LinearAlgebra/SquareMatrix.h"
#include "../LinearAlgebra/UpperTriangleMatrix.h"
#include <iostream>

using namespace std;

#define max_matrices 50
static Matrix<double> *matrices[max_matrices];
static int matrices_count = 0;

void print_menu(void) {
  cout << "Выберите действие:\n";
  cout << "1 – создать квадратную матрицу\n";
  cout << "2 – создать прямоугольную матрицу\n";
  cout << "3 – создать диагональную матрицу\n";
  cout << "4 – создать верхнетреугольную матрицу\n";
  cout << "5 – создать нижнетреугольную матрицу\n";
  cout << "6 – вывести матрицу\n";
  cout << "7 – сетнуть элемент матрицы\n";
  cout << "8 – вывести элемент матрицы\n";
  cout << "9 – сложить две матрицы\n";
  cout << "10 – вычесть две матрицы\n";
  cout << "11 – умножить две матрицы\n";
  cout << "12 – умножить матрицу на скаляр\n";
  cout << "13 – найти обратную матрицу\n";
  cout << "14 – вычилить определитель\n";
  cout << "15 – вычислить норму\n";
  cout << "16 – решить СЛАУ (Гаусс)\n";
  cout << "17 – решить СЛАУ (LU-разложение)\n";
  cout << "18 - вывести номера матриц, доступных для взаимодействия\n";
  cout << "19 – очистить матрицу\n";
  cout << "-1 – выйти из программы\n";
}

int actions(int flag) {
  int rows, cols, num, i, j;
  int num1, num2;
  double val;

  switch (flag) {
  case -1: {
    for (int k = 0; k < matrices_count; k++) {
      if (matrices[k] != nullptr) {
        delete matrices[k];
      }
    }
    return 1;
  }

  case 1: {
    cout << "Введите размер матрицы (N): ";
    cin >> rows;
    if (matrices_count < max_matrices) {
      matrices[matrices_count] = new SquareMatrix<double>(rows);
      cout << "Матрица создана под номером " << matrices_count << "\n";
      matrices_count++;
    } else {
      cout << "Нет места для новой матрицы\n";
    }
    break;
  }

  case 2: {
    cout << "Введите количество строк и столбцов в матрице: ";
    cin >> rows >> cols;
    if (matrices_count < max_matrices) {
      matrices[matrices_count] = new RectangularMatrix<double>(rows, cols);
      cout << "Матрица создана под номером " << matrices_count << "\n";
      matrices_count++;
    } else {
      cout << "Нет места для новой матрицы\n";
    }
    break;
  }

  case 3: {
    cout << "Введите размер матрицы (N): ";
    cin >> rows;
    if (matrices_count < max_matrices) {
      matrices[matrices_count] = new DiagonalMatrix<double>(rows);
      cout << "Матрица создана под номером " << matrices_count << "\n";
      matrices_count++;
    } else {
      cout << "Нет места для новой матрицы\n";
    }
    break;
  }

  case 4: {
    cout << "Введите размер матрицы (N): ";
    cin >> rows;
    if (matrices_count < max_matrices) {
      matrices[matrices_count] = new UpperTriangularMatrix<double>(rows);
      cout << "Матрица создана под номером " << matrices_count << "\n";
      matrices_count++;
    } else {
      cout << "Нет места для новой матрицы\n";
    }
    break;
  }

  case 5: {
    cout << "Введите размер матрицы (N): ";
    cin >> rows;
    if (matrices_count < max_matrices) {
      matrices[matrices_count] = new LowerTriangularMatrix<double>(rows);
      cout << "Матрица создана под номером " << matrices_count << "\n";
      matrices_count++;
    } else {
      cout << "Нет места для новой матрицы\n";
    }
    break;
  }

  case 6: {
    cout << "Введите номер матрицы, которую вы хотите вывести на экран: ";
    cin >> num;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером нет\n";
    } else {
      Matrix<double> *m = matrices[num];
      for (int r = 0; r < m->GetRows(); r++) {
        for (int c = 0; c < m->GetCols(); c++) {
          cout << m->GetIJ(r, c) << " ";
        }
        cout << "\n";
      }
    }
    break;
  }

  case 7: {
    cout << "Введите номер матрицы и номер строки и столбца элемента, который "
            "вы хотите сетнуть: ";
    cin >> num >> i >> j;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером нет\n";
    } else {
      Matrix<double> *m = matrices[num];
      if (i >= m->GetRows() || j >= m->GetCols() || i < 0 || j < 0) {
        cout << "Элементов под таким номером в матрице нет\n";
      } else {
        cout << "Введите значение (double): ";
        cin >> val;
        try {
          m->SetIJ(i, j, val);
        } catch (const exception &e) {
          cout << "Ошибка записи: " << e.what() << "\n";
        }
      }
    }
    break;
  }

  case 8: {
    cout << "Введите номер матрицы и номер строки и столбца элемента, который "
            "вы хотите вывести на экран: ";
    cin >> num >> i >> j;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером нет\n";
    } else {
      Matrix<double> *m = matrices[num];
      if (i >= m->GetRows() || j >= m->GetCols() || i < 0 || j < 0) {
        cout << "Элементов под таким номером в матрице нет\n";
      } else {
        cout << "Элемент [" << i << "][" << j << "]: " << m->GetIJ(i, j)
             << "\n";
      }
    }
    break;
  }

  case 9: {
    cout << "Введите номера матриц, которые вы хотите сложить: ";
    cin >> num1 >> num2;
    if (num1 < 0 || num2 < 0 || num1 >= matrices_count ||
        num2 >= matrices_count || matrices[num1] == nullptr ||
        matrices[num2] == nullptr) {
      cout << "Матриц с таким номером не существует\n";
    } else {
      Matrix<double> *a = matrices[num1];
      Matrix<double> *b = matrices[num2];
      if (a->GetRows() != b->GetRows() || a->GetCols() != b->GetCols()) {
        cout << "Размеры матриц не совпадают\n";
        return 0;
      }
      if (matrices_count >= max_matrices) {
        cout << "Нет места для новой матрицы\n";
        return 0;
      }
      matrices[matrices_count] = *a + *b;
      cout << "Результат сложения записан в матрицу под номером "
           << matrices_count << "\n";
      matrices_count++;
    }
    break;
  }

  case 10: {
    cout << "Введите номера матриц, которые вы хотите вычесть (первая минус "
            "вторая): ";
    cin >> num1 >> num2;
    if (num1 < 0 || num2 < 0 || num1 >= matrices_count ||
        num2 >= matrices_count || matrices[num1] == nullptr ||
        matrices[num2] == nullptr) {
      cout << "Матриц с таким номером не существует\n";
    } else {
      Matrix<double> *a = matrices[num1];
      Matrix<double> *b = matrices[num2];
      if (a->GetRows() != b->GetRows() || a->GetCols() != b->GetCols()) {
        cout << "Размеры матриц не совпадают\n";
        return 0;
      }
      if (matrices_count >= max_matrices) {
        cout << "Нет места для новой матрицы\n";
        return 0;
      }
      matrices[matrices_count] = *a - *b;
      cout << "Результат вычитания записан в матрицу под номером "
           << matrices_count << "\n";
      matrices_count++;
    }
    break;
  }

  case 11: {
    cout << "Введите номера матриц, которые вы хотите перемножить: ";
    cin >> num1 >> num2;
    if (num1 < 0 || num2 < 0 || num1 >= matrices_count ||
        num2 >= matrices_count || matrices[num1] == nullptr ||
        matrices[num2] == nullptr) {
      cout << "Матриц с таким номером не существует\n";
      return 0;
    }
    Matrix<double> *a = matrices[num1];
    Matrix<double> *b = matrices[num2];
    if (a->GetCols() != b->GetRows()) {
      cout << "Нельзя перемножить матрицы с такими размерами\n";
      return 0;
    }
    if (matrices_count >= max_matrices) {
      cout << "Нет места для новой матрицы\n";
      return 0;
    }
    matrices[matrices_count] = *a * *b;
    cout << "Результат умножения записан в матрицу под номером "
         << matrices_count << "\n";
    matrices_count++;
    break;
  }

  case 12: {
    cout << "Введите номер матрицы, которую вы хотите умножить на скаляр: ";
    cin >> num;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    if (matrices_count >= max_matrices) {
      cout << "Нет места для новой матрицы\n";
      return 0;
    }
    cout << "Введите скаляр: ";
    cin >> val;
    matrices[matrices_count] = *matrices[num] * val;
    cout << "Результат умножения на скаляр записан в матрицу под номером "
         << matrices_count << "\n";
    matrices_count++;
    break;
  }

  case 13: {
    cout << "Введите номер матрицы, для которой нужно найти обратную: ";
    cin >> num;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    if (matrices[num]->GetRows() != matrices[num]->GetCols()) {
      cout << "Матрица должна быть квадратной\n";
      return 0;
    }
    if (matrices_count >= max_matrices) {
      cout << "Нет места для новой матрицы\n";
      return 0;
    }
    try {
      matrices[matrices_count] = matrices[num]->GetInverseMatrix();
      cout << "Обратная матрица записана под номером " << matrices_count
           << "\n";
      matrices_count++;
    } catch (const exception &e) {
      cout << "Ошибка: " << e.what() << "\n";
    }
    break;
  }

  case 14: {
    cout << "Введите номер матрицы, определитель которой нужно вычислить: ";
    cin >> num;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    if (matrices[num]->GetRows() != matrices[num]->GetCols()) {
      cout << "Матрица должна быть квадратной\n";
      return 0;
    }
    cout << "Определитель = " << matrices[num]->GetDet() << "\n";
    break;
  }

  case 15: {
    cout << "Введите номер матрицы для вычисления нормы Фробениуса: ";
    cin >> num;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    cout << "Норма = " << matrices[num]->GetFrobeniusNorm() << "\n";
    break;
  }

  case 16:
  case 17: {
    cout << "Введите номер матрицы (квадратной): ";
    cin >> num;
    if (num < 0 || num >= matrices_count || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    Matrix<double> *m = matrices[num];
    if (m->GetRows() != m->GetCols()) {
      cout << "Матрица должна быть квадратной\n";
      return 0;
    }
    int n = m->GetRows();
    MutableArraySequence<double> b;
    cout << "Введите " << n << " элементов вектора (через пробел):\n";
    for (int k = 0; k < n; k++) {
      cin >> val;
      b.Append(val);
    }
    try {
      Sequence<double> *x = nullptr;
      if (flag == 16)
        x = m->SolveSlauGauss(&b);
      else
        x = m->SolveSlauLU(&b);

      cout << "Вектор решений (X):\n";
      for (int k = 0; k < x->GetLength(); k++) {
        cout << "X[" << k << "] = " << x->Get(k) << "\n";
      }
      delete x;
    } catch (const exception &e) {
      cout << "Ошибка решения: " << e.what() << "\n";
    }
    break;
  }

  case 18: {
    if (matrices_count == 0)
      cout << "Нет доступных для взаимодействия матриц\n";
    else {
      cout << "Доступные для взаимодействия матрицы: ";
      for (int k = 0; k < matrices_count; k++) {
        if (matrices[k] != nullptr)
          cout << k << " ";
      }
      cout << "\n";
    }
    break;
  }

  case 19: {
    cout << "Введите номер матрицы, которую вы хотите удалить: ";
    cin >> num;
    if (num >= matrices_count || num < 0 || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    delete matrices[num];
    matrices[num] = nullptr;
    cout << "Матрица успешно удалена\n";
    break;
  }

  default:
    cout << "Команда, которую вы ввели не существует\n";
    break;
  }
  return 0;
}

void run_ui(void) {
  int flag = 0;
  while (true) {
    print_menu();
    cin >> flag;
    if (actions(flag) == 1) {
      break;
    }
  }
}
