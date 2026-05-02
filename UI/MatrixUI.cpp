#include "MatrixUI.h"
#include "../LinearAlgebra/DiagonalMatrix.h"
#include "../LinearAlgebra/GaussMethod.h"
#include "../LinearAlgebra/LUDecompozition.h"
#include "../LinearAlgebra/LowerTriangleMatrix.h"
#include "../LinearAlgebra/RectangularMatrix.h"
#include "../LinearAlgebra/SquareMatrix.h"
#include "../LinearAlgebra/UpperTriangleMatrix.h"
#include <iostream>
#include <limits>
#include <string>

using namespace std;

// Безопасное чтение целого числа с повторным запросом при ошибке ввода
static int ReadInt(const char *prompt) {
  int value;
  while (true) {
    cout << prompt;
    cin >> value;
    if (cin.fail()) {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Ошибка: введите целое число\n";
    } else {
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      return value;
    }
  }
}

// Безопасное чтение вещественного числа с повторным запросом при ошибке ввода
static double ReadDouble(const char *prompt) {
  double value;
  while (true) {
    cout << prompt;
    cin >> value;
    if (cin.fail()) {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Ошибка: введите число\n";
    } else {
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      return value;
    }
  }
}

#define MaxMatrices 50
static Matrix<double> *matrices[MaxMatrices];
static int MatricesCount = 0;

void PrintMenu(void) {
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
  cout << "20 – выполнить Map (возвести элементы в квадрат)\n";
  cout << "21 – выполнить Reduce (сумма элементов)\n";
  cout << "-1 – выйти из программы\n";
}

int actions(int flag) {
  int rows, cols, num, i, j;
  int num1, num2;
  double val;

  switch (flag) {
  case -1: {
    for (int k = 0; k < MatricesCount; k++) {
      if (matrices[k] != nullptr) {
        delete matrices[k];
      }
    }
    return 1;
  }

  case 1: {
    rows = ReadInt("Введите размер матрицы (N): ");
    if (rows <= 0) {
      cout << "Размер матрицы должен быть положительным\n";
      break;
    }
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot != -1) {
        matrices[slot] = new SquareMatrix<double>(rows);
        cout << "Матрица создана под номером " << slot << "\n";
      } else {
        cout << "Нет места для новой матрицы\n";
      }
    }
    break;
  }

  case 2: {
    rows = ReadInt("Введите количество строк: ");
    cols = ReadInt("Введите количество столбцов: ");
    if (rows <= 0 || cols <= 0) {
      cout << "Размеры матрицы должны быть положительными\n";
      break;
    }
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot != -1) {
        matrices[slot] = new RectangularMatrix<double>(rows, cols);
        cout << "Матрица создана под номером " << slot << "\n";
      } else {
        cout << "Нет места для новой матрицы\n";
      }
    }
    break;
  }

  case 3: {
    rows = ReadInt("Введите размер матрицы (N): ");
    if (rows <= 0) {
      cout << "Размер матрицы должен быть положительным\n";
      break;
    }
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot != -1) {
        matrices[slot] = new DiagonalMatrix<double>(rows);
        cout << "Матрица создана под номером " << slot << "\n";
      } else {
        cout << "Нет места для новой матрицы\n";
      }
    }
    break;
  }

  case 4: {
    rows = ReadInt("Введите размер матрицы (N): ");
    if (rows <= 0) {
      cout << "Размер матрицы должен быть положительным\n";
      break;
    }
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot != -1) {
        matrices[slot] = new UpperTriangularMatrix<double>(rows);
        cout << "Матрица создана под номером " << slot << "\n";
      } else {
        cout << "Нет места для новой матрицы\n";
      }
    }
    break;
  }

  case 5: {
    rows = ReadInt("Введите размер матрицы (N): ");
    if (rows <= 0) {
      cout << "Размер матрицы должен быть положительным\n";
      break;
    }
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot != -1) {
        matrices[slot] = new LowerTriangularMatrix<double>(rows);
        cout << "Матрица создана под номером " << slot << "\n";
      } else {
        cout << "Нет места для новой матрицы\n";
      }
    }
    break;
  }

  case 6: {
    num =
        ReadInt("Введите номер матрицы, которую вы хотите вывести на экран: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
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
    num = ReadInt("Введите номер матрицы: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером нет\n";
      break;
    }
    i = ReadInt("Введите номер строки: ");
    j = ReadInt("Введите номер столбца: ");
    {
      Matrix<double> *m = matrices[num];
      if (i >= m->GetRows() || j >= m->GetCols() || i < 0 || j < 0) {
        cout << "Элементов под таким номером в матрице нет\n";
      } else {
        val = ReadDouble("Введите значение (double): ");
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
    num = ReadInt("Введите номер матрицы: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером нет\n";
      break;
    }
    i = ReadInt("Введите номер строки: ");
    j = ReadInt("Введите номер столбца: ");
    {
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
    num1 = ReadInt("Введите номер первой матрицы: ");
    num2 = ReadInt("Введите номер второй матрицы: ");
    if (num1 < 0 || num2 < 0 || num1 >= MatricesCount ||
        num2 >= MatricesCount || matrices[num1] == nullptr ||
        matrices[num2] == nullptr) {
      cout << "Матриц с таким номером не существует\n";
    } else {
      Matrix<double> *a = matrices[num1];
      Matrix<double> *b = matrices[num2];
      if (a->GetRows() != b->GetRows() || a->GetCols() != b->GetCols()) {
        cout << "Размеры матриц не совпадают\n";
        return 0;
      }
      {
        int slot = -1;
        for (int k = 0; k < MatricesCount; k++) {
          if (matrices[k] == nullptr) { slot = k; break; }
        }
        if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
        if (slot == -1) {
          cout << "Нет места для новой матрицы\n";
          return 0;
        }
        matrices[slot] = *a + *b;
        cout << "Результат сложения записан в матрицу под номером "
             << slot << "\n";
      }
    }
    break;
  }

  case 10: {
    num1 = ReadInt("Введите номер первой матрицы: ");
    num2 = ReadInt("Введите номер второй матрицы (вычитаемая): ");
    if (num1 < 0 || num2 < 0 || num1 >= MatricesCount ||
        num2 >= MatricesCount || matrices[num1] == nullptr ||
        matrices[num2] == nullptr) {
      cout << "Матриц с таким номером не существует\n";
    } else {
      Matrix<double> *a = matrices[num1];
      Matrix<double> *b = matrices[num2];
      if (a->GetRows() != b->GetRows() || a->GetCols() != b->GetCols()) {
        cout << "Размеры матриц не совпадают\n";
        return 0;
      }
      {
        int slot = -1;
        for (int k = 0; k < MatricesCount; k++) {
          if (matrices[k] == nullptr) { slot = k; break; }
        }
        if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
        if (slot == -1) {
          cout << "Нет места для новой матрицы\n";
          return 0;
        }
        matrices[slot] = *a - *b;
        cout << "Результат вычитания записан в матрицу под номером "
             << slot << "\n";
      }
    }
    break;
  }

  case 11: {
    num1 = ReadInt("Введите номер первой матрицы: ");
    num2 = ReadInt("Введите номер второй матрицы: ");
    if (num1 < 0 || num2 < 0 || num1 >= MatricesCount ||
        num2 >= MatricesCount || matrices[num1] == nullptr ||
        matrices[num2] == nullptr) {
      cout << "Матриц с таким номером не существует\n";
      return 0;
    }
    {
      Matrix<double> *a = matrices[num1];
      Matrix<double> *b = matrices[num2];
      if (a->GetCols() != b->GetRows()) {
        cout << "Нельзя перемножить матрицы с такими размерами\n";
        return 0;
      }
      {
        int slot = -1;
        for (int k = 0; k < MatricesCount; k++) {
          if (matrices[k] == nullptr) { slot = k; break; }
        }
        if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
        if (slot != -1) {
          matrices[slot] = *a * *b;
          cout << "Результат умножения записан в матрицу под номером "
               << slot << "\n";
        } else {
          cout << "Нет места для новой матрицы\n";
        }
      }
    }
    break;
  }

  case 12: {
    num = ReadInt(
        "Введите номер матрицы, которую вы хотите умножить на скаляр: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    val = ReadDouble("Введите скаляр: ");
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot == -1) {
        cout << "Нет места для новой матрицы\n";
        return 0;
      }
      matrices[slot] = *matrices[num] * val;
      cout << "Результат умножения на скаляр записан в матрицу под номером "
           << slot << "\n";
    }
    break;
  }

  case 13: {
    num = ReadInt("Введите номер матрицы, для которой нужно найти обратную: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    if (matrices[num]->GetRows() != matrices[num]->GetCols()) {
      cout << "Матрица должна быть квадратной\n";
      return 0;
    }
    try {
      Matrix<double> *inv = matrices[num]->GetInverseMatrix();
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot == -1) {
        delete inv;
        cout << "Нет места для новой матрицы\n";
        return 0;
      }
      matrices[slot] = inv;
      cout << "Обратная матрица записана под номером " << slot << "\n";
    } catch (const exception &e) {
      cout << "Ошибка: " << e.what() << "\n";
    }
    break;
  }

  case 14: {
    num = ReadInt(
        "Введите номер матрицы, определитель которой нужно вычислить: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
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
    num = ReadInt("Введите номер матрицы для вычисления нормы Фробениуса: ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    cout << "Норма = " << matrices[num]->GetFrobeniusNorm() << "\n";
    break;
  }

  case 16:
  case 17: {
    num = ReadInt("Введите номер матрицы (квадратной): ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    {
      Matrix<double> *m = matrices[num];
      if (m->GetRows() != m->GetCols()) {
        cout << "Матрица должна быть квадратной\n";
        return 0;
      }
      int n = m->GetRows();
      MutableArraySequence<double> b;
      cout << "Введите " << n << " элементов вектора:\n";
      for (int k = 0; k < n; k++) {
        string prompt = "  b[" + to_string(k) + "] = ";
        val = ReadDouble(prompt.c_str());
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
    }
    break;
  }

  case 18: {
    if (MatricesCount == 0)
      cout << "Нет доступных для взаимодействия матриц\n";
    else {
      cout << "Доступные для взаимодействия матрицы: ";
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] != nullptr)
          cout << k << " ";
      }
      cout << "\n";
    }
    break;
  }

  case 19: {
    num = ReadInt("Введите номер матрицы, которую вы хотите удалить: ");
    if (num >= MatricesCount || num < 0 || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    delete matrices[num];
    matrices[num] = nullptr;
    cout << "Матрица успешно удалена\n";
    break;
  }

  case 20: {
    num = ReadInt("Введите номер матрицы для Map (x^2): ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    {
      int slot = -1;
      for (int k = 0; k < MatricesCount; k++) {
        if (matrices[k] == nullptr) { slot = k; break; }
      }
      if (slot == -1 && MatricesCount < MaxMatrices) { slot = MatricesCount++; }
      if (slot == -1) {
        cout << "Нет места для новой матрицы\n";
        return 0;
      }
      matrices[slot] = matrices[num]->map([](const double &x) -> double { return x * x; });
      cout << "Результат Map записан в матрицу под номером " << slot << "\n";
    }
    break;
  }

  case 21: {
    num = ReadInt("Введите номер матрицы для Reduce (сумма): ");
    if (num < 0 || num >= MatricesCount || matrices[num] == nullptr) {
      cout << "Матрицы с таким номером не существует\n";
      return 0;
    }
    {
      double result = matrices[num]->reduce(
          [](const double &a, const double &b) -> double { return a + b; },
          0.0);
      cout << "Сумма всех элементов матрицы: " << result << "\n";
    }
    break;
  }

  default:
    cout << "Команда, которую вы ввели не существует\n";
    break;
  }
  return 0;
}

void RunUI(void) {
  int flag = 0;
  while (true) {
    PrintMenu();
    cout << ">>> ";
    cin >> flag;
    if (cin.eof()) {
      break;
    }
    if (cin.fail()) {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Ошибка: введите номер команды (целое число)\n";
      continue;
    }
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    if (flag != -1 && (flag < 1 || flag > 21)) {
      cout << "Неизвестная команда. Введите число от 1 до 21 или -1 для выхода\n";
      continue;
    }
    if (actions(flag) == 1) {
      break;
    }
  }
}
