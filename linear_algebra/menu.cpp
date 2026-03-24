#include <cmath>
#include <iostream>

#include "LU_decompozition.h"
#include "gauss_method.h"
#include "matrixes.h"
#include "menu.h"

const int MAX_MATRICES = 20;

using namespace std;

static matrixes<double> *matrices[MAX_MATRICES];
static int matrices_count = 0;

void PrintInfo() {
    cout << "1 - Создать матрицу" << endl;
    cout << "2 - Заполнить матрицей гильберта" << endl;
    cout << "3 - Заполнить единичной матрицей" << endl;
    cout << "4 - Решить СЛАУ методом Гаусса с выбором ведущего элемента  "
         << endl;
    cout << "5 - Решить СЛАУ LU разложением " << endl;
    cout << "6 - Сетнуть элемент матрицы" << endl;
    cout << "7 - Вывести матрицу" << endl;
    cout << "8 - Выйти" << endl;
}

void Print_Matrix(const matrixes<double> *matrix) {
    if (!matrix)
        return;
    cout << "| ";
    for (int index = 0; index < matrix->GetRows() * matrix->GetCols();
         index++) {
        cout << matrix->Get(index) << " ";
        if ((index + 1) % matrix->GetCols() == 0) {
            cout << "|" << endl;
            if (index + 1 < matrix->GetRows() * matrix->GetCols()) {
                cout << "| ";
            }
        }
    }
}

void open_menu() {
    int choice = -1;
    while (choice != 8) {
        PrintInfo();
        cout << "Введите выбор: ";
        if (!(cin >> choice)) {
            cin.clear();
            cin.ignore(10000, '\n');
            continue;
        }

        switch (choice) {
        case 1: {
            if (matrices_count >= MAX_MATRICES) {
                cout << "Достигнуто максимальное количество матриц" << endl;
                break;
            }
            int r, c;
            cout << "Введите количество строк и столбцов через пробел: ";
            cin >> r >> c;
            if (r > 0 && c > 0) {
                matrices[matrices_count] = new matrixes<double>(r, c);
                cout << "Матрица создана под индексом " << matrices_count
                     << endl;
                matrices_count++;
            } else {
                cout << "Неверные размеры" << endl;
            }
            break;
        }
        case 2: {
            int index;
            cout << "Введите индекс матрицы (от 0 до " << matrices_count - 1
                 << "): ";
            cin >> index;
            if (index >= 0 && index < matrices_count && matrices[index]) {
                matrices[index]->MakeGilbert();
                cout << "Матрица заполнена матрицей Гильберта" << endl;
            } else {
                cout << "Неверный индекс" << endl;
            }
            break;
        }
        case 3: {
            int index;
            cout << "Введите индекс матрицы: ";
            cin >> index;
            if (index >= 0 && index < matrices_count && matrices[index]) {
                matrices[index]->MakeOnes();
                cout << "Матрица заполнена единичной матрицей" << endl;
            } else {
                cout << "Неверный индекс" << endl;
            }
            break;
        }
        case 4: {
            int index;
            cout << "Введите индекс матрицы: ";
            cin >> index;
            if (index >= 0 && index < matrices_count && matrices[index]) {
                if (matrices[index]->GetRows() != matrices[index]->GetCols()) {
                    cout << "ОШИБКА: Матрица должна быть квадратной" << endl;
                    break;
                }
                int n = matrices[index]->GetRows();

                cout << "Введите " << n
                     << " элементов вектора свободных членов b: ";
                double *b = new double[n];
                for (int i = 0; i < n; i++) {
                    cin >> b[i];
                }

                gauss_method<double> gm(matrices[index]->GetMatrix(), n);

                double *x = gm.Solve(b);
                if (x) {
                    cout << "Решение СЛАУ методом Гаусса (вектор x): " << endl;
                    for (int i = 0; i < n; i++) {
                        cout << x[i] << " ";
                    }
                    cout << endl;
                    delete[] x;
                } else {
                    cout << "ОШИБКА: Не удалось получить решение СЛАУ." << endl;
                }
                delete[] b;
            } else {
                cout << "Неверный индекс" << endl;
            }
            break;
        }
        case 5: {
            int index;
            cout << "Введите индекс матрицы: ";
            cin >> index;
            if (index >= 0 && index < matrices_count && matrices[index]) {
                if (matrices[index]->GetRows() != matrices[index]->GetCols()) {
                    cout << "ОШИБКА: Матрица должна быть квадратной" << endl;
                    break;
                }
                int n = matrices[index]->GetRows();

                cout << "Введите " << n
                     << " элементов вектора свободных членов b: ";
                double *b = new double[n];
                for (int i = 0; i < n; i++) {
                    cin >> b[i];
                }

                LU_Decomposition<double> lu(matrices[index]->GetMatrix(), n);
                double det = lu.GetDet();
                if (std::abs(det) < 1e-9) {
                    cout << "ОШИБКА: Матрица вырождена (определитель близок к "
                            "0)."
                         << endl;
                    delete[] b;
                    break;
                }

                double *x = lu.Solve(b);

                if (x) {
                    cout << "Решение СЛАУ (вектор x): " << endl;
                    for (int i = 0; i < n; i++) {
                        cout << x[i] << " ";
                    }
                    cout << endl;
                    delete[] x; // предполагается, что метод возвращает
                                // динамический массив
                } else {
                    cout << "ОШИБКА: Не удалось получить решение СЛАУ." << endl;
                }

                delete[] b;

            } else {
                cout << "Неверный индекс" << endl;
            }
            break;
        }
        case 6: {
            int index;
            cout << "Введите индекс матрицы: ";
            cin >> index;
            if (index >= 0 && index < matrices_count && matrices[index]) {
                int r, c;
                double val;
                cout << "Введите номер строки (с 0), столбца (с 0) и новое "
                        "значение: ";
                cin >> r >> c >> val;
                if (r >= 0 && r < matrices[index]->GetRows() && c >= 0 &&
                    c < matrices[index]->GetCols()) {
                    matrices[index]
                        ->GetMatrix()[r * matrices[index]->GetCols() + c] = val;
                    cout << "Элемент изменен" << endl;
                } else {
                    cout << "ОШИБКА: Неверные координаты" << endl;
                }
            } else {
                cout << "Неверный индекс" << endl;
            }
            break;
        }
        case 7: {
            int index;
            cout << "Введите индекс матрицы: ";
            cin >> index;
            if (index >= 0 && index < matrices_count && matrices[index]) {
                cout << "Матрица [" << index << "]:" << endl;
                Print_Matrix(matrices[index]);
            } else {
                cout << "Неверный индекс" << endl;
            }
            break;
        }
        case 8: {
            cout << "Завершение работы программы" << endl;
            for (int i = 0; i < matrices_count; i++) {
                delete matrices[i];
                matrices[i] = nullptr;
            }
            matrices_count = 0;
            break;
        }
        default:
            cout << "Неверный выбор" << endl;
            break;
        }
    }
}