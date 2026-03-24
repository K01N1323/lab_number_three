#include "gauss_method.h"
#include <cmath>
#include <vector>

// Конструктор метода Гаусса с выбором ведущего элемента
template <class T>
gauss_method<T>::gauss_method(T *m, int n) {
    matrix = new T[n * n];
    conmatrix = new T[n * n];
    this->n = n;
    this->det = T(1.0);
    this->swaps = 0;

    for (int i = 0; i < (this->n) * (this->n); i++)
        this->matrix[i] = m[i];

    for (int i = 0; i < (this->n) * (this->n); i++) {
        conmatrix[i] = T(0.0);
    }

    for (int i = 0; i < (this->n); i++) {
        conmatrix[n * i + i] = T(1.0);
    }
}

// Деструктор
template <class T>
gauss_method<T>::~gauss_method() {
    if (matrix) {
        delete[] matrix;
        matrix = nullptr;
    }

    delete[] conmatrix;
    conmatrix = nullptr;
}

// Возвращает знак определителя в зависимости от перестановок
template <class T>
T gauss_method<T>::GetDet() const {
    if (swaps % 2 == 0)
        return det;
    else
        return T(-1.0) * det;
}

// Ищет максимальный по модулю элемент в столбце для устойчивости
template <class T>
int gauss_method<T>::neededrow(int skip) const { // по дефолту скип равен нулю
    T loc, mx = std::abs(matrix[skip * n + skip]);
    int numrow = skip;

    for (int i = skip; i < n; i++) {
        loc = std::abs(matrix[n * i + skip]);
        if (loc > mx) {
            mx = loc;
            numrow = i;
        }
    }

    return numrow;
}

// Переставляет строки так, чтобы ведущий элемент был сверху
template <class T>
void gauss_method<T>::swaprows(int we_swap1) {
    int we_swap = neededrow(we_swap1);
    std::vector<T> loc(n);
    std::vector<T> loccon(n);

    swaps++;

    if (we_swap1 == we_swap)
        return;

    for (int i = 0; i < n; i++) {
        loc[i] = matrix[we_swap1 * n + i];
        loccon[i] = conmatrix[we_swap1 * n + i];
    }

    for (int i = 0; i < n; i++) {
        matrix[we_swap1 * n + i] = matrix[we_swap * n + i];
        matrix[we_swap * n + i] = loc[i];
        conmatrix[we_swap1 * n + i] = conmatrix[we_swap * n + i];
        conmatrix[we_swap * n + i] = loccon[i];
    }
}

// Делит текущую строку на ведущий элемент, делая его равным 1
template <class T>
void gauss_method<T>::divisionrow(int num) {
    T el = matrix[n * num + num];

    if (el == T(0)) {
        det = T(0);
        return;
    }

    det *= el;

    for (int i = 0; i < n; i++) {
        matrix[n * num + i] /= el;
        conmatrix[n * num + i] /= el;
    }
}

// Вычитает текущую строку из всех последующих, обнуляя элементы под ведущим
template <class T>
void gauss_method<T>::subtraction(int current) {
    std::vector<T> loc(n);
    std::vector<T> loccon(n);

    for (int i = 0; i < n; i++) {
        loc[i] = matrix[current * n + i];
        loccon[i] = conmatrix[current * n + i];
    }

    for (int k = (current + 1); k < n; k++) {
        T mnozh = matrix[k * n + current];
        for (int i = 0; i < n; i++) {
            matrix[k * n + i] -= mnozh * loc[i];
            conmatrix[k * n + i] -= mnozh * loccon[i];
        }
    }
}

// Прямой ход Гаусса: приведение матрицы к верхнетреугольному виду
template <class T>
void gauss_method<T>::triangle() {
    int skip = -1;
    for (int i = 0; i < n; i++) {
        skip++;
        swaprows(skip);
        divisionrow(skip);
        subtraction(skip);
    }
}

// Обратный ход Гаусса: нахождение корней (или обратной матрицы)
template <class T>
void gauss_method<T>::obrat() {
    T koof;

    for (int i = (n - 1); i >= 0; i--) {
        for (int k = (i - 1); k >= 0; k--) {
            koof = matrix[n * k + i];
            for (int j = 0; j < n; j++) {
                conmatrix[n * k + j] -= koof * conmatrix[n * i + j];
                matrix[n * k + j] -= koof * matrix[n * i + j];
            }
        }
    }
}

// Полный процесс метода Гаусса (прямой + обратный ход)
template <class T>
void gauss_method<T>::reverse() {
    triangle();
    obrat();
}

// Решает СЛАУ для вектора b и возвращает вектор x
template <class T>
T *gauss_method<T>::Solve(T *b) {
    for (int i = 0; i < n; i++) {
        conmatrix[n * i] = b[i];
    }

    triangle();
    obrat();

    T *result = new T[n];

    for (int i = 0; i < n; i++) {
        result[i] = conmatrix[n * i];
    }

    return result;
}

// Решение СЛАУ без выделения памяти под ответ (для ускорения тестов)
template <class T>
void gauss_method<T>::SolveForTests(T *b) {
    for (int i = 0; i < n; i++) {
        conmatrix[n * i] = b[i];
    }

    triangle();
    obrat();
}

// Инкапсулирует вызов reverse()
template <class T>
void gauss_method<T>::TakeReverse() { reverse(); }

// Геттеры для матрицы и её элементов
template <class T>
T *gauss_method<T>::GetMatrix() const { return this->conmatrix; }

template <class T>
T gauss_method<T>::GetElement(int i, int j) const {
    return matrix[n * i + j];
}