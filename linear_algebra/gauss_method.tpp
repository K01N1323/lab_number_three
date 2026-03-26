#include "MutableArraySequence.h" // Для возврата вектора ответов
#include "gauss_method.h"
#include <cmath>

template <class T> T clone_mapper(const T &val) { return val; }
template <class T> T zero_mapper(const T &val) { return T(0); }

template <class T> gauss_method<T>::gauss_method(const matrix<T> *input) {
  this->n = input->GetRows();
  this->det = T(1.0);
  this->swaps = 0;

  this->matrix_ptr = input->Map(clone_mapper<T>);

  this->conmatrix_ptr = input->Map(zero_mapper<T>);
  for (int i = 0; i < n; i++) {
    this->conmatrix_ptr->SetIJ(i, i, T(1.0));
  }
}

template <class T> gauss_method<T>::~gauss_method() {
  delete matrix_ptr;
  delete conmatrix_ptr;
}

template <class T> T gauss_method<T>::GetDet() const {
  if (swaps % 2 == 0)
    return det;
  else
    return T(-1.0) * det;
}

template <class T> int gauss_method<T>::neededrow(int skip) const {
  T mx = std::abs(matrix_ptr->GetIJ(skip, skip));
  int numrow = skip;

  for (int i = skip + 1; i < n; i++) {
    T loc = std::abs(matrix_ptr->GetIJ(i, skip));
    if (loc > mx) {
      mx = loc;
      numrow = i;
    }
  }
  return numrow;
}

template <class T> void gauss_method<T>::swaprows(int we_swap1) {
  int we_swap = neededrow(we_swap1);

  if (we_swap1 != we_swap) {
    swaps++;
    for (int i = 0; i < n; i++) {

      T loc = matrix_ptr->GetIJ(we_swap1, i);
      T loccon = conmatrix_ptr->GetIJ(we_swap1, i);

      matrix_ptr->SetIJ(we_swap1, i, matrix_ptr->GetIJ(we_swap, i));
      matrix_ptr->SetIJ(we_swap, i, loc);

      conmatrix_ptr->SetIJ(we_swap1, i, conmatrix_ptr->GetIJ(we_swap, i));
      conmatrix_ptr->SetIJ(we_swap, i, loccon);
    }
  }
}

template <class T> void gauss_method<T>::divisionrow(int num) {
  T el = matrix_ptr->GetIJ(num, num);

  if (el == T(0)) {
    det = T(0);
    return;
  }

  det *= el;

  for (int i = 0; i < n; i++) {
    matrix_ptr->SetIJ(num, i, matrix_ptr->GetIJ(num, i) / el);
    conmatrix_ptr->SetIJ(num, i, conmatrix_ptr->GetIJ(num, i) / el);
  }
}

template <class T> void gauss_method<T>::subtraction(int current) {
  for (int k = current + 1; k < n; k++) {
    T mnozh = matrix_ptr->GetIJ(k, current);
    for (int i = 0; i < n; i++) {
      matrix_ptr->SetIJ(k, i,
                        matrix_ptr->GetIJ(k, i) -
                            mnozh * matrix_ptr->GetIJ(current, i));
      conmatrix_ptr->SetIJ(k, i,
                           conmatrix_ptr->GetIJ(k, i) -
                               mnozh * conmatrix_ptr->GetIJ(current, i));
    }
  }
}

template <class T> void gauss_method<T>::triangle() {
  for (int i = 0; i < n; i++) {
    swaprows(i);
    divisionrow(i);
    subtraction(i);
  }
}

template <class T> void gauss_method<T>::obrat() {
  for (int i = n - 1; i >= 0; i--) {
    for (int k = i - 1; k >= 0; k--) {
      T koof = matrix_ptr->GetIJ(k, i);
      for (int j = 0; j < n; j++) {
        conmatrix_ptr->SetIJ(k, j,
                             conmatrix_ptr->GetIJ(k, j) -
                                 koof * conmatrix_ptr->GetIJ(i, j));
        matrix_ptr->SetIJ(
            k, j, matrix_ptr->GetIJ(k, j) - koof * matrix_ptr->GetIJ(i, j));
      }
    }
  }
}

template <class T> void gauss_method<T>::reverse() {
  triangle();
  obrat();
}

template <class T> Sequence<T> *gauss_method<T>::Solve(const Sequence<T> *b) {

  for (int i = 0; i < n; i++) {
    conmatrix_ptr->SetIJ(i, 0, b->Get(i));
  }

  triangle();
  obrat();

  Sequence<T> *result = new MutableArraySequence<T>();
  for (int i = 0; i < n; i++) {
    result->Append(conmatrix_ptr->GetIJ(i, 0));
  }

  return result;
}

template <class T> void gauss_method<T>::SolveForTests(const Sequence<T> *b) {
  for (int i = 0; i < n; i++) {
    conmatrix_ptr->SetIJ(i, 0, b->Get(i));
  }
  triangle();
  obrat();
}

template <class T> void gauss_method<T>::TakeReverse() { reverse(); }

template <class T> matrix<T> *gauss_method<T>::GetInverseMatrix() const {
  // Создаем копию матрицы ответов и отдаем её наружу
  return conmatrix_ptr->Map(clone_mapper<T>);
}

template <class T> T gauss_method<T>::GetElement(int i, int j) const {
  return matrix_ptr->GetIJ(i, j);
}