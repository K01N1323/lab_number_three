#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdexcept>

#include "IEnumerator.h"

// Абстрактный базовый класс для всех последовательностей
template <class T> class Sequence {
public:
  virtual ~Sequence() {}

  virtual const T &GetFirst() const = 0;
  virtual const T &GetLast() const = 0;
  virtual const T &Get(int index) const = 0;
  virtual void Set(int index, const T &value) = 0;
  virtual int GetLength() const = 0;
  virtual Sequence<T> *GetSubsequence(int startIndex, int endIndex) const = 0;
  virtual Sequence<T> *Append(const T &item) = 0;
  virtual Sequence<T> *Prepend(const T &item) = 0;
  virtual Sequence<T> *InsertAt(const T &item, int index) = 0;
  virtual Sequence<T> *Concat(Sequence<T> *list) = 0;

  const T &operator[](int index) const { return this->Get(index); }

  virtual Sequence<T> *Instance() = 0;
  virtual Sequence<T> *CreateEmpty() const = 0;
  virtual IEnumerator<T> *GetEnumerator() const = 0;

  // Применяет функцию mapper к каждому элементу, возвращает новую последовательность
  Sequence<T> *Map(T (*mapper)(const T &)) const {
    Sequence<T> *result = this->CreateEmpty();

    IEnumerator<T> *it = this->GetEnumerator();
    if (it != nullptr) {
      while (it->HasNext()) {
        Sequence<T> *old_ptr = result;
        result = result->Append(mapper(it->GetCurrent()));
        if (result != old_ptr)
          delete old_ptr;
        it->MoveNext();
      }
      delete it;
    } else {
      for (int index = 0; index < this->GetLength(); index++) {
        Sequence<T> *old_ptr = result;
        result = result->Append(mapper(this->Get(index)));
        if (result != old_ptr)
          delete old_ptr;
      }
    }

    return result;
  }

  // Фильтрует элементы по предикату, возвращает новую последовательность
  Sequence<T> *Where(bool (*where)(const T &)) const {
    Sequence<T> *result = this->CreateEmpty();

    IEnumerator<T> *it = this->GetEnumerator();
    if (it != nullptr) {
      while (it->HasNext()) {
        if (where(it->GetCurrent())) {
          Sequence<T> *old_ptr = result;
          result = result->Append(it->GetCurrent());
          if (result != old_ptr)
            delete old_ptr;
        }
        it->MoveNext();
      }
      delete it;
    } else {
      for (int index = 0; index < this->GetLength(); index++) {
        if (where(this->Get(index))) {
          Sequence<T> *old_ptr = result;
          result = result->Append(this->Get(index));
          if (result != old_ptr)
            delete old_ptr;
        }
      }
    }

    return result;
  }

  // Сворачивает последовательность в одно значение (fold/reduce)
  template <typename T2>
  T2 Reduce(T2 (*reduce_func)(const T2 &, const T &),
            const T2 &start_value) const {
    T2 result = start_value;

    IEnumerator<T> *it = this->GetEnumerator();
    if (it != nullptr) {
      while (it->HasNext()) {
        result = reduce_func(result, it->GetCurrent());
        it->MoveNext();
      }
      delete it;
    } else {
      for (int index = 0; index < this->GetLength(); index++) {
        result = reduce_func(result, this->Get(index));
      }
    }

    return result;
  }
};

#endif // SEQUENCE_H