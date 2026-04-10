#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdexcept>

#include "IEnumerator.h"

template <class T> class Sequence {
public:
  virtual ~Sequence() {}

  virtual const T &GetFirst() const = 0;

  virtual const T &GetLast() const = 0;

  virtual const T &Get(int index) const = 0;

  virtual int GetLength() const = 0;

  virtual Sequence<T> *GetSubsequence(int startIndex, int endIndex) const = 0;

  virtual Sequence<T> *Append(T item) = 0;

  virtual Sequence<T> *Prepend(T item) = 0;

  virtual Sequence<T> *InsertAt(T item, int index) = 0;

  virtual void Set(int index, const T &) = 0;

  virtual Sequence<T> *Concat(Sequence<T> *list) = 0;

  const T &operator[](int index) const { return this->Get(index); }

  virtual Sequence<T> *Instance() = 0;

  virtual IEnumerator<T> *GetEnumerator() const = 0;

  virtual Sequence<T> *Map(T (*mapper)(const T &)) const = 0;

  virtual Sequence<T> *Where(bool (*where)(const T &)) const = 0;

  // Сворачивает последовательность в одно значение
  template <typename T2>
  T2 Reduce(T2 (*reduce_func)(const T2 &, const T &),
            const T2 &start_value) const {
    T2 result = start_value;

    for (int index = 0; index < this->GetLength(); index++) {
      result = reduce_func(result, this->Get(index));
    }

    return result;
  }
};

#endif // SEQUENCE_H