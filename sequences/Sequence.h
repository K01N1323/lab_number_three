#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdexcept>

#include "IEnumerator.h"

template <class T> class Sequence {
public:
  virtual ~Sequence() {}

  virtual const T &GetFirst() const = 0;

  virtual const T &GetLast() const = 0;

  virtual const T &get(int index) const = 0;

  virtual int GetLength() const = 0;

  virtual Sequence<T> *GetSubsequence(int StartIndex, int EndIndex) const = 0;

  virtual Sequence<T> *append(T item) = 0;

  virtual Sequence<T> *prepend(T item) = 0;

  virtual Sequence<T> *InsertAt(T item, int index) = 0;

  virtual void set(int index, const T &) = 0;

  virtual Sequence<T> *concat(Sequence<T> *list) = 0;

  const T &operator[](int index) const { return this->get(index); }

  virtual Sequence<T> *instance() = 0;

  virtual Sequence<T> *CreateEmpty() const = 0;

  virtual IEnumerator<T> *GetEnumerator() const = 0;

  virtual Sequence<T> *map(T (*mapper)(const T &)) const = 0;

  virtual Sequence<T> *where(bool (*where)(const T &)) const = 0;

  // Сворачивает последовательность в одно значение
  template <typename T2>
  T2 reduce(T2 (*ReduceFunc)(const T2 &, const T &),
            const T2 &StartValue) const {
    T2 result = StartValue;

    for (int index = 0; index < this->GetLength(); index++) {
      result = ReduceFunc(result, this->get(index));
    }

    return result;
  }
};

#endif // SEQUENCE_H