#ifndef ARRAYSEQUENCE_H
#define ARRAYSEQUENCE_H

#include <stdexcept>

#include "../BaseStructures/DynamicArray.h"
#include "Sequence.h"

template <class T> class ArraySequence : public Sequence<T> {
protected:
  DynamicArray<T> *items;

public:
  ArraySequence() { this->items = new DynamicArray<T>(0); }

  ArraySequence(T *items, int count) {
    if (items == nullptr || count == 0) {
      this->items = new DynamicArray<T>(0);
    } else {
      this->items = new DynamicArray<T>(items, count);
    }
  }

  ArraySequence(const ArraySequence<T> &other) {
    this->items = new DynamicArray<T>(*other.items);
  }

  virtual ~ArraySequence() override { delete this->items; }

  virtual ArraySequence<T> *instance() override = 0;

  virtual ArraySequence<T> *CreateEmpty() const override = 0;

  const T &GetFirst() const override { return this->items->get(0); }

  const T &GetLast() const override {
    return this->items->get(this->items->GetSize() - 1);
  }

  const T &get(int index) const override { return this->items->get(index); }

  int GetLength() const override { return this->items->GetSize(); }

  IEnumerator<T> *GetEnumerator() const override { return nullptr; }

  Sequence<T> *append(T item) override {
    ArraySequence<T> *target = this->instance();

    target->items->resize(target->items->GetSize() + 1);
    target->items->set(target->items->GetSize() - 1, item);

    return target;
  }

  Sequence<T> *prepend(T item) override {
    ArraySequence<T> *target = this->instance();
    int size = target->items->GetSize();

    target->items->resize(size + 1);

    for (int i = size; i > 0; i--) {
      target->items->set(i, target->items->get(i - 1));
    }

    target->items->set(0, item);
    return target;
  }

  Sequence<T> *InsertAt(T item, int index) override {
    if (index < 0 || index > this->GetLength()) {
      throw std::out_of_range("Индекс невалиден");
    }

    ArraySequence<T> *target = this->instance();
    int size = target->items->GetSize();

    target->items->resize(size + 1);

    for (int i = size; i > index; i--) {
      target->items->set(i, target->items->get(i - 1));
    }

    target->items->set(index, item);
    return target;
  }
  // сеттер по индексу
  void set(int index, const T &item) override { this->items->set(index, item); }

  Sequence<T> *concat(Sequence<T> *list) override {
    Sequence<T> *result = this->CreateEmpty();

    for (int i = 0; i < this->GetLength(); i++) {
      Sequence<T> *OldPtr = result;
      result = result->append(this->get(i));
      if (result != OldPtr)
        delete OldPtr;
    }

    for (int i = 0; i < list->GetLength(); i++) {
      Sequence<T> *OldPtr = result;
      result = result->append(list->get(i));
      if (result != OldPtr)
        delete OldPtr;
    }

    return result;
  }

  Sequence<T> *GetSubsequence(int StartIndex, int EndIndex) const override {
    if (StartIndex < 0 || StartIndex >= this->GetLength() || EndIndex < 0 ||
        EndIndex >= this->GetLength() || StartIndex > EndIndex) {
      throw std::out_of_range("Индексы невалидны");
    }

    Sequence<T> *result = this->CreateEmpty();

    for (int index = StartIndex; index <= EndIndex; index++) {
      Sequence<T> *OldPtr = result;

      result = result->append(this->get(index));

      if (result != OldPtr) {
        delete OldPtr;
      }
    }

    return result;
  }

  Sequence<T> *map(T (*mapper)(const T &)) const override {
    Sequence<T> *result = this->CreateEmpty();

    for (int index = 0; index < this->GetLength(); index++) {
      Sequence<T> *OldPtr = result;

      result = result->append(mapper(this->get(index)));

      if (result != OldPtr) {
        delete OldPtr;
      }
    }

    return result;
  }

  Sequence<T> *where(bool (*where)(const T &)) const override {
    Sequence<T> *result = this->CreateEmpty();

    for (int index = 0; index < this->GetLength(); index++) {
      if (where(this->get(index))) {
        Sequence<T> *OldPtr = result;

        result = result->append(this->get(index));

        if (result != OldPtr) {
          delete OldPtr;
        }
      }
    }

    return result;
  }
};

#endif // ARRAYSEQUENCE_H