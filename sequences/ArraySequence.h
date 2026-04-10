#ifndef ARRAYSEQUENCE_H
#define ARRAYSEQUENCE_H

#include <stdexcept>

#include "DynamicArray.h"
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

  virtual ArraySequence<T> *Instance() = 0;

  virtual ArraySequence<T> *CreateEmpty() const = 0;

  const T &GetFirst() const override { return this->items->Get(0); }

  const T &GetLast() const override {
    return this->items->Get(this->items->GetSize() - 1);
  }

  const T &Get(int index) const override { return this->items->Get(index); }

  int GetLength() const override { return this->items->GetSize(); }

  IEnumerator<T> *GetEnumerator() const override { return nullptr; }

  Sequence<T> *Append(T item) override {
    ArraySequence<T> *target = this->Instance();

    target->items->Resize(target->items->GetSize() + 1);
    target->items->Set(target->items->GetSize() - 1, item);

    return target;
  }

  Sequence<T> *Prepend(T item) override {
    ArraySequence<T> *target = this->Instance();
    int size = target->items->GetSize();

    target->items->Resize(size + 1);

    for (int i = size; i > 0; i--) {
      target->items->Set(i, target->items->Get(i - 1));
    }

    target->items->Set(0, item);
    return target;
  }

  Sequence<T> *InsertAt(T item, int index) override {
    if (index < 0 || index > this->GetLength()) {
      throw std::out_of_range("Индекс невалиден");
    }

    ArraySequence<T> *target = this->Instance();
    int size = target->items->GetSize();

    target->items->Resize(size + 1);

    for (int i = size; i > index; i--) {
      target->items->Set(i, target->items->Get(i - 1));
    }

    target->items->Set(index, item);
    return target;
  }
  // сеттер по индексу
  void Set(int index, const T &item) override { this->items->Set(index, item); }

  Sequence<T> *Concat(Sequence<T> *list) override {
    Sequence<T> *result = this->CreateEmpty();

    for (int i = 0; i < this->GetLength(); i++) {
      Sequence<T> *old_ptr = result;
      result = result->Append(this->Get(i));
      if (result != old_ptr)
        delete old_ptr;
    }

    for (int i = 0; i < list->GetLength(); i++) {
      Sequence<T> *old_ptr = result;
      result = result->Append(list->Get(i));
      if (result != old_ptr)
        delete old_ptr;
    }

    return result;
  }

  Sequence<T> *GetSubsequence(int start_index, int end_index) const override {
    if (start_index < 0 || start_index >= this->GetLength() || end_index < 0 ||
        end_index >= this->GetLength() || start_index > end_index) {
      throw std::out_of_range("Индексы невалидны");
    }

    Sequence<T> *result = this->CreateEmpty();

    for (int index = start_index; index <= end_index; index++) {
      Sequence<T> *old_ptr = result;

      result = result->Append(this->Get(index));

      if (result != old_ptr) {
        delete old_ptr;
      }
    }

    return result;
  }

  Sequence<T> *Map(T (*mapper)(const T &)) const override {
    Sequence<T> *result = this->CreateEmpty();

    for (int index = 0; index < this->GetLength(); index++) {
      Sequence<T> *old_ptr = result;

      result = result->Append(mapper(this->Get(index)));

      if (result != old_ptr) {
        delete old_ptr;
      }
    }

    return result;
  }

  Sequence<T> *Where(bool (*where)(const T &)) const override {
    Sequence<T> *result = this->CreateEmpty();

    for (int index = 0; index < this->GetLength(); index++) {
      if (where(this->Get(index))) {
        Sequence<T> *old_ptr = result;

        result = result->Append(this->Get(index));

        if (result != old_ptr) {
          delete old_ptr;
        }
      }
    }

    return result;
  }
};

#endif // ARRAYSEQUENCE_H