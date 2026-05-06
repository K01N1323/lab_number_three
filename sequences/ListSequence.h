#ifndef LIST_SEQUENCE_H
#define LIST_SEQUENCE_H

#include <stdexcept>

#include "LinkedList.h"
#include "Sequence.h"

// Абстрактный базовый класс для последовательностей на основе связного списка
template <class T> class ListSequence : public Sequence<T> {
protected:
  LinkedList<T> *items;

public:
  ListSequence() { this->items = new LinkedList<T>(); }

  ListSequence(T *items, int count) {
    if (items == nullptr) {
      this->items = new LinkedList<T>();
    } else {
      this->items = new LinkedList<T>(items, count);
    }
  }

  ListSequence(const ListSequence<T> &other) {
    this->items = new LinkedList<T>(*other.items);
  }

  ListSequence(const LinkedList<T> &list) {
    this->items = new LinkedList<T>(list);
  }

  virtual ~ListSequence() override { delete this->items; }

  virtual ListSequence<T> *Instance() override = 0;
  virtual ListSequence<T> *CreateEmpty() const override = 0;

  const T &GetFirst() const override { return this->items->GetFirst(); }
  const T &GetLast() const override { return this->items->GetLast(); }
  const T &Get(int index) const override { return this->items->Get(index); }
  void Set(int index, const T &value) override { this->items->Set(index, value); }
  int GetLength() const override { return this->items->GetLength(); }

  IEnumerator<T> *GetEnumerator() const override {
    return this->items->GetEnumerator();
  }

  Sequence<T> *Append(const T &item) override {
    ListSequence<T> *target = this->Instance();
    target->items->Append(item);
    return target;
  }

  Sequence<T> *Prepend(const T &item) override {
    ListSequence<T> *target = this->Instance();
    target->items->Prepend(item);
    return target;
  }

  Sequence<T> *InsertAt(const T &item, int index) override {
    ListSequence<T> *target = this->Instance();
    target->items->InsertAt(item, index);
    return target;
  }

  Sequence<T> *Concat(Sequence<T> *list) override {
    ListSequence<T> *new_list = this->CreateEmpty();

    // Обход через итератор — O(N) вместо O(N^2)
    IEnumerator<T> *it_self = this->GetEnumerator();
    if (it_self != nullptr) {
      while (it_self->HasNext()) {
        new_list->items->Append(it_self->GetCurrent());
        it_self->MoveNext();
      }
      delete it_self;
    }

    IEnumerator<T> *it = list->GetEnumerator();
    if (it != nullptr) {
      while (it->HasNext()) {
        new_list->items->Append(it->GetCurrent());
        it->MoveNext();
      }
      delete it;
    } else {
      for (int i = 0; i < list->GetLength(); i++) {
        new_list->items->Append(list->Get(i));
      }
    }

    return new_list;
  }

  Sequence<T> *GetSubsequence(int start_index, int end_index) const override {
    LinkedList<T> *raw_list = this->items->GetSubList(start_index, end_index);
    ListSequence<T> *new_list = this->CreateEmpty();

    delete new_list->items;
    new_list->items = raw_list;

    return new_list;
  }
};

#endif // LISTSEQUENCE_H