#ifndef LIST_SEQUENCE_H
#define LIST_SEQUENCE_H
#ifndef LISTSEQUENCE_H
#define LISTSEQUENCE_H

#include <stdexcept>

#include "../BaseStructures/LinkedList.h"
#include "Sequence.h"

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

  virtual ListSequence<T> *instance() override = 0;

  virtual ListSequence<T> *CreateEmpty() const override = 0;

  const T &GetFirst() const override { return this->items->GetFirst(); }

  const T &GetLast() const override { return this->items->GetLast(); }

  const T &get(int index) const override { return this->items->get(index); }

  int GetLength() const override { return this->items->GetLength(); }

  IEnumerator<T> *GetEnumerator() const override {
    return this->items->GetEnumerator();
  }

  Sequence<T> *append(T item) override {
    ListSequence<T> *target = this->instance();

    target->items->append(item);

    return target;
  }

  Sequence<T> *prepend(T item) override {
    ListSequence<T> *target = this->instance();

    target->items->prepend(item);

    return target;
  }

  Sequence<T> *InsertAt(T item, int index) override {
    ListSequence<T> *target = this->instance();

    target->items->InsertAt(item, index);

    return target;
  }
  // сеттер по индексу
  void set(int index, const T &item) override { this->items->set(index, item); }

  Sequence<T> *concat(Sequence<T> *list) override {
    ListSequence<T> *NewList = this->CreateEmpty();

    for (int i = 0; i < this->GetLength(); i++) {
      NewList->items->append(this->get(i));
    }

    IEnumerator<T> *it = list->GetEnumerator();
    if (it != nullptr) {
      while (it->HasNext()) {
        NewList->items->append(it->GetCurrent());
        it->MoveNext();
      }
      delete it;
    } else {
      for (int i = 0; i < list->GetLength(); i++) {
        NewList->items->append(list->get(i));
      }
    }

    return NewList;
  }

  Sequence<T> *GetSubsequence(int StartIndex, int EndIndex) const override {
    LinkedList<T> *RawList = this->items->GetSubList(StartIndex, EndIndex);
    ListSequence<T> *NewList = this->CreateEmpty();

    delete NewList->items;
    NewList->items = RawList;

    return NewList;
  }

  Sequence<T> *map(T (*mapper)(const T &)) const override {
    ListSequence<T> *NewList = this->CreateEmpty();

    for (int index = 0; index < this->GetLength(); index++) {
      NewList->items->append(mapper(this->get(index)));
    }

    return NewList;
  }

  Sequence<T> *where(bool (*where)(const T &)) const override {
    ListSequence<T> *NewList = this->CreateEmpty();

    for (int index = 0; index < this->GetLength(); index++) {
      if (where(this->get(index))) {
        NewList->items->append(this->get(index));
      }
    }

    return NewList;
  }
};

#endif // LISTSEQUENCE_H
#endif // LIST_SEQUENCE_H