#ifndef MUTABLELISTSEQUENCE_H
#define MUTABLELISTSEQUENCE_H

#include "LinkedList.h"
#include "ListSequence.h"

// Изменяемая последовательность на основе связного списка (Instance возвращает this)
template <class T> class MutableListSequence : public ListSequence<T> {
public:
  MutableListSequence() : ListSequence<T>() {}
  MutableListSequence(T *items, int count) : ListSequence<T>(items, count) {}
  MutableListSequence(const MutableListSequence<T> &other)
      : ListSequence<T>(other) {}
  MutableListSequence(const LinkedList<T> &list) : ListSequence<T>(list) {}

  ListSequence<T> *Instance() override { return this; }

  ListSequence<T> *CreateEmpty() const override {
    return new MutableListSequence<T>();
  }
};

#endif // MUTABLELISTSEQUENCE_H