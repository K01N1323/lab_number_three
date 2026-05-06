#ifndef IMMUTABLELISTSEQUENCE_H
#define IMMUTABLELISTSEQUENCE_H

#include "LinkedList.h"
#include "ListSequence.h"

// Неизменяемая последовательность на основе связного списка (Instance возвращает копию)
template <class T> class ImmutableListSequence : public ListSequence<T> {
public:
  ImmutableListSequence() : ListSequence<T>() {}
  ImmutableListSequence(T *items, int count) : ListSequence<T>(items, count) {}
  ImmutableListSequence(const ImmutableListSequence<T> &other)
      : ListSequence<T>(other) {}
  ImmutableListSequence(const LinkedList<T> &list) : ListSequence<T>(list) {}

  ListSequence<T> *Instance() override {
    return new ImmutableListSequence<T>(*this);
  }

  ListSequence<T> *CreateEmpty() const override {
    return new ImmutableListSequence<T>();
  }
};

#endif // IMMUTABLELISTSEQUENCE_H