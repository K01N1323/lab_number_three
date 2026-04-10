#ifndef IMMUTABLEARRAYSEQUENCE_H
#define IMMUTABLEARRAYSEQUENCE_H

#include "ArraySequence.h"

template <class T> class ImmutableArraySequence : public ArraySequence<T> {
public:
    ImmutableArraySequence() : ArraySequence<T>() {}
    
    ImmutableArraySequence(T *items, int count) : ArraySequence<T>(items, count) {}
    
    ImmutableArraySequence(const ImmutableArraySequence<T> &other) : ArraySequence<T>(other) {}

    ImmutableArraySequence<T> *Instance() override {
        return new ImmutableArraySequence<T>(*this);
    }

    ImmutableArraySequence<T> *CreateEmpty() const override {
        return new ImmutableArraySequence<T>();
    }
};

#endif // IMMUTABLEARRAYSEQUENCE_H