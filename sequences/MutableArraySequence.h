#ifndef MUTABLEARRAYSEQUENCE_H
#define MUTABLEARRAYSEQUENCE_H

#include "ArraySequence.h"

template <class T> class MutableArraySequence : public ArraySequence<T> {
public:
    MutableArraySequence() : ArraySequence<T>() {}
    
    MutableArraySequence(T *items, int count) : ArraySequence<T>(items, count) {}
    
    MutableArraySequence(const MutableArraySequence<T> &other) : ArraySequence<T>(other) {}

    MutableArraySequence<T> *instance() override { 
        return this; 
    }

    MutableArraySequence<T> *CreateEmpty() const override {
        return new MutableArraySequence<T>();
    }
};

#endif // MUTABLEARRAYSEQUENCE_H