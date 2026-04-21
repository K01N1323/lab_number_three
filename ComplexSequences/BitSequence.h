#ifndef BIT_SEQUENCE_H
#define BIT_SEQUENCE_H

#include <iostream>
#include "../Sequences/MutableArraySequence.h"

class Bit {
private:
    bool value;

public:
    Bit(bool v = false);
    Bit(int v);

    bool GetValue() const;

    Bit operator&(const Bit &other) const;
    Bit operator|(const Bit &other) const;
    Bit operator^(const Bit &other) const;
    Bit operator~() const;
    bool operator==(const Bit &other) const;

    friend std::ostream &operator<<(std::ostream &os, const Bit &b);
};

class BitSequence : public MutableArraySequence<Bit> {
public:
    BitSequence();

    BitSequence *BitAnd(Sequence<Bit> *other) const;
    BitSequence *BitOr(Sequence<Bit> *other) const;
    BitSequence *BitXor(Sequence<Bit> *other) const;
    BitSequence *BitNot() const;
};

#endif // BIT_SEQUENCE_H