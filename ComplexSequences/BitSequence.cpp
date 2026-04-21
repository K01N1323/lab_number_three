#include "BitSequence.h"

Bit::Bit(bool v) : value(v) {}

Bit::Bit(int v) : value(v != 0) {}

bool Bit::GetValue() const { 
    return value; 
}

Bit Bit::operator&(const Bit &other) const {
    return Bit(this->value && other.value);
}

Bit Bit::operator|(const Bit &other) const {
    return Bit(this->value || other.value);
}

Bit Bit::operator^(const Bit &other) const {
    return Bit(this->value != other.value);
}

Bit Bit::operator~() const { 
    return Bit(!this->value); 
}

bool Bit::operator==(const Bit &other) const { 
    return this->value == other.value; 
}

std::ostream &operator<<(std::ostream &os, const Bit &b) { 
    os << b.value;
    return os;
}


BitSequence::BitSequence() : MutableArraySequence<Bit>() {}

BitSequence *BitSequence::BitAnd(Sequence<Bit> *other) const {
    BitSequence *NewBitSequence = new BitSequence();

    int ElementsToChange = (this->GetLength() < other->GetLength())
                                 ? this->GetLength()
                                 : other->GetLength();

    for (int index = 0; index < ElementsToChange; index++) {
        NewBitSequence->append(this->get(index) & other->get(index));
    }

    return NewBitSequence;
}

BitSequence *BitSequence::BitOr(Sequence<Bit> *other) const {
    BitSequence *NewBitSequence = new BitSequence();

    int ElementsToChange = (this->GetLength() < other->GetLength())
                                 ? this->GetLength()
                                 : other->GetLength();

    for (int index = 0; index < ElementsToChange; index++) {
        NewBitSequence->append(this->get(index) | other->get(index));
    }

    return NewBitSequence;
}

BitSequence *BitSequence::BitXor(Sequence<Bit> *other) const {
    BitSequence *NewBitSequence = new BitSequence();

    int ElementsToChange = (this->GetLength() < other->GetLength())
                                 ? this->GetLength()
                                 : other->GetLength();

    for (int index = 0; index < ElementsToChange; index++) {
        NewBitSequence->append(this->get(index) ^ other->get(index));
    }

    return NewBitSequence;
}

BitSequence *BitSequence::BitNot() const {
    BitSequence *NewBitSequence = new BitSequence();

    for (int index = 0; index < this->GetLength(); index++) {
        NewBitSequence->append(~this->get(index));
    }

    return NewBitSequence;
}
