#include <iostream>

#include "../Sequences/MutableArraySequence.h"

class Bit {
private:
    bool value;

public:
    Bit(bool v = false) : value(v) {}
    
    Bit(int v) : value(v != 0) {}

    bool GetValue() const { 
        return value; 
    }

    Bit operator&(const Bit &other) const {
        return Bit(this->value && other.value);
    }
    
    Bit operator|(const Bit &other) const {
        return Bit(this->value || other.value);
    }
    
    Bit operator^(const Bit &other) const {
        return Bit(this->value != other.value);
    }
    
    Bit operator~() const { 
        return Bit(!this->value); 
    }

    bool operator==(const Bit &other) const { 
        return this->value == other.value; 
    }

    friend std::ostream &operator<<(std::ostream &os, const Bit &b) { 
        os << b.value;
        return os;
    }
};

class BitSequence : public MutableArraySequence<Bit> {
public:
    BitSequence() : MutableArraySequence<Bit>() {}

    BitSequence *BitAnd(Sequence<Bit> *other) const {
        BitSequence *NewBitSequence = new BitSequence();

        int ElementsToChange = (this->GetLength() < other->GetLength())
                                     ? this->GetLength()
                                     : other->GetLength();

        for (int index = 0; index < ElementsToChange; index++) {
            NewBitSequence->append(this->get(index) & other->get(index));
        }

        return NewBitSequence;
    }
    
    BitSequence *BitOr(Sequence<Bit> *other) const {
        BitSequence *NewBitSequence = new BitSequence();

        int ElementsToChange = (this->GetLength() < other->GetLength())
                                     ? this->GetLength()
                                     : other->GetLength();

        for (int index = 0; index < ElementsToChange; index++) {
            NewBitSequence->append(this->get(index) | other->get(index));
        }

        return NewBitSequence;
    }
    
    BitSequence *BitXor(Sequence<Bit> *other) const {
        BitSequence *NewBitSequence = new BitSequence();

        int ElementsToChange = (this->GetLength() < other->GetLength())
                                     ? this->GetLength()
                                     : other->GetLength();

        for (int index = 0; index < ElementsToChange; index++) {
            NewBitSequence->append(this->get(index) ^ other->get(index));
        }

        return NewBitSequence;
    }
    
    BitSequence *BitNot() const {
        BitSequence *NewBitSequence = new BitSequence();

        for (int index = 0; index < this->GetLength(); index++) {
            NewBitSequence->append(~this->get(index));
        }

        return NewBitSequence;
    }
};
