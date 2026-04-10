#include <iostream>

#include "MutableArraySequence.h"

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

    BitSequence *And(Sequence<Bit> *other) const {
        BitSequence *new_bit_sequence = new BitSequence();

        int elements_to_change = (this->GetLength() < other->GetLength())
                                     ? this->GetLength()
                                     : other->GetLength();

        for (int index = 0; index < elements_to_change; index++) {
            new_bit_sequence->Append(this->Get(index) & other->Get(index));
        }

        return new_bit_sequence;
    }
    
    BitSequence *Or(Sequence<Bit> *other) const {
        BitSequence *new_bit_sequence = new BitSequence();

        int elements_to_change = (this->GetLength() < other->GetLength())
                                     ? this->GetLength()
                                     : other->GetLength();

        for (int index = 0; index < elements_to_change; index++) {
            new_bit_sequence->Append(this->Get(index) | other->Get(index));
        }

        return new_bit_sequence;
    }
    
    BitSequence *Xor(Sequence<Bit> *other) const {
        BitSequence *new_bit_sequence = new BitSequence();

        int elements_to_change = (this->GetLength() < other->GetLength())
                                     ? this->GetLength()
                                     : other->GetLength();

        for (int index = 0; index < elements_to_change; index++) {
            new_bit_sequence->Append(this->Get(index) ^ other->Get(index));
        }

        return new_bit_sequence;
    }
    
    BitSequence *Not() const {
        BitSequence *new_bit_sequence = new BitSequence();

        for (int index = 0; index < this->GetLength(); index++) {
            new_bit_sequence->Append(~this->Get(index));
        }

        return new_bit_sequence;
    }
};
