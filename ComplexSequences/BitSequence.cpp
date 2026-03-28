#include <iostream>

#include "MutableArraySequence.h"

// Класс, представляющий отдельный бит (0 или 1)
class Bit {
private:
    bool value;

public:
    // Конструктор по умолчанию или из значения bool
    Bit(bool v = false) : value(v) {}
    
    // Конструктор из int (преобразуется в логическое значение)
    Bit(int v) : value(v != 0) {}

    // Возвращает значение бита
    bool GetValue() const { 
        return value; 
    }

    // Перегрузка оператора побитового И (&)
    Bit operator&(const Bit &other) const {
        return Bit(this->value && other.value);
    }
    
    // Перегрузка оператора побитового ИЛИ (|)
    Bit operator|(const Bit &other) const {
        return Bit(this->value || other.value);
    }
    
    // Перегрузка оператора побитового ИСКЛЮЧАЮЩЕГО ИЛИ (^)
    Bit operator^(const Bit &other) const {
        return Bit(this->value != other.value);
    }
    
    // Перегрузка оператора логического отрицания (~)
    Bit operator~() const { 
        return Bit(!this->value); 
    }

    // Перегрузка оператора сравнения на равенство (==)
    bool operator==(const Bit &other) const { 
        return this->value == other.value; 
    }

    // Перегрузка оператора вывода в поток (<<)
    friend std::ostream &operator<<(std::ostream &os, const Bit &b) { 
        os << b.value;
        return os;
    }
};

// Класс для работы с последовательностью битов
class BitSequence : public MutableArraySequence<Bit> {
public:
    // Конструктор по умолчанию
    BitSequence() : MutableArraySequence<Bit>() {}

    // Выполняет побитовое И с другой последовательностью битов
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
    
    // Выполняет побитовое ИЛИ с другой последовательностью битов
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
    
    // Выполняет побитовое ИСКЛЮЧАЮЩЕЕ ИЛИ с другой последовательностью битов
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
    
    // Инвертирует все биты в последовательности
    BitSequence *Not() const {
        BitSequence *new_bit_sequence = new BitSequence();

        for (int index = 0; index < this->GetLength(); index++) {
            new_bit_sequence->Append(~this->Get(index));
        }

        return new_bit_sequence;
    }
};
