#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

#include <stdexcept>

template <class T> class DynamicArray {
private:
    T *items;
    int size;

public:
    DynamicArray(int size) : DynamicArray(nullptr, size) {}

    DynamicArray(T *items, int count) {
        this->size = count;
        this->items = new T[size];

        if (items == nullptr) {
            for (int element = 0; element < size; element++) {
                this->items[element] = T();
            }
        } else {
            for (int element = 0; element < size; element++) {
                this->items[element] = items[element];
            }
        }
    }

    DynamicArray(const DynamicArray<T> &dynamic_array) {
        this->size = dynamic_array.size;
        this->items = new T[size];

        for (int element = 0; element < size; element++) {
            items[element] = dynamic_array.items[element];
        }
    }

    int GetSize() const { 
        return this->size; 
    }

    const T &Get(int index) const {
        if (index < 0 || index >= size) {
            throw std::out_of_range("Индекс невалиден");
        }
        
        return items[index];
    }

    void Set(int index, T value) {
        if (index < 0 || index >= size) {
            throw std::out_of_range("Индекс невалиден");
        }
        
        items[index] = value;
    }

    void Resize(int new_size) {
        T *new_items = new T[new_size];
        int elements_to_copy = (new_size < size) ? new_size : size;
        
        for (int element = 0; element < elements_to_copy; element++) {
            new_items[element] = items[element];
        }
        
        delete[] items;
        this->size = new_size;
        this->items = new_items;
    }

    ~DynamicArray() { 
        delete[] items; 
    }
};

#endif // DYNAMICARRAY_H