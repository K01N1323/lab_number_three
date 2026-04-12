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

  DynamicArray(const DynamicArray<T> &DynamicArraySource) {
    this->size = DynamicArraySource.size;
    this->items = new T[size];

    for (int element = 0; element < size; element++) {
      items[element] = DynamicArraySource.items[element];
    }
  }

  int GetSize() const { return this->size; }

  const T &get(int index) const {
    if (index < 0 || index >= size) {
      throw std::out_of_range("Индекс невалиден");
    }

    return items[index];
  }

  void set(int index, T value) {
    if (index < 0 || index >= size) {
      throw std::out_of_range("Индекс невалиден");
    }

    items[index] = value;
  }

  void resize(int NewSize) {
    T *NewItems = new T[NewSize];
    int ElementsToCopy = (NewSize < size) ? NewSize : size;

    for (int element = 0; element < ElementsToCopy; element++) {
      NewItems[element] = items[element];
    }

    delete[] items;
    this->size = NewSize;
    this->items = NewItems;
  }

  ~DynamicArray() { delete[] items; }
};

#endif // DYNAMICARRAY_H