#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <stdexcept>

#include "../Sequences/IEnumerator.h"

template <class T> class LinkedList {
private:
  struct Node {
    T data;
    Node *next;

    Node(T value) {
      this->data = value;
      this->next = nullptr;
    }
  };

  Node *head;
  Node *tail;
  int size;

public:
  LinkedList() : LinkedList(nullptr, 0) {}

  LinkedList(T *items, int count) {
    this->head = nullptr;
    this->tail = nullptr;
    this->size = 0;

    for (int index = 0; index < count; index++) {
      Node *NewNode = new Node(items[index]);

      if (head == nullptr) {
        head = NewNode;
        tail = NewNode;
      } else {
        tail->next = NewNode;
        tail = NewNode;
      }
      size++;
    }
  }

  LinkedList(const LinkedList<T> &list) : LinkedList() {
    Node *current = list.head;

    while (current != nullptr) {
      Node *NewNode = new Node(current->data);

      if (head == nullptr) {
        head = NewNode;
        tail = NewNode;
      } else {
        tail->next = NewNode;
        tail = NewNode;
      }

      current = current->next;
    }

    this->size = list.size;
  }

  const T &GetFirst() const {
    if (head == nullptr) {
      throw std::out_of_range("Первого элемента не существует");
    }

    return head->data;
  }

  const T &GetLast() const {
    if (tail == nullptr) {
      throw std::out_of_range("Последнего элемента не существует");
    }

    return tail->data;
  }

  const T &get(int index) const {
    if (index >= size || index < 0) {
      throw std::out_of_range("Элемента с таким индексом не существует");
    }

    Node *current = head;

    for (int ind = 0; ind < index; ind++) {
      current = current->next;
    }

    return current->data;
  }

  LinkedList<T> *GetSubList(int StartIndex, int EndIndex) {
    if (StartIndex < 0 || EndIndex < 0 || EndIndex >= size ||
        StartIndex >= size || EndIndex < StartIndex) {
      throw std::out_of_range("Индексы невалидны для данного списка");
    }

    LinkedList<T> *SubList = new LinkedList<T>();
    Node *current = this->head;

    for (int index = 0; index < StartIndex; index++) {
      current = current->next;
    }

    Node **CurrentPtr = &(SubList->head);

    for (int index = StartIndex; index <= EndIndex; index++) {
      Node *NewNode = new Node(current->data);

      *CurrentPtr = NewNode;
      CurrentPtr = &(NewNode->next);

      SubList->tail = NewNode;
      SubList->size++;

      current = current->next;
    }

    return SubList;
  }

  int GetLength() const { return size; }

  void append(T item) {
    Node *NewNode = new Node(item);

    if (head == nullptr) {
      head = NewNode;
      tail = NewNode;
    } else {
      tail->next = NewNode;
      tail = NewNode;
    }

    size++;
  }

  void prepend(T item) {
    Node *NewNode = new Node(item);

    if (this->head == nullptr) {
      this->head = NewNode;
      this->tail = NewNode;
    } else {
      NewNode->next = this->head;
      this->head = NewNode;
    }

    size++;
  }

  void InsertAt(T item, int index) {
    if (index >= size || index < 0) {
      throw std::out_of_range("Индекс вне списка");
    }

    if (index == 0) {
      prepend(item);
      return;
    }

    Node *NewNode = new Node(item);
    Node *current = head;

    for (int ind = 0; ind < (index - 1); ind++) {
      current = current->next;
    }

    NewNode->next = current->next;
    current->next = NewNode;
    size++;
  }

  LinkedList<T> *concat(LinkedList<T> *list) {
    LinkedList<T> *SubList = new LinkedList<T>(*this);
    Node *current = list->head;

    for (int index = 0; index < list->size; index++) {
      SubList->append(current->data);
      current = current->next;
    }

    return SubList;
  }

  ~LinkedList() {
    Node *current = this->head;

    while (current != nullptr) {
      Node *LocalCurrent = current->next;
      delete current;
      current = LocalCurrent;
    }
  }

  class LinkedListEnumerator : public IEnumerator<T> {
  private:
    Node *current;

  public:
    LinkedListEnumerator(Node *head) { current = head; }

    const T &GetCurrent() const override { return current->data; }

    void MoveNext() override {
      if (current != nullptr) {
        current = current->next;
      }
    }

    bool HasNext() const override { return current != nullptr; }
  };

  IEnumerator<T> *GetEnumerator() const {
    return new LinkedListEnumerator(this->head);
  }
};

#endif // LINKEDLIST_H