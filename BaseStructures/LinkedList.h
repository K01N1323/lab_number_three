#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <stdexcept>

#include "IEnumerator.h"

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
            Node *new_node = new Node(items[index]);
            
            if (head == nullptr) {
                head = new_node;
                tail = new_node;
            } else {
                tail->next = new_node;
                tail = new_node;
            }
            size++;
        }
    }

    LinkedList(const LinkedList<T> &list) : LinkedList() {
        Node *current = list.head;
        
        while (current != nullptr) {
            Node *new_node = new Node(current->data);
            
            if (head == nullptr) {
                head = new_node;
                tail = new_node;
            } else {
                tail->next = new_node;
                tail = new_node;
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

    const T &Get(int index) const {
        if (index >= size || index < 0) {
            throw std::out_of_range("Элемента с таким индексом не существует");
        }
            
        Node *current = head;

        for (int ind = 0; ind < index; ind++) {
            current = current->next;
        }
        
        return current->data;
    }

    LinkedList<T> *GetSubList(int start_index, int end_index) {
        if (start_index < 0 || end_index < 0 || end_index >= size ||
            start_index >= size || end_index < start_index) {
            throw std::out_of_range("Индексы невалидны для данного списка");
        }

        LinkedList<T> *sub_list = new LinkedList<T>();
        Node *current = this->head;

        for (int index = 0; index < start_index; index++) {
            current = current->next;
        }

        Node **current_ptr = &(sub_list->head);

        for (int index = start_index; index <= end_index; index++) {
            Node *new_node = new Node(current->data);

            *current_ptr = new_node;
            current_ptr = &(new_node->next);

            sub_list->tail = new_node;
            sub_list->size++;

            current = current->next;
        }

        return sub_list;
    }

    int GetLength() const { 
        return size; 
    }

    void Append(T item) {
        Node *new_node = new Node(item);
        
        if (head == nullptr) {
            head = new_node;
            tail = new_node;
        } else {
            tail->next = new_node;
            tail = new_node;
        }
        
        size++;
    }

    void Prepend(T item) {
        Node *new_node = new Node(item);
        
        if (this->head == nullptr) {
            this->head = new_node;
            this->tail = new_node;
        } else {
            new_node->next = this->head;
            this->head = new_node;
        }
        
        size++;
    }

    void InsertAt(T item, int index) {
        if (index >= size || index < 0) {
            throw std::out_of_range("Индекс вне списка");
        }
            
        if (index == 0) {
            Prepend(item);
            return;
        }
        
        Node *new_node = new Node(item);
        Node *current = head;
        
        for (int ind = 0; ind < (index - 1); ind++) {
            current = current->next;
        }
        
        new_node->next = current->next;
        current->next = new_node;
        size++;
    }

    LinkedList<T> *Concat(LinkedList<T> *list) {
        LinkedList<T> *sub_list = new LinkedList<T>(*this);
        Node *current = list->head;

        for (int index = 0; index < list->size; index++) {
            sub_list->Append(current->data);
            current = current->next;
        }
        
        return sub_list;
    }

    ~LinkedList() {
        Node *current = this->head;
        
        while (current != nullptr) {
            Node *local_current = current->next;
            delete current;
            current = local_current;
        }
    }

    class LinkedListEnumerator : public IEnumerator<T> {
    private:
        Node *current;
        
    public:
        LinkedListEnumerator(Node *head) { 
            current = head; 
        }
        
        const T &GetCurrent() const override { 
            return current->data; 
        }
        
        void MoveNext() override {
            if (current != nullptr) {
                current = current->next;
            }
        }
        
        bool HasNext() const override { 
            return current != nullptr; 
        }
    };

    IEnumerator<T> *GetEnumerator() const { 
        return new LinkedListEnumerator(this->head); 
    }
};

#endif // LINKEDLIST_H