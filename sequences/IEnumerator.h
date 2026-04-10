#ifndef IENUMERATOR_H
#define IENUMERATOR_H

// Интерфейс для итератора (перечислителя) элементов
template <class T> class IEnumerator {
public:
    virtual ~IEnumerator() {}

    virtual const T &GetCurrent() const = 0;

    virtual void MoveNext() = 0;

    virtual bool HasNext() const = 0;
};

#endif // IENUMERATOR_H