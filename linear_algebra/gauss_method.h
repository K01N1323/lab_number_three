#ifndef SECONDFILE_H
#define SECONDFILE_H

template <class T>
class matrixes;

template <class T>
class gauss_method {
  public:
    friend class matrixes<T>;
    gauss_method(T *, int);
    T *GetMatrix() const;
    T GetDet() const;
    void TakeReverse();
    ~gauss_method();
    T *Solve(T *b);
    void SolveForTests(T *b);

  private:
    int n;
    T *matrix;
    T *conmatrix;
    T det;
    int swaps;
    int neededrow(int skip) const;
    void swaprows(int);
    void divisionrow(int);
    void subtraction(int);
    void triangle();
    void obrat();
    void reverse();
    T GetElement(int, int) const;
};

#include "../src/gauss_method.tpp"

#endif // secondfile
