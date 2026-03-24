#ifndef DEFAULT_GAUSS_H
#define DEFAULT_GAUSS_H

template <class T> class matrixes;

template <class T>
class default_gauss {
  public:
    friend class matrixes<T>;
    default_gauss(T *, int);
    T *GetMatrix() const;
    T GetDet() const;
    void TakeReverse();
    ~default_gauss();
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

#include "../src/default_gauss.tpp"

#endif // default gauss h