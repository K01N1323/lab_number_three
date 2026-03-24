#ifndef LU_DECOMPOZITION_H
#define LU_DECOMPOZITION_H

template <class T>
class LU_Decomposition {
  public:
    LU_Decomposition(int n) : LU_Decomposition(nullptr, n) {}
    LU_Decomposition(T *matrix_data, int n);

    ~LU_Decomposition();

    T *GetL() const;
    T *GetU() const;

    int *GetP() const;

    T GetDet() const;

    T *Solve(T *b) const;
    void SolveForTests(T *b) const;

  private:
    int n;

    T *L;
    T *U;
    int *P;
    int swaps;

    T *x;
    T *y;

    int neededrow(int skip) const;

    void swaprows(int step);

    void subtraction(int current);

    void decompose();
};

#include "../src/LU_decompozition.tpp"

#endif // lu decompozition h