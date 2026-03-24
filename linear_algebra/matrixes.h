#ifndef MATRIXES_H
#define MATRIXES_H

template <class T>
class matrixes {
  public:
    matrixes(int, int);
    matrixes(T *, int, int);
    T GetIJ(int, int) const;
    T Get(int) const;
    T *GetMatrix() const;
    T GetDet() const;
    T *GetInverseMatrix() const;
    int GetRows() const;
    int GetCols() const;
    void MakeGilbert();
    void MakeOnes();
    void MakeRandomNormal();
    const matrixes<T> operator+(const matrixes<T> &rv) const;
    const matrixes<T> operator-(const matrixes<T> &rv) const;
    const matrixes<T> operator*(const matrixes<T> &rv) const;
    ~matrixes();

  private:
    int rows;
    int cols;
    T *matrix;
    const T SubstrRowCol(const matrixes<T> &rv, int, int) const;
};

#include "../src/matrixes.tpp"

#endif // matrixes_h