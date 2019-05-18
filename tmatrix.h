#ifndef TMATRIX_H
#define TMATRIX_H
#include <tvector.h>


class TMatrix
{
protected:
    std::vector<TVector> data;
    int n, m;
public:
    //дефолтный конструктор
    TMatrix();
    //создание матрицы по размеру
    TMatrix(int n, int m);
    //конструктор копирования
    TMatrix(const TMatrix& rval);
    //количество строк
    inline int row_count() const { return data.size(); }
    //количество столбцов
    inline int col_count() const { return data.back().size(); }
    //изменение размера матрицы
    void resize(int n, int m);
    TMatrix& operator = (const TMatrix& rval);
    //получение строк
    inline TVector operator[](int i) const { return data[i]; }
    inline TVector& operator[](int i) { return data[i]; }
    //операторы
    TMatrix operator - () const;
    TMatrix operator - (const TMatrix& rval) const;
    TMatrix operator + (const TMatrix& rval) const;
    TMatrix operator * (const long double rval) const;
    //матричное умножение
    TMatrix operator * (const TMatrix& rval) const;
    //умножение матрицы на вектор
    TVector operator * (const TVector& rval) const;
    //обращение матрицы
    TMatrix operator !() const;
    //транспонирование матрицы
    TMatrix t() const;
    TMatrix& swap_rows(int i, int j);
    static TMatrix E(int n);
    ~TMatrix();
    void clear();
    void show_matrix();
};

#endif // TMATRIX_H
