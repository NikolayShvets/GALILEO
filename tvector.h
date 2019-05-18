#ifndef VECTOR_ALG_H
#define VECTOR_ALG_H
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;
class TMatrix;
class TVector
{
protected:
    int n;
    std::vector<long double> data;
public:
    //конструктор дефолтный
    TVector();
    //создание вектора определённого размера
    TVector(int n);
    //конструктор копирования
    TVector(const TVector& rval);
    //получение размера вектора
    inline int size() const {return data.size();}
    //изменение размера вектора
    void resize(int n);
    //получение элемента
    int row_count(){return this->data.size();}
    inline long double operator[](int i) const {return data[i];}
    inline long double& operator[](int i) {return data[i];}
    //присвоение вектору
    TVector& operator = (const TVector& rval);
    //модуль вектора
    long double length() const;
    //операторы
    TVector operator - () const;
    TVector operator - (const long double rval) const;
    TVector operator - (const TVector& rval) const;
    TVector operator + (const long double rval) const;
    TVector operator + (const TVector& rval) const;
    //векторное умнжение векторов
    TVector operator ^ (const TVector& rval) const;
    TVector operator * (const long double rval) const;
    //скалярное умножение векторов
    long double operator * (const TVector& rval) const;
    void show_vector();
    void vector_clear();
    ~TVector();
};

#endif VECTOR_ALG_H
