#pragma once//---------------------------------------------------------------------------
#include "model.h"
#include <tuple>
#include <math.h>
#include <vector>
#include <fstream>
#include "generator.h"
using namespace std;
using namespace arma;
//---------------------------------------------------------------------------

class TArenstorfModel  : public TModel
{
    protected:
        static const long double m;
        long double D1, D2;
    public:
        TArenstorfModel(  );
        void getRight( const vec& X, long double t, vec& Y );
};

class TSatellite  : public TModel
{
protected:
    const long double Re{6371000.0L};
    const long double nu{398600.436e9}; //грав параметр
    const long double Rorb{23222000.0L}; //m
    const long double W{7.292115E-5}; //угловая скорость вращения Земли
    int satellite_number{0}, satellite_count{27};
    //параметры орибты
    long double orb_height{Re + Rorb}; //высота от центра Земли (м)
    const int planes_number{3}; // количество плоскостей
    long double period{50685.0L}; //период обращения (с)
    long double i{0.0L/*0.977384L*/};//наклонение (рад)
    long double eccentricity{0.1L};//безразмерная (доли единицы)
    long double u{M_PI/*0.174533L*/};//аргумент широты (рад)
    long double omega{3.0L*M_PI/2.0L/*3.75246L*/};//долгота восходящего узла (рад)
    vec X_oscul;//вектор положения в оскулирующих элементах
    vec V_oscul;//вектор скоростей в оскулирующих элементах
    //матрицы перехода
    mat A; //матрица перехода между подвижной орбитальной СК к инерц геоцентр СК
    //файл с дистанциями
public:
    mat R;
    long double ch_lat, ch_lon;
    TSatellite(int satellite_number, long double d_u, long double d_omega, long double i);//по номеру спуткника (кратность 9) задаем одну из трех орбит, широту
    void getRight(const vec& X, long double t, vec& Y);
};








