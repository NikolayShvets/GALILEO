//---------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include "custom.h"
//#include "LA.h"

/*#include "tvector.h"
#include "tmatrix.h"*/
#include "iomanip"
using namespace std;
int rev_sgn(long double val)
{
   return (val<0)?(1):((val>0)?(-1):(0));
}

//---------------------------------------------------------------------------
const long double TArenstorfModel::m  = 0.012277471;
TArenstorfModel::TArenstorfModel() : TModel()
{
    X0.resize(4);
    X0[0] = 0.994;
    X0[1] = 0;
    X0[2] = 0;
    X0[3] = -2.0015851063790825224053786222;
}

//---------------------------------------------------------------------------

void TArenstorfModel::getRight( const TVector& X, long double t, TVector& Y )
{
    Y.resize(4);
    D1 = pow( pow( X[0] + m, 2 ) + pow( X[1], 2 ), 1.5 );
    D2 = pow( pow( X[0] + m - 1, 2 ) + pow( X[1], 2 ), 1.5 );
    Y[0] = X[2];
    Y[1] = X[3];
    Y[2] = X[0] + 2 * X[3] - (1 - m)*(X[0] + m) / D1 - m * (X[0] + m - 1) / D2;
    Y[3] = X[1] - 2 * X[2] - (1 - m)*X[1] / D1 - m * X[1] / D2;
}
//---------------------------------------------------------------------------

TSatellite::TSatellite(int satellite_number,long double d_u, long double d_omega, long double i) : TModel()
{
    cout<<"\nSatellite number "<<satellite_number<<": \n"<<endl;
    /*cout<<"Input chronometr coordinates: "<<endl;
    cout<< "Сhronometr latitude: ";
    cin >> ch_lat;
    cout<< "Сhronometr longtitude: ";
    cin >> ch_lon;*/
    ch_lat = 0.95993108859;//M_PI/2.50L;//1.291L;//Москва
    ch_lon = 0.0L;//0.64577182324//2.523L;
    long double latitude{0.0L}, longtitude{0.0L}, r_vector{0.0L}, v_vector{0.0L};
    this->satellite_number = satellite_number;
    this->u += d_u;
    this->i = i;
    this->omega += d_omega;
    this->A.resize(3,3);
    this->A[0][0] = cos(u)*cos(omega) - sin(u)*sin(omega)*cos(i); this->A[0][1] = -sin(u)*cos(omega) - cos(u)*sin(omega)*cos(i); this->A[0][2] = sin(i)*sin(omega);
    this->A[1][0] = cos(u)*sin(omega) + sin(u)*cos(omega)*cos(i); this->A[1][1] = -sin(u)*sin(omega) + cos(u)*cos(omega)*cos(i); this->A[1][2] = -sin(i)*cos(omega);
    this->A[2][0] = sin(u)*sin(i);                                this->A[2][1] = cos(u)*sin(i) ;                                 this->A[2][2] = cos(i);

    cout<<"Matrix A: "<<endl;
    cout<< "\t" << A[0][0] << "\t" << A[0][1] << "\t" << A[0][2] << endl;
    cout<< "\t" << A[1][0] << "\t" << A[1][1] << "\t" << A[1][2] << endl;
    cout<< "\t" << A[2][0] << "\t" << A[2][1] << "\t" << A[2][2] << endl;

    this->X_oscul.resize(3);
    X_oscul[0] = orb_height;
    X_oscul[1] = 0.0L;
    X_oscul[2] = 0.0L;

    this->V_oscul.resize(3);
    V_oscul[0] = 0.0L;//pow((nu/orb_height), 0.5L);
    V_oscul[1] = pow((nu/orb_height), 0.5L);
    V_oscul[2] = 0.0L;

    //пересчет координат вектора начальных условий
    TVector temp_coordinates(3);
    TVector temp_velosities(3);
    temp_coordinates = this->A*this->X_oscul;
    temp_velosities = this->A*this->V_oscul;

    cout<<fixed;
    cout<<"Osculate elements:"<<endl;
    cout<< "    u: "<<this->u<<" rad"<<endl;
    cout<< "    omega: "<<this->omega<<" rad"<<endl;
    cout<< "    i: "<<this->i<<" rad"<<endl;

    X0.resize(6);
    X0[0] = temp_coordinates[0];
    X0[1] = temp_coordinates[1];
    X0[2] = temp_coordinates[2];

    latitude = atan(X0[2]/sqrt( ( pow(X0[0],2.L) + pow(X0[1], 2.L) )));
    longtitude = atan(X0[1]/X0[0]);
    r_vector = sqrt(pow(X0[0],2.) + pow(X0[1], 2.) + pow(X0[2], 2.));
    v_vector = 3677.233;//sqrt(pow(temp_velosities[0],2.) + pow(temp_velosities[1], 2.) + pow(temp_velosities[2], 2.));

    cout<< "Sfericheskie: "<<endl;
    cout<< "    lat: "<<latitude<<endl;
    cout<< "    lon: "<<longtitude<<endl;
    cout<< "    r_v: "<<r_vector<<endl;
    cout<< "    v_v: "<<v_vector<<endl;

    X0[3] = /*v_vector*cos(latitude)*sin(longtitude);*/temp_velosities[0];
    X0[4] = /*-v_vector*cos(latitude)*cos(longtitude);*/temp_velosities[1];
    X0[5] = /*v_vector*sin(latitude);*/temp_velosities[2];

    cout<<"Coordinates: "<<endl;
    cout<< "    x: "<<X0[0]<<" m"<<endl;
    cout<< "    y: "<<X0[1]<<" m"<<endl;
    cout<< "    z: "<<X0[2]<<" m"<<endl;
    cout<< "    |r-vector| = "<<sqrt(pow(X0[0],2.) + pow(X0[1], 2.) + pow(X0[2], 2.))<<endl;

    cout<<"Velosities: "<<endl;
    cout<< "    Vx: "<<X0[3]<<" m/sec"<<endl;
    cout<< "    Vy: "<<X0[4]<<" m/sec"<<endl;
    cout<< "    Vz: "<<X0[5]<<" m/sec"<<endl;
    cout<< "    |v-vector| = "<<sqrt(pow(X0[3],2.) + pow(X0[4], 2.) + pow(X0[5], 2.))<<endl;
    /*this->latitude = latitude;
    this->longtitude = longtitude;
    //координаты
    this->X0.resize(3);
    X0[0] = this->latitude; //0
    X0[1] = this->longtitude; //0
    X0[2] = Re+Rorb;
    cout<<"Satellite namber "<<satellite_number<<endl;
    cout<<fixed;
    cout<<"geocenter: "<<endl;
    cout<<"fi: "<<X0[0]<<endl;
    cout<<"lambda: "<<X0[1]<<endl;
    cout<<"H: "<<X0[2]<<endl;
    //переход от геоцентрических координат к топоцентрическим

    X0[0] = (Re+Rorb)*cos(latitude)*sin(longtitude);
    X0[1] = (Re+Rorb)*cos(latitude)*cos(longtitude);
    X0[2] = (Re+Rorb)*sin(latitude);
    X0.resize(6);
    //скорость у и z поменялись местами
    X0[3] = 3677.233*cos(latitude)*cos(longtitude);//0.0L;
    X0[4] = 3677.233*cos(latitude)*sin(longtitude);//3677.233*cos(latitude);
    X0[5] = -3677.233*sin(latitude);//-3677.233*sin(latitude);


    cout<<fixed;
    cout<<"topocentr: "<<endl;
    cout<<"X0: "<<X0[0]<<endl;
    cout<<"Y0: "<<X0[1]<<endl;
    cout<<"Z0: "<<X0[2]<<endl;
    cout<<"VX0: "<<X0[3]<<endl;
    cout<<"VY0: "<<X0[4]<<endl;
    cout<<"VZ0: "<<X0[5]<<endl;
    cout<<"speed: "<<endl;
    cout<<"|V| = "<<sqrt(pow(X0[3],2.) + pow(X0[4],2.) + pow(X0[5],2.))<<endl;
    cout<<"|R| = "<<sqrt(pow(X0[0],2.) + pow(X0[1],2.) + pow(X0[2],2.))<<endl;
    cout<<"************************"<<endl;

    //скорости*/
}
//---------------------------------------------------------------------------

void TSatellite::getRight(const TVector &X, long double t, TVector &Y)
{
    Y.resize(6);
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    long double ro = sqrt(pow(X[0],2.) + pow(X[1], 2.) + pow(X[2], 2.));
    Y[3] = -nu*X[0]/pow(ro,3.);
    Y[4] = -nu*X[1]/pow(ro,3.);
    Y[5] = -nu*X[2]/pow(ro,3.);
}
//---------------------------------------------------------------------------
long double TSatellite::do_thing(const TVector &X, long double t)
{
    generator _generator;
    long double distance;
    bool visibility{false};
    TVector chronometr_vector(3);
    TVector temp_v(3);
    temp_v[0] = X[0];
    temp_v[1] = X[1];
    temp_v[2] = X[2];
    //доворачиваем Землю
    //long double temp_lon{ch_lon};
    ch_lon = W*(t - t0);
    //cout<<ch_lon<<endl;
    //перевод координат радиус вектора к часам на повехрности в декартову систему
    chronometr_vector[0] = Re*cos(ch_lat)*cos(ch_lon);
    chronometr_vector[1] = Re*cos(ch_lat)*sin(ch_lon);
    chronometr_vector[2] = Re*sin(ch_lat);
    long double alfa = acos(temp_v*chronometr_vector/(temp_v.length()*chronometr_vector.length()));
    //cout<<"ugol: "<<alfa*180.0L/M_PI<<" lon chronom: "<<ch_lon*180.0L/M_PI<<" time: "<<t - t0<<endl;
    if (alfa < M_PI/3.0L) visibility = true; else visibility = false;
    if (visibility == true)
    {//рассчет расстояния
        //cout<<_generator.white_noise_generator()<<endl;
        distance = sqrt(pow((X[0] - chronometr_vector[0]), 2.0L) +
                pow((X[1] - chronometr_vector[1]),2.0L) +
                pow((X[2] - chronometr_vector[2]),2.0L)) + (2.0L+0.01L*t)*1.1L + _generator.white_noise_generator(0.0L, 10.0L);
    }
    else
    {
        distance = 0.0L;
    }
    return distance;
}

