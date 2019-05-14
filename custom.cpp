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
    D1 = powl( powl( X[0] + m, 2 ) + powl( X[1], 2 ), 1.5 );
    D2 = powl( powl( X[0] + m - 1, 2 ) + powl( X[1], 2 ), 1.5 );
    Y[0] = X[2];
    Y[1] = X[3];
    Y[2] = X[0] + 2 * X[3] - (1 - m)*(X[0] + m) / D1 - m * (X[0] + m - 1) / D2;
    Y[3] = X[1] - 2 * X[2] - (1 - m)*X[1] / D1 - m * X[1] / D2;
}
//---------------------------------------------------------------------------

TSatellite::TSatellite(int satellite_number,long double d_u, long double d_omega, long double i) : TModel()
{
    cout<<"\nSatellite number "<<satellite_number<<": \n"<<endl;
    ch_lat = 0.95993108859;//M_PI/2.50L;//1.291L;//Москва
    ch_lon = 0.0L;//0.64577182324//2.523L;
    long double latitude{0.0L}, longtitude{0.0L}, r_vector{0.0L}, v_vector{0.0L};
    this->satellite_number = satellite_number;
    this->u += d_u;
    this->i = i;
    this->omega += d_omega;
    this->A.resize(3,3);
    this->A[0][0] = cosl(u)*cosl(omega) - sinl(u)*sinl(omega)*cosl(i); this->A[0][1] = -sinl(u)*cosl(omega) - cosl(u)*sinl(omega)*cosl(i); this->A[0][2] = sinl(i)*sinl(omega);
    this->A[1][0] = cosl(u)*sinl(omega) + sinl(u)*cosl(omega)*cosl(i); this->A[1][1] = -sinl(u)*sinl(omega) + cosl(u)*cosl(omega)*cosl(i); this->A[1][2] = -sinl(i)*cosl(omega);
    this->A[2][0] = sinl(u)*sinl(i);                                this->A[2][1] = cosl(u)*sinl(i) ;                                 this->A[2][2] = cosl(i);

    cout<<"Matrix A: "<<endl;
    cout<< "\t" << A[0][0] << "\t" << A[0][1] << "\t" << A[0][2] << endl;
    cout<< "\t" << A[1][0] << "\t" << A[1][1] << "\t" << A[1][2] << endl;
    cout<< "\t" << A[2][0] << "\t" << A[2][1] << "\t" << A[2][2] << endl;

    this->X_oscul.resize(3);
    X_oscul[0] = orb_height;
    X_oscul[1] = 0.0L;
    X_oscul[2] = 0.0L;

    this->V_oscul.resize(3);
    V_oscul[0] = 0.0L;//powl((nu/orb_height), 0.5L);
    V_oscul[1] = powl((nu/orb_height), 0.5L);
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

    latitude = atan(X0[2]/sqrtl( ( powl(X0[0],2.L) + powl(X0[1], 2.L) )));
    longtitude = atan(X0[1]/X0[0]);
    r_vector = sqrtl(powl(X0[0],2.) + powl(X0[1], 2.) + powl(X0[2], 2.));
    v_vector = 3677.233;//sqrtl(powl(temp_velosities[0],2.) + powl(temp_velosities[1], 2.) + powl(temp_velosities[2], 2.));

    cout<< "Sfericheskie: "<<endl;
    cout<< "    lat: "<<latitude<<endl;
    cout<< "    lon: "<<longtitude<<endl;
    cout<< "    r_v: "<<r_vector<<endl;
    cout<< "    v_v: "<<v_vector<<endl;

    X0[3] = /*v_vector*cosl(latitude)*sinl(longtitude);*/temp_velosities[0];
    X0[4] = /*-v_vector*cosl(latitude)*cosl(longtitude);*/temp_velosities[1];
    X0[5] = /*v_vector*sinl(latitude);*/temp_velosities[2];

    cout<<"Coordinates: "<<endl;
    cout<< "    x: "<<X0[0]<<" m"<<endl;
    cout<< "    y: "<<X0[1]<<" m"<<endl;
    cout<< "    z: "<<X0[2]<<" m"<<endl;
    cout<< "    |r-vector| = "<<sqrtl(powl(X0[0],2.) + powl(X0[1], 2.) + powl(X0[2], 2.))<<endl;

    cout<<"Velosities: "<<endl;
    cout<< "    Vx: "<<X0[3]<<" m/sec"<<endl;
    cout<< "    Vy: "<<X0[4]<<" m/sec"<<endl;
    cout<< "    Vz: "<<X0[5]<<" m/sec"<<endl;
    cout<< "    |v-vector| = "<<sqrtl(powl(X0[3],2.) + powl(X0[4], 2.) + powl(X0[5], 2.))<<endl;
}
//---------------------------------------------------------------------------

void TSatellite::getRight(const TVector &X, long double t, TVector &Y)
{
    Y.resize(6);
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    long double ro = sqrtl(powl(X[0],2.) + powl(X[1], 2.) + powl(X[2], 2.));
    Y[3] = -nu*X[0]/powl(ro,3.);
    Y[4] = -nu*X[1]/powl(ro,3.);
    Y[5] = -nu*X[2]/powl(ro,3.);
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
    chronometr_vector[0] = Re*cosl(ch_lat)*cosl(ch_lon);
    chronometr_vector[1] = Re*cosl(ch_lat)*sinl(ch_lon);
    chronometr_vector[2] = Re*sinl(ch_lat);
    long double alfa = acosl(temp_v*chronometr_vector/(temp_v.length()*chronometr_vector.length()));
    //cout<<"ugol: "<<alfa*180.0L/M_PI<<" lon chronom: "<<ch_lon*180.0L/M_PI<<" time: "<<t - t0<<endl;
    if (alfa < M_PI/3.0L) visibility = true; else visibility = false;
    if (visibility == true)
    {//рассчет расстояния
        //cout<<_generator.white_noise_generator()<<endl;
        distance = sqrtl(powl((X[0] - chronometr_vector[0]), 2.0L) +
                powl((X[1] - chronometr_vector[1]),2.0L) +
                powl((X[2] - chronometr_vector[2]),2.0L)) + (2.0L+0.01L*t)*1.1L + _generator.white_noise_generator(0.0L, 10.0L);
    }
    else
    {
        distance = 0.0L;
    }
    return distance;
}

