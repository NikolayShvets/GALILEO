#ifndef CONSUMER_H
#define CONSUMER_H
#include "tvector.h"
#include "tmatrix.h"
#include <fstream>
#include "generator.h"

class Consumer
{
protected:
    vector<ofstream> navigation_log_vector;
    long double latitude{0.0L}, longtitude{0.0L}, step{1.0L};
    long double light_speed{299792458.0L/*0.0L*/}, betta_0{0.0L}, betta_1{0.0L};
    long double measurement_period{0.0L};
    int measure_number{0};
    long double eps{1e-10};
    long double delta{1e-6};
public:
    int current_s_num{0}, k{-1}, s{0};
    TMatrix D;
    vector<TMatrix> derivatives;
    TMatrix true_distances;
    long double distance{0.0L};
    generator _generator;
    bool visibility{false};
    long double Re{6371000.0L}, omega{7.292115E-5};
    TVector consumer_vector;
    Consumer(long double latitude, long double longtitude,long double measurement_period, long double betta_0, long double betta_1);
    ~Consumer();
    void navigation(std::vector<TMatrix> finish_modeling, bool flag);

    long double get_lat () {return latitude;}
    void set_lat(long double latitude) {this->latitude = latitude;}

    long double get_lon() {return longtitude;}
    void set_lon(long double longtitude) {this->longtitude = longtitude;}

    long double get_eps() {return eps;}
    long double get_delta() {return delta;}

    long double get_betta_0 () {return betta_0;}
    void set_betta_0(long double betta_0) {this->betta_0 = betta_0;}

    long double get_betta_1() {return betta_1;}
    void set_betta_1(long double betta_1) {this->betta_1 = betta_1;}

    long double get_measure_number() {return measure_number;}
    void set_measure_number(long double measure_number) {this->measure_number = measure_number;}

    void get_derivatives(TVector X_navigation_spacecraft, long double t, int measure_number, int satellite_number);
};

#endif // CONSUMER_H
