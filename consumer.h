#ifndef CONSUMER_H
#define CONSUMER_H
#include "tvector.h"
#include "tmatrix.h"
#include <armadillo>
#include <fstream>
#include "generator.h"
using namespace arma;

class Consumer
{
protected:
    vector<ofstream> navigation_log_vector;
    long double latitude{0.0L}, longtitude{0.0L}, step{1.0L};
    long double light_speed{299792458.0L/*1.0L*/}, betta_0{0.0L}, betta_1{0.0L};
    long double measurement_period{0.0L};
    int measure_number{0};
    long double eps{1e-4};
    long double delta{1e-10};
    long double m{0.0L}, d{0.0L};
public:
    vec measure_vector;
    int current_s_num{0}, k{-1}, s{0};
    vec delta_x;
    vec delta_y;
    mat true_distances; //один раз, рассчетные
    mat init_distances; //каждый шаг, реальные наперед известные
    vec delta_true_distances;
    vec H_plus_eps;
    vec H_minus_eps;
    mat D;
    mat derivatives;
    vec distances;
    long double distance{0.0L};
    generator _generator;
    bool visibility{false};
    long double Re{6371000.0L}, omega{7.292115E-5};
    vec consumer_vector;
    Consumer(long double latitude, long double longtitude,long double measurement_period, long double betta_0, long double betta_1, long double m, long double d);
    ~Consumer();
    void navigation(std::vector<mat> finish_modeling, bool w_err);

    long double get_lat () {return latitude;}
    void set_lat(long double latitude) {this->latitude = latitude;}

    long double get_m () {return m;}
    void set_m(long double m) {this->m = m;}


    long double get_d () {return d;}
    void set_d(long double d) {this->d = d;}

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

    void get_derivatives(vec X_navigation_spacecraft, vec consumer_AES, long double t, long double curr_lon, int measure_number);
};

#endif // CONSUMER_H
