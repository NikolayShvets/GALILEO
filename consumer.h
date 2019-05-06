#ifndef CONSUMER_H
#define CONSUMER_H
#include "tvector.h"
#include "tmatrix.h"
#include <fstream>
#include "generator.h"

class Consumer
{
protected:
    std::vector<std::ofstream> navigation_log_vector;
    long double latitude{0.0L}, longtitude{0.0L}, step{1.0L};
    long double light_speed{299792458.0L}, betta_0{0.0L}, betta_1{0.0L};
    long double frequency{0.0L};
public:
    long double distance{0.0L};
    generator _generator;
    bool visibility{false};
    long double Re{6371000.0L}, omega{7.292115E-5};
    TVector consumer_vector;
    Consumer(long double latitude, long double longtitude,long double frequency, long double betta_0, long double betta_1);
    ~Consumer();
    void navigation(TMatrix result, int model_number);
};

#endif // CONSUMER_H
