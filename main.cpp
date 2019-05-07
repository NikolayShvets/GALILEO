#include <iostream>
//#include "LA.h"
#include "tvector.h"
#include "tmatrix.h"
#include "integrator.h"
#include "custom.h"
#include "consumer.h"
#include <fstream>
#include "stdlib.h"
#include "iomanip"
#include "math.h"

using namespace std;

void calcModel(TModel* model, /*Consumer* chronometr,*/ int i);

int main()
{
    long double latitude{0.0L}, longtitude{0.0L};
    long double betta_0{0.0L}, betta_1{0.0L}, measurement_period{0.0L};
    std::cout<<"Input Consumer coordinates: "<<std::endl;
        std::cout<< "Consumer latitude: ";
        std::cin >> latitude;
        std::cout<< "Consumer longtitude: ";
        std::cin >> longtitude;
    std::cout<<"Input time scale parametrs: "<<std::endl;
        std::cout<< "Betta_0: ";
        std::cin >> betta_0;
        std::cout<< "Betta_1: ";
        std::cin >> betta_1;
        std::cout<< "Measurement period: ";
        std::cin >> measurement_period;
    Consumer *chronometr = new Consumer(latitude, longtitude, measurement_period, betta_0, betta_1);
    int satelliteNum = 27;
    vector<TMatrix> finish_modeling;

    //создаем 27 объектов "спутник", каждый из которых является моделью
    TModel* model[satelliteNum];
    for(int i = 0; i < satelliteNum; ++i)
    {
        if ( (i >= 0) && (i < 9) )
        {
            model[i] = new TSatellite(i, 0.785398L*i, 0.0L, 0.977384L);
        }
        if ( (i >= 9) && (i < 18) )
        {
            model[i] = new TSatellite(i, (0.785398L- 0.226893L)*(i-9), 2.0944L, 0.977384L);
        }
        if ( (i >= 18) && (i < 27) )
        {
            model[i] = new TSatellite(i, (0.785398L+ 0.226893L)*(i-18) , 2.0944L*2.0L, 0.977384L);
        }
    }
    TIntegrator* Integrator = new TDormandPrinceIntegrator();
    Integrator->setPrecision(1E-16);
    //поочередно для каждого спутника вызываем интегратор
    for(int i = 0 ; i < satelliteNum; ++i)
    {
        //этап 1 моделирование
        Integrator->Run(model[i]);
        finish_modeling.push_back(model[i]->getResult());
        calcModel(model[i], i);
    }
    Integrator->~TIntegrator();
    //этап 2 подготовка к МНК
    chronometr->navigation(finish_modeling, false); //дальности по результатам моделирования
    //этап 3 начало МНК
    chronometr->set_lat(chronometr->get_lat()+chronometr->get_delta());
    chronometr->set_lon(chronometr->get_lon()+chronometr->get_delta());
    chronometr->set_betta_0(chronometr->get_betta_0()+chronometr->get_delta());
    chronometr->set_betta_1(chronometr->get_betta_1()+chronometr->get_delta()*0.001L);
    chronometr->navigation(finish_modeling, true);
    return 0;
}

void calcModel(TModel* model, int i)
{
    std::ofstream file("Integration_results_" + std::to_string(i)+ ".txt");
        TMatrix Result = model->getResult();

        for (int i=0; i<Result.row_count(); i++)
        {
            for (int j=0; j<Result.col_count(); j++)
            {
                file<<fixed;
                file << Result[i][j] << " ";
                cout<<fixed;
            }
            file << std::endl;
        }

        file.close();
        //chronometr->navigation(model->getResult(), i);
}

