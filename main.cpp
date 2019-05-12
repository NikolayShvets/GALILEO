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
    std::ofstream distance_file("distance_file.txt");
    std::ofstream true_distance_file("true_distances.txt");
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

    //этап 1 моделирование
    //поочередно для каждого спутника вызываем интегратор

    for(int i = 0 ; i < satelliteNum; ++i)
    {
        Integrator->Run(model[i]);
        finish_modeling.push_back(model[i]->getResult());
        calcModel(model[i], i);
    }
    Integrator->~TIntegrator();

    //этап 2 подготовка к МНК

    chronometr->navigation(finish_modeling, false); //дальности по результатам моделирования
    //запомнить моменты времени, номера спутников и порядок
    distance_file<<std::fixed;
    for(int i = 0; i < chronometr->true_distances.row_count(); ++i)
    {
        for(int j = 0; j < chronometr->true_distances.col_count(); ++j)
        {
            distance_file<<chronometr->true_distances[i][j]<<"\t";
        }
        distance_file<<std::endl;
    }
    distance_file.close();
    std::cout<<chronometr->true_distances.row_count()<<"X"<<chronometr->true_distances.col_count()<<std::endl;
    chronometr->set_measure_number(0);
    chronometr->true_distances.clear();

    //этап 3 начало МНК
    //смещаем начальные условия

    chronometr->set_lat(chronometr->get_lat()+chronometr->get_delta());
    chronometr->set_lon(chronometr->get_lon()+chronometr->get_delta());
    chronometr->set_betta_0(chronometr->get_betta_0()+chronometr->get_delta());
    chronometr->set_betta_1(chronometr->get_betta_1()+chronometr->get_delta());
    std::cout<<fixed;
    std::cout<<"latitude: "<<chronometr->get_lat()<<"\n"<<"longtitude: "<<chronometr->get_lon()<<std::endl;
    //в том же порядке, те же номера спутников, не смотря на видимость, в те же моменты времени считаем дальности
    chronometr->navigation(finish_modeling, false);
    std::cout<<chronometr->true_distances.row_count()<<"X"<<chronometr->true_distances.col_count()<<std::endl;
    true_distance_file<<std::fixed;
    for(int i = 0; i < chronometr->true_distances.row_count(); ++i)
    {
        for(int j = 0; j < chronometr->true_distances.col_count(); ++j)
        {
            true_distance_file << chronometr->true_distances[i][j] << "\t";
        }
        true_distance_file<<std::endl;
    }
    true_distance_file.close();
    chronometr->set_measure_number(0);

    //этап 4 матрица H
    //отклоняем на eps н.у из (2)
    chronometr->set_lat(chronometr->get_lat()-chronometr->get_delta() + chronometr->get_eps());
    chronometr->set_lon(chronometr->get_lon()-chronometr->get_delta() + chronometr->get_eps());
    chronometr->set_betta_0(chronometr->get_betta_0()-chronometr->get_delta() + chronometr->get_eps());
    chronometr->set_betta_1(chronometr->get_betta_1()-chronometr->get_delta() + chronometr->get_eps());
    chronometr->navigation(finish_modeling, true);
    //chronometr->derivatives.show_matrix();
    //std::cout<<chronometr->derivatives.row_count()<<"X"<<chronometr->derivatives.col_count()<<std::endl;
    //этап 5 матрица D
    //по главной диагонали дисперсии этта
    std::cout<<chronometr->D.row_count()<<"X"<<chronometr->D.col_count()<<std::endl;
    for(int k = 0; k < chronometr->derivatives.size(); ++k)
    { //std::cout<<"number: "<<k<<" "<<chronometr->derivatives[k].row_count()<<"X"<<chronometr->derivatives[k].col_count()<<std::endl;
        chronometr->derivatives[k].show_matrix();
        std::cout<<std::endl;
    }

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

