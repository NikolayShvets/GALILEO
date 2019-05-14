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
    long double latitude{0.0L}, longtitude{0.0L}, betta_0{0.0L}, betta_1{0.0L}, measurement_period{0.0L}, m{0.0L}, d{0.0L};
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
    std::cout<<"Input statistic: "<<std::endl;
        std::cout<< "M: ";
        std::cin >> m;
        std::cout<< "D: ";
        std::cin >> d;
    //объект потребителя
    Consumer *chronometr = new Consumer(latitude, longtitude, measurement_period, betta_0, betta_1, m, d);
    //количество спутников
    int satelliteNum = 4;
    //блочная матрица результатов интегрирования каждого спутника
    vector<TMatrix> finish_modeling;
    //создаем 27 объектов "спутник", каждый из которых является моделью
    TModel* model[satelliteNum];
    //инициализируем каждый спутник
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
    //поочередно для каждого спутника вызываем интегратор --> истиная траектория движения
    //каждого спутника в каждый момент времени
    for(int i = 0 ; i < satelliteNum; ++i)
    {
        Integrator->Run(model[i]);
        finish_modeling.push_back(model[i]->getResult());
        calcModel(model[i], i);
    }
    Integrator->~TIntegrator();

    //далее МНК

    //единожды получаем истинные результаты без ошибкии на каждый момент времени
    //--> вектор истиных дальностей
    //с шагом 1с реальные дальности всех спутников не смотря на видимость
    chronometr->navigation(finish_modeling, true, false, false);
    chronometr->distances.vector_clear();
    chronometr->set_measure_number(0);

    //инициализируем вектор поправок delta_x
    chronometr->delta_x.resize(5);
    for(int i = 0; i < chronometr->delta_x.row_count(); ++i)
        chronometr->delta_x[i] = 1.0L;
    int counter{0};
    //пока каждая из координат вектора оцениваемых параметров не меньше чем епсилон пробегаем МНК
    while(
          (  abs(chronometr->delta_x[0]) > chronometr->get_eps())
          ||(abs(chronometr->delta_x[1]) > chronometr->get_eps())
          ||(abs(chronometr->delta_x[2]) > chronometr->get_eps())
          ||(abs(chronometr->delta_x[3]) > chronometr->get_eps())
          ||(abs(chronometr->delta_x[4]) > chronometr->get_eps())
          )
    {
        ++counter;
        if(counter == 10000)
            break;
        std::cout<<"Iteration "<<counter<<": "<<std::endl;
        chronometr->set_measure_number(0);
        //на каждый момент времени для каждого спутника в области видимости
        //для каждой истинной дальности в этот момент времени, рассчитываем
        //расчетную дальность (с ошибкой)
        //с шагом 150с расчетные дальности всех спутников в области видимости
        chronometr->measure_vector[0] += chronometr->get_delta(); //lat
        chronometr->measure_vector[1] += chronometr->get_delta(); //lon
        chronometr->measure_vector[2] += chronometr->get_delta(); //H
        chronometr->measure_vector[3] += chronometr->get_delta(); //b0
        chronometr->measure_vector[4] += chronometr->get_delta(); //b1

        chronometr->navigation(finish_modeling, false, true, true);
        chronometr->distances.vector_clear();


        //формируем разности между реальными дальностями и расчетными по моментам времени расчетных дальностей
        chronometr->delta_y.resize(chronometr->get_measure_number());
        int counter {0};
        for(int i = 0; i < chronometr->true_distances.row_count(); ++i)
        {
            for(int j = 0; j < chronometr->init_distances.row_count(); ++j)
            {
                if(
                    (chronometr->true_distances[i][0] == chronometr->init_distances[j][0]) &&
                    (chronometr->true_distances[i][1] == chronometr->init_distances[j][1])
                  )
                {
                    chronometr->delta_y[counter] = chronometr->init_distances[j][2] - chronometr->true_distances[i][2];
                    ++counter;
                }
            }
        }

        // D - диагональная матрица дисперсий ошибок этта
        chronometr->D.resize(chronometr->get_measure_number(),chronometr->get_measure_number());
        chronometr->D = chronometr->D.E(chronometr->get_measure_number());
        chronometr->D = chronometr->D * powl(chronometr->get_d(),2.0L);
        std::cout<<chronometr->derivatives.row_count()<<"x"<<chronometr->derivatives.col_count()<<std::endl;
        std::cout<<chronometr->D.row_count()<<"x"<<chronometr->D.col_count()<<std::endl;
        chronometr->set_measure_number(0);
        std::cout<<chronometr->delta_y.row_count()<<"x"<<"1"<<std::endl;
        //вычисляем поправки к вектору начальных условий (p,l,h,b0,b1)T; delta_x = (Ht*!D*H)^(-1)*Ht*D^(-1)*delta_y
        TMatrix temp_test(5,5); //вспомогательная матрица
        //(Ht*!D*H)^(-1)
        temp_test = (chronometr->derivatives.t()*chronometr->D.operator !()*chronometr->derivatives).operator !();
        std::cout<<"temp_test"<<std::endl;
        std::cout<<temp_test.row_count()<<"x"<<temp_test.col_count()<<std::endl;
        //Ht*D^(-1)*delta_y
        chronometr->delta_x = chronometr->derivatives.t()*chronometr->D.operator !() * chronometr->delta_y;
        std::cout<<chronometr->delta_x.row_count()<<"x"<<"1"<<std::endl;
        //одно умножаем на другое
        chronometr->delta_x = temp_test*chronometr->delta_x;
        std::cout<<chronometr->delta_x.row_count()<<"x"<<"1"<<std::endl;
        chronometr->measure_vector = chronometr->measure_vector + chronometr->delta_x;
        std::cout<<"delta_x: "<<std::endl;
        chronometr->delta_x.show_vector();
        std::cout<<std::endl;
        std::cout<<"measure_vector: "<<std::endl;
        chronometr->measure_vector.show_vector();
/*
        //смещаем начальные условия на дельту

        chronometr->measure_vector[0] += chronometr->get_delta(); //lat
        chronometr->measure_vector[1] += chronometr->get_delta(); //lon
        chronometr->measure_vector[2] += chronometr->get_delta(); //H
        chronometr->measure_vector[3] += chronometr->get_delta(); //b0
        chronometr->measure_vector[4] += chronometr->get_delta(); //b1

        //в том же порядке, те же номера спутников, не смотря на видимость, в те же моменты времени считаем дальности
        chronometr->navigation(finish_modeling, false, false, false);
        chronometr->delta_true_distances = chronometr->distances;
        chronometr->distances.vector_clear();
        chronometr->set_measure_number(0);

        //матрица H, отклоняем на +eps смещенные на дельту н.у
        chronometr->measure_vector[0] += chronometr->get_eps(); //lat
        chronometr->measure_vector[1] += chronometr->get_eps(); //lon
        chronometr->measure_vector[2] += chronometr->get_eps(); //H
        chronometr->measure_vector[3] += chronometr->get_eps(); //b0
        chronometr->measure_vector[4] += chronometr->get_eps(); //b1
        chronometr->navigation(finish_modeling, false, false, true);
        chronometr->H_plus_eps = chronometr->distances;
        chronometr->distances.vector_clear();
        chronometr->set_measure_number(0);

        //матрица H, отклоняем на -eps смещенные на дельту н.у
        chronometr->measure_vector[0] -= 2*chronometr->get_eps(); //lat
        chronometr->measure_vector[1] -= 2*chronometr->get_eps(); //lon
        chronometr->measure_vector[2] -= 2*chronometr->get_eps(); //H
        chronometr->measure_vector[3] -= 2*chronometr->get_eps(); //b0
        chronometr->measure_vector[4] -= 2*chronometr->get_eps(); //b1
        chronometr->navigation(finish_modeling, false, false, true);
        chronometr->H_minus_eps = chronometr->distances;
        chronometr->distances.vector_clear();

        //матрица D есть квадратная матрица к-во_измер Х к-во_измер с дисперсиями этта на г.д.
        chronometr->D.resize(chronometr->get_measure_number(),chronometr->get_measure_number());
        chronometr->D = chronometr->D.E(chronometr->get_measure_number());
        chronometr->D = chronometr->D * chronometr->get_d();
        chronometr->set_measure_number(0);*/


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
}

