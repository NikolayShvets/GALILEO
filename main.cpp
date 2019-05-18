#include <iostream>
//#include
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
using namespace arma;

void calcModel(TModel* model, /*Consumer* chronometr,*/ int i);

int main()
{
    std::ofstream init_distances_file("init_dist.txt");
    std::ofstream true_distances_file("true_dist.txt");
    init_distances_file<<std::fixed;
    true_distances_file<<std::fixed;
    std::ofstream mnk_file("mnk_file.txt");
    std::ofstream deltas("deltas_file.txt");
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
    vector<mat> finish_modeling;
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
    //cout<<finish_modeling[3](18500,0)<<endl;

    //далее МНК

    //единожды получаем истинные результаты c ошибкой на каждый момент времени
    //--> вектор истиных дальностей
    //с шагом 1с реальные дальности всех спутников с ошибкой, с проверкой на видимость
    chronometr->navigation(finish_modeling, true, false, false);
    chronometr->distances.clear();
    for (int i = 0; i <chronometr->init_distances.n_rows; ++i)
    {
        init_distances_file<<std::fixed;
        init_distances_file<<chronometr->init_distances(i,2)<<std::endl;
    }
   // std::cout<<chronometr->init_distances<<endl;

    //инициализируем вектор поправок delta_x
    chronometr->delta_x.resize(5);
    for(int i = 0; i < chronometr->delta_x.n_rows; ++i)
        chronometr->delta_x[i] = 1.0L;
    //инициализируем вектор ну(оцениваемые параметры) на 0 шаг мнк (забиваем нулями)
    chronometr->measure_vector.resize(5);
    for(int i = 0; i < chronometr->measure_vector.n_rows; ++i)
        chronometr->measure_vector[i] = 0.0L;

    //счетчик итераций МНК
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
        /*if(counter == 10)
            break;*/
        std::cout<<"Iteration "<<counter<<": "<<std::endl;
        chronometr->set_measure_number(0);
        //на каждый момент времени для каждого спутника в области видимости
        //для каждой истинной дальности в этот момент времени, рассчитываем
        //расчетную дальность без ошибки
        //с шагом 150с расчетные дальности c ошибкой всех спутников в области видимости

        /*chronometr->measure_vector[0] += chronometr->get_delta(); //lat
        chronometr->measure_vector[1] += chronometr->get_delta(); //lon
        chronometr->measure_vector[2] += chronometr->get_delta(); //H
        chronometr->measure_vector[3] += chronometr->get_delta(); //b0
        chronometr->measure_vector[4] += chronometr->get_delta(); //b1*/

        chronometr->navigation(finish_modeling, false, true, true);
        chronometr->distances.clear();

        //формируем разности между реальными дальностями и расчетными по моментам времени расчетных дальностей
        chronometr->delta_y.resize(chronometr->get_measure_number());
        int c {0};
        for(int i = 0; i < chronometr->true_distances.n_rows; ++i)
        {
            for(int j = 0; j < chronometr->init_distances.n_rows; ++j)
            {
                if(
                    (chronometr->true_distances(i,0) == chronometr->init_distances(j,0)) &&
                    (chronometr->true_distances(i,1) == chronometr->init_distances(j,1))
                  )
                {
                    chronometr->delta_y[c] = chronometr->init_distances(j,2) - chronometr->true_distances(i,2);
                    ++c;
                }
            }
        }

        // D - диагональная матрица дисперсий ошибок этта
        chronometr->D.resize(chronometr->get_measure_number(),chronometr->get_measure_number());
        chronometr->D = chronometr->D.eye();
        chronometr->D = chronometr->D * powl(chronometr->get_d(),2.0L);
        chronometr->set_measure_number(0);
        //вычисляем поправки к вектору начальных условий (p,l,h,b0,b1)T; delta_x = (Ht*!D*H)^(-1)*Ht*D^(-1)*delta_y
       // TMatrix temp_test(5,5),test(5,5), _test(5,5); //вспомогательная матрица
        //(Ht*!D*H)^(-1)
        std::cout<<std::fixed;
      //  std::cout<<"D: \n"<<chronometr->D<<std::endl;
        std::cout<<chronometr->D.n_rows<<"x"<<chronometr->D.n_cols<<std::endl;

        std::cout<<std::fixed;
        std::cout<<"H: \n"<<chronometr->derivatives<<std::endl;
      //  std::cout<<chronometr->derivatives.n_rows<<"x"<<chronometr->derivatives.n_cols<<std::endl;
        std::cout<<std::fixed;

        std::cout<<"delta_y: \n"<<chronometr->delta_y<<std::endl;
      //  std::cout<<chronometr->delta_y.n_rows<<"x"<<"1"<<std::endl;

        chronometr->delta_x = (chronometr->derivatives.t()*chronometr->D.i()*chronometr->derivatives).i()*chronometr->derivatives.t()*chronometr->D.i()*chronometr->delta_y;

        std::cout<<std::fixed;
        std::cout<<"delta_x: "<<std::endl;
        std::cout<<chronometr->delta_x<<std::endl;
        std::cout<<chronometr->delta_x.n_rows<<"x"<<"1"<<std::endl;

        chronometr->measure_vector = chronometr->measure_vector + chronometr->delta_x;
        std::cout<<"measure_vector: "<<std::endl;
        std::cout<<chronometr->measure_vector<<std::endl;

        mnk_file<<std::fixed;
        for(int i = 0; i < chronometr->measure_vector.n_rows; ++i)
            mnk_file<<chronometr->measure_vector[i]<<" ";
        mnk_file<<std::endl;

        deltas<<std::fixed;
        for(int i = 0; i < chronometr->delta_x.n_rows; ++i)
                deltas<<chronometr->delta_x[i]<<" ";
        deltas<<std::endl;
        mnk_file.flush();
        deltas.flush();
    }

    mnk_file.close();
    init_distances_file.close();
    true_distances_file.close();
    return 0;
}

void calcModel(TModel* model, int i)
{
    std::ofstream file("Integration_results_" + std::to_string(i)+ ".txt");
        mat Result = model->getResult();

        for (int i=0; i<Result.n_rows; i++)
        {
            for (int j=0; j<Result.n_cols; j++)
            {
                file<<fixed;
                file << Result(i,j) << " ";
                cout<<fixed;
            }
            file << std::endl;
        }

        file.close();
}

