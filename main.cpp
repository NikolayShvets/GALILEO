#include <iostream>
//#include "LA.h"
#include "tvector.h"
#include "tmatrix.h"
#include "integrator.h"
#include "custom.h"
#include <fstream>
#include "stdlib.h"
#include "iomanip"
#include "math.h"

using namespace std;

void calcModel(TModel* model, int i);
TMatrix get_transition_matrix(long double latitude, long double longtitude, bool flag);

int main()
{
    int satelliteNum = 27;
    //создаем 27 объектов "спутник", каждый из которых является моделью
    TModel* model[satelliteNum];
    for(int i = 0; i < satelliteNum; ++i)
    {
        if ( (i >= 0) && (i < 9) )
        {
            model[i] = new TSatellite(i, 0.785398L*i, 0.0L, 0.977384L );
        }
        if ( (i >= 9) && (i < 18) )
        {
            model[i] = new TSatellite(i, 0.785398L*i, 2.0944L, 0.977384L );
        }
        if ( (i >= 18) && (i < 27) )
        {
            model[i] = new TSatellite(i, 0.785398L*i, 2.0944L*2.0L, 0.977384L );
        }
    }
    TIntegrator* Integrator = new TDormandPrinceIntegrator();
    Integrator->setPrecision(1E-16);
    //поочередно для каждого спутника вызываем интегратор
    for(int i = 0 ; i < satelliteNum; ++i)
    {
        Integrator->Run(model[i]);
        calcModel(model[i], i);
    }
    for(int i = 0 ; i < satelliteNum; ++i)
    {
         //delete[] model[i];
    }
    //delete[] model;
    Integrator->~TIntegrator();
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
                //cout<<Result(i, j) << " ";
            }
            file << std::endl;
            //cout<<std::endl;
        }

        file.close();
}
TMatrix get_transition_matrix(long double latitude, long double longtitude, bool flag)
{
    if (flag)
    {
        TMatrix R(3,3);
        R[0][0] = - sin(latitude) * cos(longtitude);
        R[0][1] = - sin(latitude) * sin(longtitude);
        R[0][2] =   cos(latitude);
        R[1][0] =   cos(latitude) * cos(longtitude);
        R[1][1] =   cos(latitude) * sin(longtitude);
        R[1][2] = - sin(latitude);
        R[2][0] = - sin(longtitude);
        R[2][1] =   cos(longtitude);
        R[2][2] =   0;
        return R;
    }
    if(!flag)
    {
        TMatrix R(3,3);
        R[0][0] = -sin(longtitude);
        R[0][1] = cos(longtitude);
        R[0][2] = 0;
        R[1][0] = -sin(latitude)*cos(longtitude);
        R[1][1] = -sin(latitude)*sin(longtitude);
        R[1][2] = cos(latitude);
        R[2][0] = cos(latitude)*cos(longtitude);
        R[2][1] = cos(latitude)*sin(longtitude);
        R[2][2] = sin(latitude);
        return R;
    }
}
