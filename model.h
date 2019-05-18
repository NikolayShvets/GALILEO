#pragma once
/*#include "tvector.h"
#include "tmatrix.h"*/
#include <armadillo>
using namespace std;
using namespace arma;
class TModel
{
    protected:
        vec X0;
        long double SamplingIncrement;
        long double t0, t1;
        mat Result;
		int N;
    public:
        TModel() : SamplingIncrement(1.), t0( 0. ), t1( 10000. ), N( 0. ) {}

        virtual void getRight( const vec& X, long double t, vec& Y ) = 0;
		
        inline vec getInitialConditions() const { return X0; }
        inline int getOrder() const { return X0.size(); }

        inline long double getSamplingIncrement() const { return SamplingIncrement; }

        inline long double getT0() const { return t0; }
        inline long double getT1() const { return t1; }
       
        inline mat getResult() const { return Result; }
        virtual void addResult( const vec& X, long double t );
		virtual void clearResult();
		virtual void prepareResult();
};

