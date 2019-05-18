//---------------------------------------------------------------------------

#include "model.h"
#include <iostream>

//---------------------------------------------------------------------------
void TModel::addResult( const vec& X, long double t )
{

    if (N == Result.n_rows)
        Result.resize(N + 1, getOrder() + 1);

    Result(N,0) = t;
    for (int i = X.size(); i > 0 ; i--)
    {
        Result(N,i) = X[i-1];
    }

	N++;
}

void TModel::clearResult() 
{
	Result.resize(0, 0); 
	N = 0; 
}


void TModel::prepareResult() 
{ 
	Result.resize((int)((t1 - t0)/SamplingIncrement) + 1, getOrder() + 1);

	N = 0; 
}
