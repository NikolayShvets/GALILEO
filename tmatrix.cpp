#include "tmatrix.h"

TMatrix::TMatrix()
{
}

TMatrix::TMatrix(int n, int m)
{
    resize(n, m);
}

TMatrix::TMatrix(const TMatrix &rval)
{
    this->data.assign(rval.data.begin(), rval.data.end());
}

void TMatrix::show_matrix()
{
    for(int i = 0; i < this->row_count(); ++i)
    {
        for(int j = 0; j < this->col_count(); ++j)
        {
            std::cout<<this->data[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }
}

void TMatrix::resize(int n, int m)
{
    this->n = n; this->m = m;
    data.resize(n);
    for (int i = 0; i < n; i++)
        data[i].resize(m);
}

TMatrix& TMatrix::operator = (const TMatrix& rval)
{
    this->data.assign(rval.data.begin(), rval.data.end());
    return (*this);
}

TMatrix TMatrix::operator - () const
{
    TMatrix temp(this->row_count(), this->col_count());
    for (int i = 0; i < temp.row_count(); i++)
        for (int j = 0; j < temp.col_count(); j++)
            temp[i][j] = -this->data[i][j];
    return temp;
}

TMatrix TMatrix::operator - (const TMatrix& rval) const
{
    TMatrix temp(this->row_count(), this->col_count());
    for (int i = 0; i < this->row_count(); i++)
        for (int j = 0; j < this->col_count(); j++)
            temp[i][j] = this->data[i][j] - rval.data[i][j];
    return temp;
}
TMatrix TMatrix::operator + (const TMatrix& rval) const
{
    TMatrix temp(this->row_count(), this->col_count());
    for (int i = 0; i < this->row_count(); i++)
        for (int j = 0; j < this->col_count(); j++)
            temp[i][j] = this->data[i][j] + rval[i][j];
    return temp;
}

TMatrix TMatrix::operator * (long double rval) const
{
    TMatrix temp(this->row_count(), this->col_count());
    for (int i = 0; i < this->row_count(); i++)
        for (int j = 0; j < this->col_count(); j++)
            temp[i][j] = this->data[i][j] * rval;
    return temp;
}

TMatrix TMatrix::operator * (const TMatrix& rval) const
{
    TMatrix temp(this->row_count(), rval.col_count());
    for (int i = 0; i < temp.row_count(); i++)
        for (int j = 0; j < temp.col_count(); j++)
        {
            temp[i][j] = 0.0L;
            for (int k = 0; k < this->row_count(); k++)
                temp[i][j] += this->data[i][j] * rval[i][j];
        }
    return temp;
}

TVector TMatrix::operator * (const TVector& rval) const
{
    TVector temp (this->row_count());
    for (int i = 0; i < row_count(); i++)
    {
        temp[i] = 0.0L;
        for (int k = 0; k < col_count(); k++)
            temp[i] += this->data[i][k] * rval[k];
    }
    return temp;
}

TMatrix::~TMatrix()
{
    this->data.erase(this->data.begin(), this->data.end());
}
void TMatrix::clear()
{
    this->data.clear();
}
