#include "tvector.h"

TVector::TVector() : n{0} {}
TVector::TVector(int n) : n{0}
{
    resize(n);
}
TVector::TVector(const TVector& rval) : n{0}
{
    this->data.assign(rval.data.begin(), rval.data.end());
}
TVector& TVector::operator = (const TVector& rval)
{
    this->data = rval.data;
    return *this;
}
TVector::~TVector()
{
    this->data.erase(this->data.begin(), this->data.end());
}
void TVector::resize(int n)
{
    this->data.resize(n);
}
long double TVector::length() const
{
    long double length{0.0L};
    for (int i = 0; i < this->size(); ++i)
        length += pow(this->data[i], 2.0L);
    return sqrt(length);
}
TVector TVector::operator - () const
{
    TVector temp_vector(this->size());
    for (unsigned int i = 0; i < this->data.size(); ++i)
        temp_vector.data[i] = -this->data[i];
    return temp_vector;
}

TVector TVector::operator - (const TVector& rval) const
{
    TVector temp_vector(this->size());
    for (unsigned int i = 0; i < this->data.size(); ++i)
        temp_vector[i] = data[i] - rval[i];
    return temp_vector;
}

TVector TVector::operator - (const long double rval) const
{
    TVector temp_vector(this->size());
    for (int i = 0; i < this->size(); i++)
        temp_vector[i] = this->data[i] - rval;
    return temp_vector;
}

TVector TVector::operator + (const TVector& rval) const
{
    TVector temp_vector(this->size());
    for (unsigned int i = 0; i < this->data.size(); ++i)
        temp_vector[i] = this->data[i] + rval[i];
    return temp_vector;
}

TVector TVector::operator + (const long double rval) const
{
    TVector temp_vector(this->size());
    for (int i = 0; i < this->size(); i++)
        temp_vector[i] = this->data[i] + rval;
    return temp_vector;
}

TVector TVector::operator * (const long double rval) const
{
    TVector temp_vector(this->size());
    for (int i = 0; i < this->size(); i++)
        temp_vector[i] = this->data[i] * rval;
    return temp_vector;
}

long double TVector::operator * (const TVector& rval) const
{
    long double temp_double = 0.0L;
    for (int i = 0; i < this->size(); i++)
        temp_double += this->data[i] * rval[i];
    return temp_double;
}

TVector TVector::operator ^ (const TVector& rval) const
{
    TVector temp_vector(3);
    temp_vector[0] = this->data[1]*rval[2] - this->data[2]*rval[1];
    temp_vector[1] = -(this->data[0]*rval[2] - this->data[2]*rval[0]);
    temp_vector[2] = this->data[0]*rval[1] - this->data[1]*rval[0];
    return temp_vector;
}
