#ifndef MNK_H
#define MNK_H
#include "tvector.h"
#include "tmatrix.h"


class MNK
{
public:
    MNK();
    TMatrix derivates(TVector parametrs);
};

#endif // MNK_H
