#include "SimRandom.h"

ClassImp(SimRandom)

SimRandom::SimRandom() : TRandom3(0),
    fSeed(0)
{

}

SimRandom::SimRandom(unsigned int seed) : TRandom3(seed),
    fSeed(seed)
{

}

SimRandom::~SimRandom()
{

}

int SimRandom::VMult3(int min, int max)
{
    //Implement sampling from custom distribution (histogram)
    return 0;
}

Point SimRandom::GausPoint(double xyS, double zS)
{
    double x = Gaus(0., xyS);
    double y = Gaus(0., xyS);
    double z = Gaus(0., zS);
    return Point(x, y, z);
}
