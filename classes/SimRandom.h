#ifndef SIMRANDOM_H
#define SIMRANDOM_H

#include <TRandom3.h>
#include <TH1.h>
#include "Point.h"

class SimRandom : public TRandom3 {
    public:
        SimRandom();
        SimRandom(unsigned int seed);

        ~SimRandom();

        int VMult1(int min, int max) const {return max;}
        int VMult2(int min, int max) {return min + Integer(max - min + 1);}
        int VMult3(int min, int max);

        Point GausPoint(double, double); 



    private:
    double fSeed;

    ClassDef(SimRandom,1);
};

#endif // SIMRANDOM_H