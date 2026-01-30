#ifndef PARTICLE_H
#define PARTICLE_H

#include <TObject.h>
#include "SimRandom.h"

class Particle : public TObject
{
    public:
        Particle() {}
        Particle(SimRandom* srnd): fSimrand(srnd) {}

        void Init(double x, double y, double z, double beta, double p, int q);


    
    private:
        SimRandom* fSimrand;

        double fX;
        double fY;
        double fZ;
        double fTheta;
        double fPhi;
        double fP;
        double fBeta;
        int fQ;
    
    ClassDef(Particle,1)
};

#endif