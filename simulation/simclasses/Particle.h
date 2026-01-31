#ifndef PARTICLE_H
#define PARTICLE_H

#include <TObject.h>
#include "SimRandom.h"
#include "MyPoint.h"

class Particle : public TObject
{
    public:
        Particle() {}
        Particle(SimRandom* srnd): fSimrand(srnd) {}

        void Init(double x, double y, double z, double beta, double p);
        void Propagation(double Rext);
        void MultScatter(double X0, double W, double R);

        double GetZ() const {return fZ;}
        MyPoint GetPoint() const {return MyPoint(fX, fY, fZ);}
    
    private:
        SimRandom* fSimrand;
        void AngleUpdate();

        double fX, fY, fZ;                      // coordinates
        double fTheta, fPhi;                    // angles
        double fP, fBeta;                       // momentum and beta
        double fC1, fC2, fC3;                   // direction cosines
    
    ClassDef(Particle,1)
};

#endif // PARTICLE_H