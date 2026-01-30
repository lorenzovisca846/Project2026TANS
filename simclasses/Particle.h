#ifndef PARTICLE_H
#define PARTICLE_H

#include <TObject.h>
#include "SimRandom.h"
#include "Point.h"

class Particle : public TObject
{
    public:
        Particle() {}
        Particle(SimRandom* srnd): fSimrand(srnd) {}

        void Init(const double& x, const double& y, const double& z, const double& beta, const double& p);
        double GetZ() const {return fZ;}

        void Propagation(const double& Rext);

        Point GetPoint() const {return Point(fX, fY, fZ);}

        void MultScatter(const double& X0, const double& W, const double& R);
    
    private:
        SimRandom* fSimrand;
        double tParam(const double& Rext);
        void AngleUpdate();

        double fX, fY, fZ;
        double fTheta, fPhi;
        double fP, fBeta;
        double fC1, fC2, fC3;
    
    ClassDef(Particle,1)
};

#endif // PARTICLE_H