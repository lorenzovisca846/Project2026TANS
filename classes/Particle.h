#ifndef PARTICLE_H
#define PARTICLE_H

#include <TObject.h>
#include "SimRandom.h"

class Particle : public TObject
{
    public:
        Particle() {}
        Particle(SimRandom* srnd): fSimrand(srnd) {}

        void Init(double x, double y, double z, double beta, double p, int trackID);
        void Propagation(double Rext);
        void MultScatter(double X0, double W, double R);

        double GetZ() const {return fZ;}
        double GetR() const {return sqrt(fX*fX + fY*fY);}
        double GetPhi() const {return fPhi;}
        int GetTrackID() const {return fTrackID;}     
    
    private:
        SimRandom* fSimrand;
        void AngleUpdate();

        double fX, fY, fZ;                      // coordinates
        double fTheta, fPhi;                    // angles
        double fP, fBeta;                       // momentum and beta
        double fC1, fC2, fC3;                   // direction cosines

        int fTrackID;                           // track ID
    
    ClassDef(Particle,1)
};

#endif // PARTICLE_H