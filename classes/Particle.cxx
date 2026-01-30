#include "Particle.h"

ClassImp(Particle)

void Particle::Init(double x, double y, double z, double beta, double p, int q)
    {
        fX = x;
        fY = y;
        fZ = z;
        fP = p;
        fBeta = beta;
        fQ = q;

        fPhi = fSimrand->PhiDist();
        fTheta = fSimrand->ThetaDist();
    }