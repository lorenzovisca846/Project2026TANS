#include "Particle.h"

ClassImp(Particle)

void Particle::Init(const double& x, const double& y, const double& z, const double& beta, const double& p)
    {
        fX = x;
        fY = y;
        fZ = z;
        fP = p;
        fBeta = beta;
        
        fPhi = fSimrand->PhiDist();
        fTheta = fSimrand->ThetaDist();

        fC1 = sin(fTheta) * cos(fPhi);
        fC2 = sin(fTheta) * sin(fPhi);
        fC3 = cos(fTheta);
    }

void Particle::Propagation(const double& Rext)
{
    double t = tParam(Rext);
    
    fX += t * fC1;
    fY += t * fC2;
    fZ += t * fC3;
}

double Particle::tParam(const double& Rext)
{
    double delta = (fX*fC1 + fY*fC2)*(fX*fC1 + fY*fC2) - (fC1*fC1 + fC2*fC2)*(fX*fX + fY*fY - Rext*Rext);
    return  (sqrt(delta) - (fX*fC1 + fY*fC2)) / (fC1*fC1 + fC2*fC2);
}

void Particle::MultScatter(const double& X0, const double& W, const double& R)
{
    double thickness = W/X0; //Correggere con angolo di incidenza per spessore effettivo attraversato

    double dTheta = fSimrand->ScatterDist(fP, fBeta, thickness); //Controlla se va messo sqrt2
    double dPhi = fSimrand->PhiDist();

    double ms[3][3];
    ms[0][0] = -sin(fPhi);
    ms[0][1] = -cos(fTheta)*cos(fPhi);
    ms[0][2] = sin(fTheta)*cos(fPhi);
    ms[1][0] = cos(fPhi);
    ms[1][1] = -cos(fTheta)*sin(fPhi);
    ms[1][2] = sin(fTheta)*sin(fPhi);
    ms[2][0] = 0.;
    ms[2][1] = sin(fTheta);
    ms[2][2] = cos(fTheta);

    double dC[3]={sin(dTheta)*cos(dPhi), sin(dTheta)*sin(dPhi), cos(dTheta)};

    fC1 = 0;
    fC2 = 0;
    fC3 = 0;

    for(int i=0;i<3;i++)
    {
        fC1 += ms[0][i]*dC[i];
        fC2 += ms[1][i]*dC[i];
        fC3 += ms[2][i]*dC[i];
    }

    AngleUpdate();
}

void Particle::AngleUpdate()
{
    fTheta = acos(fC3);
    fPhi = atan2(fC2, fC1);
}