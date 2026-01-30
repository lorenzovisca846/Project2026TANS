#include "SimRandom.h"

ClassImp(SimRandom)

SimRandom::SimRandom() : TRandom3(0),
    fSeed(0),
    fEtaHist(nullptr),
    fMultHist(nullptr)
{

}

SimRandom::SimRandom(unsigned int seed, TH1F* multHist, TH1F* etaHist) : TRandom3(seed),
    fSeed(seed),
    fMultHist(nullptr),
    fEtaHist(nullptr)
{
    fMultHist = (TH1F*)multHist->Clone("fMultHist_Sim");
    fMultHist->SetDirectory(nullptr);
    
    fEtaHist = (TH1F*)etaHist->Clone("fEtaHist_Sim");
    fEtaHist->SetDirectory(nullptr);
}

SimRandom::~SimRandom()
{
    delete fMultHist;
    delete fEtaHist;
}

int SimRandom::VMult2(int min, int max)
{
    return fMultHist->GetRandom();
}

void SimRandom::GausPoint(double& x, double& y, double& z, double xyS, double zS)
{
    x = Gaus(0., xyS);
    y = Gaus(0., xyS);
    z = Gaus(0., zS);
}

void SimRandom::OriginPoint(double& x, double& y, double& z)
{
    x = 0.;
    y = 0.;
    z = 0.;
}

void SimRandom::UnifPoint(double& x, double& y, double& z, double xyS, double zS)
{
    x = Gaus(0., xyS);
    y = Gaus(0., xyS);
    z = -zS + Rndm()*2.*zS;
}

double SimRandom::ThetaDist()
{
    double eta = fEtaHist->GetRandom();
    return 2.*atan(exp(-eta));
}

double SimRandom::ScatterDist(double p, double beta, double len)
{
    double theta0 = (0.0136/(beta*p)) * sqrt(len) * (1. + 0.038*log(len));
    return Gaus(0., theta0);
}

int SimRandom::RateDist1(double rate, int max)
{
    int result;
    do
    {
        result = Poisson(rate);
    } while (result>max);
    return result;
}

int SimRandom::RateDist2(int max)
{
    return Integer(max + 1);
}