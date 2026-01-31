#ifndef SIMRANDOM_H
#define SIMRANDOM_H
#include <TRandom3.h>
#include <TH1F.h>
#include <TMath.h>

class SimRandom : public TRandom3 {

    public:
        SimRandom():TRandom3(0),fSeed(0),fEtaHist(nullptr),fMultHist(nullptr){}
        SimRandom(unsigned int seed, TH1F* multHist, TH1F* etaHist);

        void VertGaus(double&, double&, double&, double, double);
        void VertOrig(double&, double&, double&, double, double);
        void VertUnif(double&, double&, double&, double, double);

        int MultUnif(int min, int max) {return min + Integer(max - min + 1);}
        int MultHisto(int min, int max) {return fMultHist->GetRandom();}
        int MultFixed(int min, int max) {return max;}

        double PhiDist() {return Rndm()*2.*TMath::Pi();}
        double ThetaDist();
        double ScatterDist(double p, double beta, double len);

        int NoisePois(double rate, int max);
        int NoiseUnif(double rate, int max) {return Integer(max + 1);}

        double ZDist(const double& zLen) {return (Rndm()-0.5)*zLen;}

    private:
    unsigned int fSeed;
    TH1F* fMultHist;
    TH1F* fEtaHist;

    ClassDef(SimRandom,1);
};

#endif // SIMRANDOM_H