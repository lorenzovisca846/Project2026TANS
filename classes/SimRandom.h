#ifndef SIMRANDOM_H
#define SIMRANDOM_H

#include <TRandom3.h>
#include <TH1F.h>
#include <TMath.h>

class SimRandom : public TRandom3 {
    public:
        SimRandom();
        SimRandom(unsigned int seed, TH1F* multHist, TH1F* etaHist);

        ~SimRandom();

        int VMult1(int min, int max) {return min + Integer(max - min + 1);}
        int VMult2(int min, int max);

        void GausPoint(double&, double&, double&, double, double);
        void OriginPoint(double&, double&, double&);
        void UnifPoint(double&, double&, double&, double, double);

        double PhiDist() {return Rndm()*2.*TMath::Pi();}
        double ThetaDist();

    private:
    unsigned int fSeed;
    TH1F* fMultHist;
    TH1F* fEtaHist;

    ClassDef(SimRandom,1);
};

#endif // SIMRANDOM_H