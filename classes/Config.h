#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <TObject.h>
#include <TEnv.h>
#include "SimRandom.h"
using namespace std;

class Config : public TObject
{

    public:

    Config(SimRandom* srnd, TEnv* configEnv);

    double beamPipeRadius;
    double beamPipeThickness;
    
    double layer1Radius;
    double layer1Thickness;
    
    double layer2Radius;
    double layer2Thickness;
    
    double detectorLength;

    double vertexZSigma;
    double vertexXYSigma;

    double smearZ;
    double smearRPhi;

    int nEvents;
    int multiplicityMin;
    int multiplicityMax;
    double noiseRateLayer;
    int noiseMaxLayer;

    string inputFileName;

    double deltaPhiCut;              // rad (massima differenza in phi per matching hits)
    double runningWindowSize;

    SimRandom* MyRandom() const {return fSimrand;}

    void Print();

    private:
        SimRandom* fSimrand;

    ClassDef(Config,1)
};



#endif // CONFIG_H