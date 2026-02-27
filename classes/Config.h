#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <TObject.h>
#include <TEnv.h>
#include <string>
#include "SimRandom.h"
using namespace std;

class Config : public TObject
{

    public:

    Config(SimRandom* srnd, TEnv* configEnv);

    double beamPipeRadius;
    double beamPipeThickness;
    string beamPipeMaterial;
    
    double layer1Radius;
    double layer1Thickness;
    
    double layer2Radius;
    double layer2Thickness;
    
    double detectorLength;
    string layerMaterial;

    double vertexZSigma;
    double vertexXYSigma;
    double vertexZedges;

    double smearZ;
    double smearRPhi;

    int nEvents;
    int multiplicityMin;
    int multiplicityMax;
    string gentypes;

    bool msEnabled;

    bool noiseEnabled;
    double noiseRateLayer;
    int noiseMaxLayer;

    string inputFileName;

    double deltaPhiCut;
    double runningWindowSize;

    int multminZoom;
    int multmaxZoom;

    bool displayerrfull;
    bool displayerrselect;
    double errZlimit;

    bool displayefffull;
    bool displayeff1sigma;
    bool displayeff3sigma;

    bool displayresfull;
    bool displayres1sigma;
    bool displayres3sigma;

    bool displayeffZvrt;
    bool displayresZvrt;

    SimRandom* MyRandom() const {return fSimrand;}

    void Print();

    private:
        SimRandom* fSimrand;

    ClassDef(Config,1)
};



#endif // CONFIG_H