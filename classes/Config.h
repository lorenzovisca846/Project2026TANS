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
    int multiplicityMax ;
    double noiseRateLayer;

    // Parametri per l'algoritmo Metropolis
    double metropolisStepSize;      // cm (dimensione del passo per Metropolis)
    int metropolisNSteps;          // Numero di passi per Metropolis

    // Taglio per formare i tracklets
    double deltaPhiCut;              // rad (massima differenza in phi per matching hits)

    void Print();

    private:
        SimRandom* fSimrand;

    ClassDef(Config,1)
};



#endif // CONFIG_H