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

    SimRandom* MyRandom() const {return fSimrand;}

    double GetBPRadius() const {return beamPipeRadius;}
    double GetBPThickness() const {return beamPipeThickness;}
    string GetBPMaterial() const {return beamPipeMaterial;}

    double GetL1Radius() const {return layer1Radius;}
    double GetL1Thickness() const {return layer1Thickness;}

    double GetL2Radius() const {return layer2Radius;}
    double GetL2Thickness() const {return layer2Thickness;}

    double GetDetectorLength() const {return detectorLength;}
    string GetLayerMaterial() const {return layerMaterial;}

    double GetVTXZSigma() const {return vertexZSigma;}
    double GetVTXXYSigma() const {return vertexXYSigma;}
    double GetVTXZEdges() const {return vertexZedges;}

    double GetSmearZ() const {return smearZ;}
    double GetSmearRPhi() const {return smearRPhi;}

    int GetNEvents() const {return nEvents;}
    int GetMultMin() const {return multiplicityMin;}
    int GetMultMax() const {return multiplicityMax;}
    string GetGenTypes() const {return gentypes;}

    bool IsMSEnabled() const {return msEnabled;}

    bool IsNoiseEnabled() const {return noiseEnabled;}
    double GetNoiseRateLayer() const {return noiseRateLayer;}
    int GetNoiseMaxLayer() const {return noiseMaxLayer;}

    string GetInputFileName() const {return inputFileName;}

    double GetDeltaPhiCut() const {return deltaPhiCut;}
    double GetRunningWindowSize() const {return runningWindowSize;}

    int GetMultMinZoom() const {return multminZoom;}
    int GetMultMaxZoom() const {return multmaxZoom;}

    bool DisplayErrFull() const {return displayerrfull;}
    bool DisplayErrSelect() const {return displayerrselect;}
    double GetErrZLimit() const {return errZlimit;}

    bool DisplayEffFull() const {return displayefffull;}
    bool DisplayEff1Sigma() const {return displayeff1sigma;}
    bool DisplayEff3Sigma() const {return displayeff3sigma;}
    
    bool DisplayResFull() const {return displayresfull;}
    bool DisplayRes1Sigma() const {return displayres1sigma;}
    bool DisplayRes3Sigma() const {return displayres3sigma;}

    bool DisplayEffZvrt() const {return displayeffZvrt;}
    bool DisplayResZvrt() const {return displayresZvrt;}
    
    void Print();

    private:
        SimRandom* fSimrand;
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

    ClassDef(Config,1)
};



#endif // CONFIG_H