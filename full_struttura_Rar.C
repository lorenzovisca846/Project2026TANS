#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

// ============================================================
// configurazione
// ============================================================ 

struct config
{
    // geometria
    double beamPipeRadius = 3.0;           // cm // berillio    
    double beamPipeThickness = 0.08;       // cm
    double layer1Radius = 4.0;             // cm
    double layer1Thickness = 0.02;         // cm
    double layer2Radius = 7.0;             // cm
    double layer2Thickness = 0.02;         // cm
    double detectorLength = 27.0;          // cm

    // parametri vertice
    double vertexZSigma = 5.3;             // cm
    double vertexXYSigma = 0.01;           // cm

    // smearing
    double smearZ = 0.0120;                // cm
    double smearRPhi = 0.0030;             // cm
    
    // Materiali
    double X0_Be = 35.28;                  // cm
    double X0_Si = 9.37;                   // cm
    double mspConstant = 0.0136;           // MeV


    // parametri simulazione
    int nEvents = 1000;
    int multiplicityMin = 20;
    int multiplicityMax = 80;
    double noiseRateLayer = 5.0;           // %  

    // Metropolis
    double metropolisStepSize = 0.01;
    int metropolisNSteps = 10000;

    // Taglio in phi per formare tracklets
    double deltaPhiCut = 0.1;              // rad

    void Print()
    {
        cout << "=== CONFIGURAZIONE ===" << endl;
        cout << "Geometria: BeamPipe R=" << beamPipeRadius
             << "cm, Layer1 R=" << layer1Radius
             << "cm, Layer2 R=" << layer2Radius << "cm" << endl;
        cout << "Vertice: σ_z=" << vertexZSigma
             << "cm, σ_xy=" << vertexXYSigma << "cm" << endl;
        cout << "Delta Phi Cut: " << deltaPhiCut << " rad" << endl;
        cout << "Eventi: " << nEvents << endl;
    }
};

// ============================================================
// CLASSi dei dati
// ============================================================

class Particle
{
public: 
    int trackID;
    TVector3 position;
    TVector3 direction;
    double phi;
    double theta;
    double p;
    double pt;
    int charge;
    double beta;    

    // Costruttore di Particle
    Particle() : trackID(-1), theta(0), phi(0), p(0), pt(0), charge(0), beta(1.0) {}
};

class Hit
{
public:
    int layer;
    double x;
    double y;
    double z;
    double r;
    double phi;
    double z_smeared;
    double phi_smeared;
    bool isSignal;
    int trackID;   

    // Costruttore di Hit
    Hit() : layer(0), x(0), y(0), z(0), r(0), phi(0), z_smeared(0), phi_smeared(0), isSignal(true), trackID(-1) {}
};

class Tracklet
{
public:
    Hit hit1;
    Hit hit2;
    double z_intersection;
    double slope;
    double chi2;
    bool valid;

    // Costruttore di Tracklet
    Tracklet() : z_intersection(0), slope(0), chi2(0), valid(false) {}
};

class Event
{
public:
    int eventID;
    int multiplicity;
    double trueVertexZ;
    double recoVertexZ;
    double recoVertexZMetropolis;
    bool vertexFound;
    int nTracklets;

    // Costruttore di Event
    Event() : eventID(0), multiplicity(0), trueVertexZ(0), recoVertexZ(0), recoVertexZMetropolis(0), vertexFound(false), nTracklets(0) {}

};

// ============================================================
// CLASSI DETECTOR
// ============================================================

class Detector
{
private:
    config fConfig;

public:
    // Costruttore di Detector, usa una lista di inizializzazione
    Detector(const config& config) : fConfig(config) {} 

    bool acceptance(const TVector3& punto , int layer) const
    {
      sesso caldo bollente omosessuale
      poi ci penso   
    }
};

// ============================================================
// CLASSE SCATTERING
// ============================================================

class MultipleScattering {
private:
    Config fConfig;
    TRandom3 fRandom;  

public:
    // Costruttore di MultipleScattering, usa una lista di inizializzazione
    MultipleScattering(const Config& config) : fConfig(config), fRandom(0)

    double ScatterAngle(double p, double beta, double thickness, double charge , const string& material) const
    {
        double X0;
        if (material == "Be") {
            X0 = fConfig.X0_Be;
        } else if (material == "Si") {
            X0 = fConfig.X0_Si;
        } else {
            cerr << "Material not recognized!" << endl;
            return 0.0;
        }
        //  codice per double theta0
    }

    void DajeCoStoScatter(Particle& particle, double thickness, const string& \material)
    {
        //qualcoda che prende theta0 e applica lo scattering alla particella in maniera gaussiana
        // mentre phi di scattering e' uniforme perchè simmetria cilindrica
    }
    
};

// ============================================================
// E CLASSI PAA TRACCIA
// ============================================================

class TrackIntersector {
private:
    Config fConfig; 

public:
    // Costruttore di TrackIntersector, usa una lista di inizializzazione   
    TrackIntersector(const Config& config) : fConfig(config) {}

    //CODICE PE CARCOLA LE NTERSEZIONI DAE SLIDE DE MASSIMO7
};

// ============================================================
// CLASSE METROPOLIS
// ============================================================

class MetropolisVertexReconstructor {
private:
    Config fConfig;
    TRandom3 fRandom;
    double fVertexZ;
    double fCurrentLLH; // log-likelihood corrente questo mea suggerito copailot
    double fStepSize;
    int fIterations;
    double fBestZ; 
    double fBestLLH;

public:
    MetropolisVertexReconstructor(double stepSize = 0.01, int iterations = 5000) : fRandom(0), fStepSize(stepSize), fIterations(iterations),  fBestZ(0), fBestLLH(-1e9) {}

    double FitVertex(stocazzo)
    {
        // codice per fit del vertice con metropolis
    }

    double CalculateLLH(stocazzo)
    {
        // codice per calcolo log-likelihood
    }
};

// ============================================================
// CLASSE RICOSTRUTTORE
// ============================================================

class VertexReconstructor {
private:
    Config fConfig;
    TRandom3 fRandom;

public:
    VertexReconstructor(const Config& config) : fConfig(config), fRandom(0)

    //qualcosa 
    double SimpleFit(const vector<Tracklet>& tracklets)
    {
        // codice per fit semplice del vertice
    }

    void TrackletIntersection(Tracklet& tracklet) const
    {
        // codice per calcolo intersezione del tracklet
    }

};