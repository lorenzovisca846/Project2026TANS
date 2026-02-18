#include <iostream>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TEnv.h>

#include "classes/Cylinder.h"
#include "classes/SimRandom.h"
#include "classes/Particle.h"
#include "classes/MyPoint.h"

#define DISPLAY false

using namespace std;

/*
====================================================================================
POSSIBILI SCELTE PER LA GENERAZIONE CASUALE:

1) Vertice:
    - o         fissa (nell'origine)
    - u         gaussiana xy e uniforme z
    - g         gaussiana xyz

2) Molteplicità:
    - f         fissa (max)
    - h         da istogramma
    - u         uniforme

3) Rumore:
    - p         poissoniano
    - u         uniforme
    - f         fisso (rate)
====================================================================================
*/

typedef void (SimRandom::*vtxGen)(double&, double&, double&, double, double);
typedef int  (SimRandom::*mGen)(int, int);
typedef void (*nGen)(int, double, const Cylinder&, TClonesArray&, int&, SimRandom*);

void Transport(Particle* part, const Cylinder& layer, TClonesArray& hits, int& counter, bool detector, bool msEnabled);

void NoiseU(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand);
void NoiseP(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand);
void NoiseF(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand);

void FunctionAssignment(vtxGen& vptr, mGen& mptr, nGen& nptr, const string& gentypes);

int main(int argc, char** argv)
{
    //================================= Config parameters =================================
    string cFile = "fullConfig.txt";
    if (argc > 1) cFile = argv[1];

    #if DISPLAY
        cFile = "simDisplay.txt";
    #endif

    string configFile = "../config/" + cFile;

    TEnv *config = new TEnv(configFile.c_str());

    int Nevents         = config->GetValue("Events", 10000);
    string gentypes     = config->GetValue("Generation", "ghp");
    bool msEnabled      = config->GetValue("MultScattering", false);
    unsigned int seed   = config->GetValue("Seed", 0);

    int multMin         = config->GetValue("Minimum", 1);
    int multMax         = config->GetValue("Maximum", 69);

    bool noiseEnabled   = config->GetValue("NoiseEnabled", true);
    int noiseMax        = config->GetValue("MaxNoise", 20);
    double noiseRate    = config->GetValue("NoiseRate", 5.0);

    double vtxXYsigma   = config->GetValue("SigmaXY", 0.01);
    double vtxZsigma    = config->GetValue("SigmaZ", 5.3);
    double unifZedges   = config->GetValue("Zedges_uniform", 5.3);

    double bpR          = config->GetValue("BeamPipeRadius", 3.0);
    double bpL          = config->GetValue("Length", 27.0);
    double bpW          = config->GetValue("BeamPipeWidth", 0.08);

    double l1R          = config->GetValue("Layer1Radius", 4.0);
    double l1L          = config->GetValue("Length", 27.0);
    double l1W          = config->GetValue("LayerWidth", 0.02);

    double l2R          = config->GetValue("Layer2Radius", 7.0);
    double l2L          = config->GetValue("Length", 27.0);

    string bpMat        = config->GetValue("BeamPipeMaterial", "Be");
    string lMat         = config->GetValue("LayerMaterial", "Si");
    
    string inputN0      = config->GetValue("inputName", "kinem.root");
    string inputName    = "../config/" + inputN0; 
    string inputHM      = config->GetValue("inputHistoMult", "multHist");
    string inputHE      = config->GetValue("inputHistoEta", "etaHist");

    string outputN0     = config->GetValue("outputName", "simulation_output.root");
    string outputName   = "../outputs/" + outputN0;

    delete config;

    Cylinder beamPipe(bpR, bpL, bpW, bpMat, 0);
    Cylinder Layer1(l1R, l1L, l1W, lMat, 1);
    Cylinder Layer2(l2R, l2L, 0., lMat, 2);

    typedef struct{
        double X, Y, Z;
        int mult;} VTX;
    static VTX vertex;

    TFile *inputFile = new TFile(inputName.c_str(),"READ");
    TH1F *multHist= (TH1F*)inputFile->Get(inputHM.c_str()); 
    TH1F *etaHist= (TH1F*)inputFile->Get(inputHE.c_str()); 

    delete gRandom;
    SimRandom *simrand = new SimRandom(seed, multHist, etaHist, unifZedges);
    gRandom = simrand;

    inputFile->Close();
    delete inputFile;

    vtxGen VertGen;
    mGen MultGen;
    nGen NoiseGen;
    FunctionAssignment(VertGen, MultGen, NoiseGen, gentypes);


    TFile hfile(outputName.c_str(),"RECREATE");
    TTree *tree = new TTree("Tree_SimOut","Vertex-Hits TTree");

    int arrdim = multMax + noiseMax + 3; // safety margin

    TClonesArray *ptrhits1 = new TClonesArray("MyPoint",arrdim);
    TClonesArray &hits1 = *ptrhits1;

    TClonesArray *ptrhits2 = new TClonesArray("MyPoint",arrdim);
    TClonesArray &hits2 = *ptrhits2;

    tree->Branch("Vertex",&vertex.X,"X/D:Y:Z:mult/I");
    tree->Branch("Hits_L1",&ptrhits1);
    tree->Branch("Hits_L2",&ptrhits2);

    #if DISPLAY
        TClonesArray *ptrhitsBP = new TClonesArray("MyPoint",arrdim);
        TClonesArray &hitsBP = *ptrhitsBP;
        tree->Branch("Hits_BP",&ptrhitsBP);
    #endif

    //================================= Status print =================================
    string Vertexstr = "Gauss (default)";
    if(gentypes.length() > 0)
    {
        if(gentypes[0]=='o' || gentypes[0]=='O') Vertexstr = "Origin";
        if(gentypes[0]=='u' || gentypes[0]=='U') Vertexstr = "Uniform";
        if(gentypes[0]=='g' || gentypes[0]=='G') Vertexstr = "Gauss";
    }

    string Multstr = "From Histogram (default)";
    if(gentypes.length() > 1)
    {
        if(gentypes[1]=='f' || gentypes[1]=='F') Multstr = "Fixed";
        if(gentypes[1]=='u' || gentypes[1]=='U') Multstr = "Uniform";
        if(gentypes[1]=='h' || gentypes[1]=='H') Multstr = "From Histogram";
    }

    string Noisestr = "Poisson (default)";
    if(gentypes.length() > 2)
    {
        if(gentypes[2]=='u' || gentypes[2]=='U') Noisestr = "Uniform";
        if(gentypes[2]=='p' || gentypes[2]=='P') Noisestr = "Poisson";
        if(gentypes[2]=='f' || gentypes[2]=='F') Noisestr = "Fixed";
    }
    

    cout << "\n";
    cout << "=============== Simulation parameters ===============" << endl;
    cout << "  Number of events:         " << Nevents << endl;
    cout << "  Multiple scattering:      " << (msEnabled ? "YES" : "NO") << endl;
    cout << "  Vertex generation:        " << Vertexstr << endl;
    cout << "  Multiplicity generation:  " << Multstr << endl;
    cout << "  Noise enabled:            " << (noiseEnabled ? "YES" : "NO") << endl;
    if(noiseEnabled)
    cout << "  Noise generation:         " << Noisestr << endl;
    cout << "=====================================================" << endl;
    cout << "\n";

    #if DISPLAY
        cout << "================ Display mode enabled ================" << endl;
        cout << "\n";
    #endif

    //================================= Event loop =================================
    Particle *ptrPart = new Particle(simrand);
    int counter1, counter2, counterBP;

    TStopwatch timer;
    timer.Start();

    for(int i=0; i<Nevents; i++)
    {
        //================================= Vertex generation =================================
        if(i%10000==0) cout << "Simulating event " << i << "/" << Nevents << endl;
        
        vertex.mult = (simrand->*MultGen)(multMin, multMax);
        (simrand->*VertGen)(vertex.X, vertex.Y, vertex.Z, vtxXYsigma, vtxZsigma);

        counter1 = 0;
        counter2 = 0;
        counterBP = 0;

        for(int j=0; j<vertex.mult; j++)
        {
            //================================= Particle propagation =================================
            ptrPart->Init(vertex.X, vertex.Y, vertex.Z, 1.0, 0.7, j);

            #if DISPLAY
                Transport(ptrPart, beamPipe, hitsBP, counterBP, true, msEnabled);
            #else
                Transport(ptrPart, beamPipe, hits1, counter1, false, msEnabled);
            #endif
            Transport(ptrPart, Layer1, hits1, counter1, true, msEnabled);
            Transport(ptrPart, Layer2, hits2, counter2, true, false);
        }

        //================================= Noise generation =================================
        if(noiseEnabled)
        {
            NoiseGen(noiseMax, noiseRate, Layer1, hits1, counter1, simrand);
            NoiseGen(noiseMax, noiseRate, Layer2, hits2, counter2, simrand);
        }

        //================================= Tree fill =================================
        tree->Fill();
        ptrhits1->Clear();
        ptrhits2->Clear();
        #if DISPLAY
            ptrhitsBP->Clear();
        #endif
    }

    timer.Stop(); 
    timer.Print();

    hfile.Write();
    hfile.Close();

    delete ptrPart;
    delete simrand;
    delete ptrhits1;
    delete ptrhits2;
    #if DISPLAY
        delete ptrhitsBP;
    #endif

    return 0;
}

void Transport(Particle* part, const Cylinder& layer, TClonesArray& hits, int& counter, bool detector, bool msEnabled)
{
    part->Propagation(layer.GetR());

    if(fabs(part->GetZ()) < layer.GetL()*0.5)
    {
        if(detector)
        {
            new(hits[counter])MyPoint(layer.GetR(), part->GetPhi(), part->GetZ(), part->GetTrackID());
            counter++;
        }

        if(msEnabled)
            part->MultScatter(layer.GetX0(), layer.GetW(), layer.GetR());
    
    }
}

void NoiseU(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand)
{
    int nNoise = simrand->NoiseUnif(noiseRate, noiseMax);

    for(int i=0;i<nNoise;i++)
    {
        double phiNoise = simrand->PhiDist();
        double zNoise = simrand->ZDist(layer.GetL());

        new(hits[counter])MyPoint(layer.GetR(), phiNoise, zNoise);
        counter++;
    }
}

void NoiseP(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand)
{
    int nNoise = simrand->NoisePois(noiseRate, noiseMax);

    for(int i=0;i<nNoise;i++)
    {
        double phiNoise = simrand->PhiDist();
        double zNoise = simrand->ZDist(layer.GetL());

        new(hits[counter])MyPoint(layer.GetR(), phiNoise, zNoise);
        counter++;
    }
}

void NoiseF(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand)
{
    int nNoise = (int)noiseRate;

    for(int i=0;i<nNoise;i++)
    {
        double phiNoise = simrand->PhiDist();
        double zNoise = simrand->ZDist(layer.GetL());

        new(hits[counter])MyPoint(layer.GetR(), phiNoise, zNoise);
        counter++;
    }
}

void FunctionAssignment(vtxGen& vptr, mGen& mptr, nGen& nptr, const string& gentypes)
{
    char vt = (gentypes.length() > 0) ? gentypes[0] : ' ';
    char mt = (gentypes.length() > 1) ? gentypes[1] : ' ';
    char nt = (gentypes.length() > 2) ? gentypes[2] : ' ';

    switch(vt)
    {
        case 'o':
            vptr = &SimRandom::VertOrig;
            break;
        case 'u':
            vptr = &SimRandom::VertUnif;
            break;
        default:
            vptr = &SimRandom::VertGaus;
            break;
    }

    switch(mt)
    {
        case 'f':
            mptr = &SimRandom::MultFixed;
            break;
        case 'u':
            mptr = &SimRandom::MultUnif;
            break;
        default:
            mptr = &SimRandom::MultHisto;
            break;
    }

    switch(nt)
    {
        case 'u':
            nptr = &NoiseU;
            break;
        case 'f':
            nptr = &NoiseF;
            break;
        default:
            nptr = &NoiseP;
            break;
    }
}