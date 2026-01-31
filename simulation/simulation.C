#include <Riostream.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TStopwatch.h>

#include "simclasses/Cylinder.h"
#include "simclasses/SimRandom.h"
#include "simclasses/Particle.h"
#include "simclasses/MyPoint.h"

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
====================================================================================
*/

typedef void (SimRandom::*vtxGen)(double&, double&, double&, double, double);
typedef int  (SimRandom::*mGen)(int, int);
typedef void (*nGen)(int, double, const Cylinder&, TClonesArray&, int&, SimRandom*);

void Transport(Particle* part, const Cylinder& layer, TClonesArray& hits, int& counter, bool detector, bool msEnabled);

void NoiseU(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand);
void NoiseP(int noiseMax, double noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand);

void FunctionAssignment(vtxGen& vptr, mGen& mptr, nGen& nptr, const string& gentypes);

void simulation(double Nevents = 10000, bool msEnabled = false, string gentypes = "ghp", unsigned int seed = 0)
{
    //================================= Config parameters =================================
    int multMin = 0;
    int multMax = 60;
    const int noiseMax = 20;
    const double noiseRate = 5.;

    const double vtxXYsigma = 0.01;
    const double vtxZsigma = 5.3;

    const double bpR = 3.0, bpL = 27.0, bpW = 0.08;
    const double l1R = 4.0, l1L = 27.0, l1W = 0.02;
    const double l2R = 7.0, l2L = 27.0;

    Cylinder beamPipe(bpR, bpL, bpW, "Be");
    Cylinder Layer1(l1R, l1L, l1W, "Si");
    Cylinder Layer2(l2R, l2L, 0., "Si");

    typedef struct{
        double X, Y, Z;
        int mult;} VTX;
    static VTX vertex;

    //================================= Input file for distributions =================================
    TFile *inputFile = new TFile("inputDistributions.root","READ");
    TH1F *multHist= (TH1F*)inputFile->Get("multHist"); 
    TH1F *etaHist= (TH1F*)inputFile->Get("etaHist"); 

    delete gRandom;
    SimRandom *simrand = new SimRandom(seed, multHist, etaHist);
    gRandom = simrand;

    inputFile->Close();
    delete inputFile;

    vtxGen VertGen;
    mGen MultGen;
    nGen NoiseGen;
    FunctionAssignment(VertGen, MultGen, NoiseGen, gentypes);

    //================================= Output file and tree =================================
    TFile hfile("htree.root","RECREATE");
    TTree *tree = new TTree("Tree","Vertex-Hits TTree");

    int arrdim = multMax + noiseMax + 5;

    TClonesArray *ptrhits1 = new TClonesArray("MyPoint",arrdim);
    TClonesArray &hits1 = *ptrhits1;

    TClonesArray *ptrhits2 = new TClonesArray("MyPoint",arrdim);
    TClonesArray &hits2 = *ptrhits2;

    tree->Branch("Vertex",&vertex.X,"X/D:Y:Z:mult/I");
    tree->Branch("Hits_L1",&ptrhits1);
    tree->Branch("Hits_L2",&ptrhits2);

    //================================= Status print =================================
    string Vertexstr = "Gaus (default)";
    if(gentypes.length() > 0)
    {
        if(gentypes[0]=='o' || gentypes[0]=='O') Vertexstr = "Origin";
        if(gentypes[0]=='u' || gentypes[0]=='U') Vertexstr = "Uniform";
        if(gentypes[0]=='g' || gentypes[0]=='G') Vertexstr = "Gaus";
    }

    string Multstr = "Histo (default)";
    if(gentypes.length() > 1)
    {
        if(gentypes[1]=='f' || gentypes[1]=='F') Multstr = "Fixed";
        if(gentypes[1]=='u' || gentypes[1]=='U') Multstr = "Uniform";
        if(gentypes[1]=='h' || gentypes[1]=='H') Multstr = "Histo";
    }

    string Noisestr = "Poisson (default)";
    if(gentypes.length() > 2)
    {
        if(gentypes[2]=='u' || gentypes[2]=='U') Noisestr = "Uniform";
        if(gentypes[2]=='p' || gentypes[2]=='P') Noisestr = "Poisson";
    }
    

    cout << "\n";
    cout << "=============== Simulation parameters ===============" << endl;
    cout << "  Number of events:         " << Nevents << endl;
    cout << "  Multiple scattering:      " << (msEnabled ? "YES" : "NO") << endl;
    cout << "  Vertex generation:        " << Vertexstr << endl;
    cout << "  Multiplicity generation:  " << Multstr << endl;
    cout << "  Noise generation:         " << Noisestr << endl;
    cout << "=====================================================" << endl;
    cout << "\n";


    //================================= Event loop =================================
    Particle *ptrPart = new Particle(simrand);
    int counter1, counter2;

    TStopwatch timer;
    timer.Start();

    for(int i=0; i<Nevents; i++)
    {
        //================================= Vertex generation =================================
        if(i%10000==0) cout << "Processing event " << i << "/" << Nevents << endl;
        
        vertex.mult = (simrand->*MultGen)(multMin, multMax);
        (simrand->*VertGen)(vertex.X, vertex.Y, vertex.Z, vtxXYsigma, vtxZsigma);

        counter1 = 0;
        counter2 = 0;

        for(int j=0; j<vertex.mult; j++)
        {
            //================================= Particle propagation =================================
            ptrPart->Init(vertex.X, vertex.Y, vertex.Z, 1.0, 0.7, j);

            Transport(ptrPart, beamPipe, hits1, counter1, false, msEnabled);
            Transport(ptrPart, Layer1, hits1, counter1, true, msEnabled);
            Transport(ptrPart, Layer2, hits2, counter2, true, false);
        }

        //================================= Noise generation =================================
        if(noiseMax>0)
        {
            NoiseGen(noiseMax, noiseRate, Layer1, hits1, counter1, simrand);
            NoiseGen(noiseMax, noiseRate, Layer2, hits2, counter2, simrand);
        }

        //================================= Tree fill =================================

        tree->Fill();
        ptrhits1->Clear();
        ptrhits2->Clear();
    }

    timer.Stop(); 
    timer.Print();

    hfile.Write();
    hfile.Close();

    //Check if these are needed
    delete ptrPart;
    delete simrand;
    delete ptrhits1;
    delete ptrhits2;
}

void Transport(Particle* part, const Cylinder& layer, TClonesArray& hits, int& counter, bool detector, bool msEnabled)
{
    part->Propagation(layer.GetR());

    if(fabs(part->GetZ()) < layer.GetL()*0.5)
    {
        if(detector)
        {
            new(hits[counter])MyPoint(part->GetPoint());
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

        double xNoise = layer.GetR() * cos(phiNoise);
        double yNoise = layer.GetR() * sin(phiNoise);

        new(hits[counter])MyPoint(xNoise, yNoise, zNoise);
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

        double xNoise = layer.GetR() * cos(phiNoise);
        double yNoise = layer.GetR() * sin(phiNoise);

        new(hits[counter])MyPoint(xNoise, yNoise, zNoise);
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
        default:
            nptr = &NoiseP;
            break;
    }
}
