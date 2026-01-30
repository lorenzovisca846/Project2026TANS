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
    - VertOrig          fissa (nell'origine)
    - VertUnif          gaussiana xy e uniforme z
    - VertGaus          gaussiana xyz

2) Molteplicità:
    - MultFixed         fissa (max)
    - MultHisto         da istogramma
    - MultUnif          uniforme

3) Rumore:
    - NoisePois         poissoniano
    - NoiseUnif         uniforme

====================================================================================
*/

void Transport(Particle* part, const Cylinder& layer, TClonesArray& hits, int& counter, bool detector, bool msEnabled);
void Noise(const int& noiseMax, const double& noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand);

void simulation(double Nevents = 100, bool msEnabled = false, unsigned int seed = 0)
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

    //================================= Output file and tree =================================
    TFile hfile("htree.root","RECREATE");
    TTree *tree = new TTree("Tree","Vertex-Hits TTree");

    int arrdim = Nevents * (multMax + noiseMax + 5);

    TClonesArray *ptrhits1 = new TClonesArray("MyPoint",arrdim);
    TClonesArray &hits1 = *ptrhits1;

    TClonesArray *ptrhits2 = new TClonesArray("MyPoint",arrdim);
    TClonesArray &hits2 = *ptrhits2;

    tree->Branch("Vertex",&vertex.X,"X/D:Y:Z:mult/I");
    tree->Branch("Hits_L1",&ptrhits1);
    tree->Branch("Hits_L2",&ptrhits2);


    cout << "================================ Simulation parameters =================================" << endl;
    cout << "Number of events:      " << Nevents << endl;
    cout << "Multiple scattering:   " << (msEnabled ? "YES" : "NO") << endl;

    //================================= Event loop =================================

    Particle *ptrPart = new Particle(simrand);
    int counter1, counter2;

    TStopwatch timer;
    timer.Start();

    for(unsigned int i=0; i<Nevents; i++)
    {
        //================================= Vertex generation =================================

        if(i%1000==0) cout << "Processing event " << i << "/" << Nevents << endl;
        
        vertex.mult = simrand->MultHisto(multMin, multMax);
        simrand->VertUnif(vertex.X, vertex.Y, vertex.Z, vtxXYsigma, vtxZsigma);

        counter1 = 0;
        counter2 = 0;

        for(unsigned int j=0; j<vertex.mult; j++)
        {
            //================================= Particle propagation =================================

            ptrPart->Init(vertex.X, vertex.Y, vertex.Z, 1.0, 0.7);

            Transport(ptrPart, beamPipe, hits1, counter1, false, msEnabled);
            Transport(ptrPart, Layer1, hits1, counter1, true, msEnabled);
            Transport(ptrPart, Layer2, hits2, counter2, true, false);
        }

        //================================= Noise generation =================================

        if(noiseMax>0)
        {
            Noise(noiseMax, noiseRate, Layer1, hits1, counter1, simrand);
            Noise(noiseMax, noiseRate, Layer2, hits2, counter2, simrand);
        }

        //================================= Tree fill =================================

        tree->Fill();
        ptrhits1->Clear();
        ptrhits2->Clear();
    }

    timer.Stop(); 
    timer.Print();

    //================================= Cleanup =================================

    hfile.Write();
    hfile.Close();
    delete ptrhits1;
    delete ptrhits2;

}

void Transport(Particle* part, const Cylinder& layer, TClonesArray& hits, int& counter, bool detector, bool msEnabled)
{
    part->Propagation(layer.GetR());

    if(fabs(part->GetZ()) < layer.GetL()/2.)
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

void Noise(const int& noiseMax, const double& noiseRate, const Cylinder& layer, TClonesArray& hits, int& counter, SimRandom* simrand)
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
