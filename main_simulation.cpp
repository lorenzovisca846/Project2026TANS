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
#include "classes/Config.h"

#define DISPLAY false

using namespace std;

typedef void (SimRandom::*vtxGen)(double&, double&, double&, double, double);
typedef int  (SimRandom::*mGen)(int, int);
typedef int (SimRandom::*nGen)(double, int);

void Transport(Particle*, const Cylinder&, TClonesArray&, int&, bool, bool);
void Noise(const Cylinder&, TClonesArray&, int&, nGen&, Config&);

void FunctionAssignment(vtxGen&, mGen&, nGen&, const string&);

int main(int argc, char** argv)
{
    //================================= Config parameters =================================
    string cFile = "fullConfig.txt";
    if (argc > 1) cFile = argv[1];

    #if DISPLAY
        cFile = "simDisplay.txt";
    #endif

    string configFile = "../config/" + cFile;

    TEnv *configEnv = new TEnv(configFile.c_str());
    unsigned int seed   = configEnv->GetValue("Seed", 0);
    double unifZedges   = configEnv->GetValue("Zedges_uniform", 5.3);
    
    string inputN0      = configEnv->GetValue("inputName", "kinem.root");
    string inputName    = "../config/" + inputN0; 
    string inputHM      = configEnv->GetValue("inputHistoMult", "multHist");
    string inputHE      = configEnv->GetValue("inputHistoEta", "etaHist");
    string outputN0     = configEnv->GetValue("outputName", "simulation_output.root");
    string outputName   = "../outputs/" + outputN0;

    TFile *inputFile = new TFile(inputName.c_str(),"READ");
    TH1F *multHist= (TH1F*)inputFile->Get(inputHM.c_str()); 
    TH1F *etaHist= (TH1F*)inputFile->Get(inputHE.c_str()); 

    delete gRandom;
    SimRandom *simrand = new SimRandom(seed, multHist, etaHist, unifZedges);
    gRandom = simrand;

    inputFile->Close();
    delete inputFile;

    Config config(simrand, configEnv);

    delete configEnv;

    Cylinder beamPipe(config.beamPipeRadius, config.detectorLength, config.beamPipeThickness, config.beamPipeMaterial, 0);
    Cylinder Layer1(config.layer1Radius, config.detectorLength, config.layer1Thickness, config.layerMaterial, 1);
    Cylinder Layer2(config.layer2Radius, config.detectorLength, 0., config.layerMaterial, 2);

    typedef struct{
        double X, Y, Z;
        int mult;} VTX;
    static VTX vertex;

    vtxGen VertGen;
    mGen MultGen;
    nGen NoiseGen;
    FunctionAssignment(VertGen, MultGen, NoiseGen, config.gentypes);


    TFile hfile(outputName.c_str(),"RECREATE");
    TTree *tree = new TTree("Tree_SimOut","Vertex-Hits TTree");

    int arrdim = config.multiplicityMax + config.noiseMaxLayer + 3;

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
    string gt = config.gentypes;
    if(gt.length() > 0)
    {
        if(gt[0]=='o' || gt[0]=='O') Vertexstr = "Origin";
        if(gt[0]=='u' || gt[0]=='U') Vertexstr = "Uniform";
        if(gt[0]=='g' || gt[0]=='G') Vertexstr = "Gauss";
    }

    string Multstr = "From Histogram (default)";
    if(gt.length() > 1)
    {
        if(gt[1]=='f' || gt[1]=='F') Multstr = "Fixed";
        if(gt[1]=='u' || gt[1]=='U') Multstr = "Uniform";
        if(gt[1]=='h' || gt[1]=='H') Multstr = "From Histogram";
    }

    string Noisestr = "Poisson (default)";
    if(gt.length() > 2)
    {
        if(gt[2]=='u' || gt[2]=='U') Noisestr = "Uniform";
        if(gt[2]=='p' || gt[2]=='P') Noisestr = "Poisson";
        if(gt[2]=='f' || gt[2]=='F') Noisestr = "Fixed";
    }
    

    cout << "\n";
    cout << "=============== Simulation parameters ===============" << endl;
    cout << "  Number of events:         " << config.nEvents << endl;
    cout << "  Multiple scattering:      " << (config.msEnabled ? "YES" : "NO") << endl;
    cout << "  Vertex generation:        " << Vertexstr << endl;
    cout << "  Multiplicity generation:  " << Multstr << endl;
    cout << "  Noise enabled:            " << (config.noiseEnabled ? "YES" : "NO") << endl;
    if(config.noiseEnabled)
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

    for(int i=0; i<config.nEvents; i++)
    {
        //================================= Vertex generation =================================
        if(i%10000==0) cout << "Simulating event " << i << "/" << config.nEvents << endl;
        
        vertex.mult = (simrand->*MultGen)(config.multiplicityMin, config.multiplicityMax);
        (simrand->*VertGen)(vertex.X, vertex.Y, vertex.Z, config.vertexXYSigma, config.vertexZSigma);

        counter1 = 0;
        counter2 = 0;
        counterBP = 0;

        for(int j=0; j<vertex.mult; j++)
        {
            //================================= Particle propagation =================================
            ptrPart->Init(vertex.X, vertex.Y, vertex.Z, 1.0, 0.7, j);

            #if DISPLAY
                Transport(ptrPart, beamPipe, hitsBP, counterBP, true, config.msEnabled);
            #else
                Transport(ptrPart, beamPipe, hits1, counter1, false, config.msEnabled);
            #endif
            Transport(ptrPart, Layer1, hits1, counter1, true, config.msEnabled);
            Transport(ptrPart, Layer2, hits2, counter2, true, false);
        }

        //================================= Noise generation =================================
        if(config.noiseEnabled)
        {
            Noise(Layer1, hits1, counter1, NoiseGen, config);
            Noise(Layer2, hits2, counter2, NoiseGen, config);
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

    delete simrand;

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

void Noise(const Cylinder& layer, TClonesArray& hits, int& counter, nGen& noiseFunc, Config& config)
{
    int nNoise = (config.MyRandom()->*noiseFunc)(config.noiseRateLayer, config.noiseMaxLayer);

    for(int i=0;i<nNoise;i++)
    {
        double phiNoise = config.MyRandom()->PhiDist();
        double zNoise = config.MyRandom()->ZDist(layer.GetL());

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
            nptr = &SimRandom::NoiseUnif;
            break;
        case 'f':
            nptr = &SimRandom::NoiseFixed;
            break;
        default:
            nptr = &SimRandom::NoisePois;
            break;
    }
}