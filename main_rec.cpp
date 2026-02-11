#include <TClonesArray.h>
#include <TEnv.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <string>

#include "classes/Config.h"
#include "classes/MyPoint.h"

int main(int argc, char** argv)
{

    //================================= Config parameters =================================
    string cFile = "inputConfig.txt";
    if (argc > 1) cFile = argv[1];
    string configFile = "../config/" + cFile;

    TEnv *configEnv = new TEnv(configFile.c_str());

    SimRandom* simrand = new SimRandom(configEnv->GetValue("Seed",0));
    delete gRandom;
    gRandom = simrand;

    Config config(simrand, configEnv);
    config.Print();

    delete configEnv;



    //================================= Input file reading =================================

    typedef struct{
        double X, Y, Z;
        int mult;} VTX;
    static VTX trueVertex;

    int arrdim = config.multiplicityMax + config.noiseMaxLayer + 3;

    TClonesArray *ptrhits1 = new TClonesArray("MyPoint",arrdim);
    TClonesArray *ptrhits2 = new TClonesArray("MyPoint",arrdim);


    TFile inputFile("../outputs/simulation_output.root","READ");

    TTree *inputTree = (TTree*)inputFile.Get("Tree_SimOut");
    TBranch *bVertex=inputTree->GetBranch("Vertex");
    TBranch *bLayer1=inputTree->GetBranch("Hits_L1");
    TBranch *bLayer2=inputTree->GetBranch("Hits_L2");

    bVertex->SetAddress(&trueVertex.X);
    bLayer1->SetAddress(&ptrhits1);
    bLayer2->SetAddress(&ptrhits2);

    //================================= Output file =================================

    typedef struct {
        double Ztrue, Zrec;
        bool succes;
        int mult;} REC; 
    static REC recVertex;

    string outName = "../outputs/reconstruction_output";
    string metroName = (config.MetropolisUsed) ? "_metropolis" : "";
    string ext = ".root";


    TFile outputFile((outName+metroName+ext).c_str(),"RECREATE");

    TTree *outputTree = new TTree("Tree_RecOut","Reconstruction Tree");
    outputTree->Branch("Vertex",&recVertex.Ztrue,"Ztrue/D:Zrec:succes/O:mult/I");


    //================================= Reconstruction ==============================

    vector<MyPoint> hitsL1;
    vector<MyPoint> hitsL2;

    int Nevents = inputTree->GetEntries();

    //================================= Normal reconstruction ==============================
    for(int i=0; i<Nevents;i++)
    {
        inputTree->GetEvent(i);

        int nHitsL1 = ptrhits1->GetEntries();
        hitsL1.clear();
        hitsL1.reserve(nHitsL1);

        for(int j=0; j<nHitsL1; j++)
        {
            hitsL1.push_back(*(MyPoint*)ptrhits1->At(j));
            config.MyRandom()->Smear(hitsL1.back(), config.smearZ, config.smearRPhi);
        }

        int nHitsL2 = ptrhits2->GetEntries();
        hitsL2.clear();
        hitsL2.reserve(nHitsL2);

        for(int j=0; j<nHitsL2; j++)
        {
            hitsL2.push_back(*(MyPoint*)ptrhits2->At(j));
            config.MyRandom()->Smear(hitsL2.back(), config.smearZ, config.smearRPhi);
        }






        
    }

    inputFile.Close();
    outputFile.Close();

    return 0;
}