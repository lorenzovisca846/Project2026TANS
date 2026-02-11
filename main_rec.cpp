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


    //================================= Reconstruction ==============================?







    inputFile.Close();
    outputFile.Close();

    return 0;
}