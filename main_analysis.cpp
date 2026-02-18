#include <iostream>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
    string cFile = "inputConfig.txt";
    if (argc > 1) cFile = argv[1];
    string configFile = "../config/" + cFile;

    //================================= Input file reading =================================

    typedef struct {
        double Ztrue, Zrec;
        bool success;
        int mult;} REC; 
    static REC recVertex;

    TFile inputFile("../outputs/reconstruction_output.root","READ");
    TTree *inputTree = (TTree*)inputFile.Get("Tree_RecOut");
    TBranch *bVertex=inputTree->GetBranch("Vertex");
    bVertex->SetAddress(&recVertex.Ztrue);

    




    inputFile.Close();
    return 0;
}