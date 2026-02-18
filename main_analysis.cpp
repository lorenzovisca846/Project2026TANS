#include <iostream>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TEnv.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <string>

using namespace std;

void DisplayHisto(TH1F* histo, string title, string filename)
{
    TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 800, 600);
    histo->Draw("PLE");
    histo->Fit("gaus", "Q");
    gStyle->SetOptStat(0001);
    c->SaveAs(("../outputs/plots" + filename).c_str());
}

int main(int argc, char** argv)
{
    string cFile = "analysis_params.txt";
    if (argc > 1) cFile = argv[1];
    string configFile = "../config/" + cFile;

    TEnv *config = new TEnv(configFile.c_str());

    int multMin             = config->GetValue("MultiplicityMin",0);
    int multMax             = config->GetValue("MultiplicityMax",70);

    bool displayerrfull     = config->GetValue("Residuals_all_mult", false);
    bool displayerrselect   = config->GetValue("Residuals_select_mult", false);
    double errZlimit        = config->GetValue("Residuals_zlim", 0.1);   

    bool displayefffull     = config->GetValue("Efficiency_mult_allZ", false);
    bool displayeff1sigma   = config->GetValue("Efficiency_mult_1sigma", false);
    bool displayeff3sigma   = config->GetValue("Efficiency_mult_3sigma", false);

    bool displayresfull     = config->GetValue("Resolution_mult_allZ", false);
    bool displayres1sigma   = config->GetValue("Resolution_mult_1sigma", false);
    bool displayres3sigma   = config->GetValue("Resolution_mult_3sigma", false);

    bool displayeffZvrt     = config->GetValue("Efficiency_Zvert", false);
    bool displayresZvrt     = config->GetValue("Resolution_Zvert", false);

    char Zdistribution      = config->GetValue("ZDist", 'g');
    char multdistribution   = config->GetValue("MultDist", 'h');

    delete config;

    //================================= Input file =================================

    typedef struct {
        double Ztrue, Zrec;
        bool success;
        int mult;} REC; 
    static REC recVertex;

    TFile inputFile("../outputs/reconstruction_output.root","READ");
    TTree *inputTree = (TTree*)inputFile.Get("Tree_RecOut");
    TBranch *bVertex=inputTree->GetBranch("Vertex");
    bVertex->SetAddress(&recVertex.Ztrue);

    string selectedmult = "(" + to_string(multMin) + " < mult < " + to_string(multMax) + ")";
    int nbinM = multMax - multMin + 1;


    TH1F* FullErrHisto = new TH1F("FullErrHisto", "Residuals;Z_{rec}-Z_{true} (cm);Entries", 200, -errZlimit, errZlimit);
    TH1F* ErrSelectHisto = new TH1F("ErrSelectHisto", ("Residuals " + selectedmult + ";Z_{rec}-Z_{true} (cm);Entries").c_str(), 200, -errZlimit, errZlimit);

    //================================= Analysis ==============================

    int Nevents = inputTree->GetEntries();

    for(int i_event=0; i_event<Nevents; i_event++)
    {
        inputTree->GetEntry(i_event);
        
        FullErrHisto->Fill(recVertex.Zrec - recVertex.Ztrue);
        if(recVertex.mult > multMin && recVertex.mult < multMax)
            ErrSelectHisto->Fill(recVertex.Zrec - recVertex.Ztrue);
    
    }

    if(displayerrfull)  DisplayHisto(FullErrHisto, "Residuals", "residuals_full_mult.png");
    if(displayerrselect) DisplayHisto(ErrSelectHisto, ("Residuals " + selectedmult).c_str(), "residuals_select_mult.png");

    inputFile.Close();
    return 0;
}