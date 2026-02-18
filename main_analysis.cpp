#include <iostream>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TEnv.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <string>

using namespace std;

void DisplayResiduals2D(TH2F* histo, string title, string filename);
void DisplayResiduals(TH1F* histo, string title, string filename);
void DisplayResolution(TH1F* histo, string title, string filename);
void DisplayEfficiency(TEfficiency* eff, string title, string filename);
vector<double> FormBinEdges(TH2F* histo);

int main(int argc, char** argv)
{
    string cFile = "fullConfig.txt";
    if (argc > 1) cFile = argv[1];
    string configFile = "../config/" + cFile;

    TEnv *config = new TEnv(configFile.c_str());

    int multMin             = config->GetValue("MultZoom_Min",5);
    int multMax             = config->GetValue("MultZoom_Max",10);

    int multMinGlobal       = config->GetValue("Minimum",1);
    int multMaxGlobal       = config->GetValue("Maximum",69);

    double sigmaZ           = config->GetValue("SigmaZ", 5.3);
    double unifZedges       = config->GetValue("Zedges_uniform", 5.3);

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
    bool displayresZvrt     = config->GetValue("Resolution_Zvert", false);              //WIP

    string gentypes         = config->GetValue("Generation", "ghp");

    char Zdistribution = (gentypes.length() > 0) ? gentypes[0] : 'g';
    char multdistribution = (gentypes.length() > 1) ? gentypes[1] : 'h';

    if(multdistribution == 'h' || multdistribution == 'H') multMinGlobal = 2;

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
    
    int Nevents = inputTree->GetEntries();

    string selectedmult = "(" + to_string(multMin) + " #leq Multiplicity #leq " + to_string(multMax) + ")";
    if(multMin == multMax) selectedmult = "(Multiplicity = " + to_string(multMin) + ")";

    int nbinMG = multMaxGlobal - multMinGlobal + 1;

    TH2F* ErrMultHisto2D    = new TH2F("ErrMultHisto",    "histo;Vertex multiplicity;Z_{rec}-Z_{true}(#mum)", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5, 200, -errZlimit*1e4, errZlimit*1e4);
    TH2F* ErrMultHisto2D_1s = new TH2F("ErrMultHisto_1s", "histo;Vertex multiplicity;Z_{rec}-Z_{true}(#mum)", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5, 200, -errZlimit*1e4, errZlimit*1e4);
    TH2F* ErrMultHisto2D_3s = new TH2F("ErrMultHisto_3s", "histo;Vertex multiplicity;Z_{rec}-Z_{true}(#mum)", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5, 200, -errZlimit*1e4, errZlimit*1e4);

    TH1F* MultEventsHisto  = new TH1F("MultEventsHisto", "histo;Vertex multiplicity;Efficiency", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);
    TH1F* MultSuccessHisto = new TH1F("MultSuccessHisto", "histo;Vertex multiplicity;Efficiency", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);

    TH1F* MultEventsHisto_1s  = new TH1F("MultEventsHisto_1s", "histo;Vertex multiplicity;Efficiency", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);
    TH1F* MultSuccessHisto_1s = new TH1F("MultSuccessHisto_1s", "histo;Vertex multiplicity;Efficiency", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);

    TH1F* MultEventsHisto_3s  = new TH1F("MultEventsHisto_3s", "histo;Vertex multiplicity;Efficiency", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);
    TH1F* MultSuccessHisto_3s = new TH1F("MultSuccessHisto_3s", "histo;Vertex multiplicity;Efficiency", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);

    double binW = 0.5;
    double zMin = -15.;
    double zMax = 15.;
    if(Zdistribution == 'u' || Zdistribution == 'U')
    {
        zMin = -unifZedges;
        zMax = unifZedges;
        cout << "Uniform Z distribution selected: setting Z range to [" << zMin << ", " << zMax << "] cm" << endl;
    }
    int nBinZ = (zMax - zMin) / binW + 1;
    if(nBinZ % 2 == 0) nBinZ++;             // Forzo un bin centrato in zero   
    TH2F* ErrZHisto2D = new TH2F("ErrZHisto2D", "histo;Z_{true}(cm);Z_{rec}-Z_{true}(#mum)", nBinZ, zMin-binW/2., zMax+binW/2., 200, -errZlimit*1e4, errZlimit*1e4);

    // ================================ Loop over events ================================

    for(int i_event=0; i_event<Nevents; i_event++)
    {
        if(i_event%10000 == 0) cout << "Simulating event " << i_event << "/" << Nevents << endl;
        inputTree->GetEntry(i_event);
        MultEventsHisto->Fill(recVertex.mult);

        if(abs(recVertex.Ztrue) < sigmaZ) MultEventsHisto_1s->Fill(recVertex.mult);
        if(abs(recVertex.Ztrue) < 3*sigmaZ) MultEventsHisto_3s->Fill(recVertex.mult);

        if(recVertex.success)
        {
            ErrMultHisto2D->Fill(recVertex.mult, (recVertex.Zrec - recVertex.Ztrue)*1e4);
            if(abs(recVertex.Ztrue) < sigmaZ) ErrMultHisto2D_1s->Fill(recVertex.mult, (recVertex.Zrec - recVertex.Ztrue)*1e4);
            if(abs(recVertex.Ztrue) < 3*sigmaZ) ErrMultHisto2D_3s->Fill(recVertex.mult, (recVertex.Zrec - recVertex.Ztrue)*1e4);

            MultSuccessHisto->Fill(recVertex.mult);
            if(abs(recVertex.Ztrue) < sigmaZ) MultSuccessHisto_1s->Fill(recVertex.mult);
            if(abs(recVertex.Ztrue) < 3*sigmaZ) MultSuccessHisto_3s->Fill(recVertex.mult);

            ErrZHisto2D->Fill(recVertex.Ztrue, (recVertex.Zrec - recVertex.Ztrue)*1e4);
        }
    }

    int binMin = ErrMultHisto2D->GetXaxis()->FindBin(multMin);
    int binMax = ErrMultHisto2D->GetXaxis()->FindBin(multMax);

    TH1F* ErrMultHistoFull      = (TH1F*)ErrMultHisto2D->ProjectionY("ErrMultHistoFull");
    ErrMultHistoFull->SetTitle("Residuals");
    ErrMultHistoFull->GetYaxis()->SetTitle("Entries");
    TH1F* ErrMultHistoSelect    = (TH1F*)ErrMultHisto2D->ProjectionY("ErrMultHistoSelect", binMin, binMax);
    ErrMultHistoSelect->SetTitle(("Residuals " + selectedmult).c_str());
    ErrMultHistoSelect->GetYaxis()->SetTitle("Entries");

    // ================================ Residuals vs multiplicity ================================

    if(displayerrfull)   DisplayResiduals(ErrMultHistoFull, "Residuals", "residuals_full.png");
    if(displayerrselect) DisplayResiduals(ErrMultHistoSelect, ("Residuals " + selectedmult).c_str(), "residuals_selected.png");


    if(displayerrfull)   DisplayResiduals2D(ErrMultHisto2D, "Residuals vs Vertex multiplicity", "residuals2D_vs_mult.png");
    if(displayerrfull)   DisplayResiduals2D(ErrMultHisto2D_1s, "Residuals vs Vertex multiplicity #left(#left|Z_{true}#right|<#sigma#right)", "residuals2D_vs_mult_1sigma.png");
    if(displayerrfull)   DisplayResiduals2D(ErrMultHisto2D_3s, "Residuals vs Vertex multiplicity #left(#left|Z_{true}#right|<3#sigma#right)", "residuals2D_vs_mult_3sigma.png");


    // ================================ Resolution vs multiplicity ================================


    TH1F* ResMultHisto    = new TH1F("ResMultHisto",    "histo;Vertex multiplicity;Resolution (#mum)", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);
    TH1F* ResMultHisto_1s = new TH1F("ResMultHisto_1s", "histo;Vertex multiplicity;Resolution (#mum)", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);
    TH1F* ResMultHisto_3s = new TH1F("ResMultHisto_3s", "histo;Vertex multiplicity;Resolution (#mum)", nbinMG, multMinGlobal-0.5, multMaxGlobal+0.5);


    TH1D* slice;
    for(int i_mult = multMinGlobal; i_mult <= multMaxGlobal; i_mult++)
    {
        int targetBin = ErrMultHisto2D->GetXaxis()->FindBin(i_mult);
        slice = ErrMultHisto2D->ProjectionY("slice", targetBin, targetBin);
        if(slice->GetEntries() > 0)
        {
            ResMultHisto->SetBinContent(targetBin, slice->GetRMS());
            ResMultHisto->SetBinError(targetBin, slice->GetRMSError());
        }
        delete slice;

        slice = ErrMultHisto2D_1s->ProjectionY("slice", targetBin, targetBin);
        if(slice->GetEntries() > 0)
        {
            ResMultHisto_1s->SetBinContent(targetBin, slice->GetRMS());
            ResMultHisto_1s->SetBinError(targetBin, slice->GetRMSError());
        }
        delete slice;

        slice = ErrMultHisto2D_3s->ProjectionY("slice", targetBin, targetBin);
        if(slice->GetEntries() > 0)
        {
            ResMultHisto_3s->SetBinContent(targetBin, slice->GetRMS());
            ResMultHisto_3s->SetBinError(targetBin, slice->GetRMSError());
        }
        delete slice;
    }

    if(displayresfull)     DisplayResolution(ResMultHisto, "Resolution vs Vertex multiplicity", "resolution_vs_mult.png");
    if(displayres1sigma)   DisplayResolution(ResMultHisto_1s, "Resolution vs Vertex multiplicity #left(#left|Z_{true}#right|<#sigma#right)", "resolution_vs_mult_1sigma.png");
    if(displayres3sigma)   DisplayResolution(ResMultHisto_3s, "Resolution vs Vertex multiplicity #left(#left|Z_{true}#right|<3#sigma#right)", "resolution_vs_mult_3sigma.png");

    // ================================ Efficiency vs multiplicity ================================

    TEfficiency* effMultHisto = new TEfficiency(*MultSuccessHisto, *MultEventsHisto);
    TEfficiency* effMultHisto_1s = new TEfficiency(*MultSuccessHisto_1s, *MultEventsHisto_1s);
    TEfficiency* effMultHisto_3s = new TEfficiency(*MultSuccessHisto_3s, *MultEventsHisto_3s);

    if(displayefffull)     DisplayEfficiency(effMultHisto, "Efficiency vs Vertex multiplicity", "efficiency_vs_mult.png");
    if(displayeff1sigma)   DisplayEfficiency(effMultHisto_1s, "Efficiency vs Vertex multiplicity #left(#left|Z_{true}#right|<#sigma#right)", "efficiency_vs_mult_1sigma.png");
    if(displayeff3sigma)   DisplayEfficiency(effMultHisto_3s, "Efficiency vs Vertex multiplicity #left(#left|Z_{true}#right|<3#sigma#right)", "efficiency_vs_mult_3sigma.png");

    // ================================ Efficiency vs Zvert ================================

    vector<double> BinEdgesZ = FormBinEdges(ErrZHisto2D);
    int nBinsZ = BinEdgesZ.size() - 1;

    double BinEdgesArray[nBinsZ+1];
    for(int i=0; i<nBinsZ+1; i++)     BinEdgesArray[i] = BinEdgesZ[i];
    
    TH1F* ZEventsHisto  = new TH1F("ZEventsHisto", "histo;Z_{true}(cm);Efficiency", nBinsZ, BinEdgesArray);
    TH1F* ZSuccessHisto = new TH1F("ZSuccessHisto", "histo;Z_{true}(cm);Efficiency", nBinsZ, BinEdgesArray);

    for(int i_event=0; i_event<Nevents; i_event++)
    {
        inputTree->GetEntry(i_event);
        ZEventsHisto->Fill(recVertex.Ztrue);

        if(recVertex.success)
            ZSuccessHisto->Fill(recVertex.Ztrue);
    }

    TEfficiency* effZ = new TEfficiency(*ZSuccessHisto, *ZEventsHisto);

    if(displayeffZvrt) DisplayEfficiency(effZ, "Efficiency vs Z_{true}", "efficiency_vs_Zvert.png");


    // ================================ Resolution vs Zvert ================================

    TH1F* ZResHisto     = new TH1F("ZResHisto", "histo;Z_{true}(cm);Resolution (#mum)", nBinsZ, BinEdgesArray);

    for(int i=1;i<=nBinsZ;i++)
    {
        int bin1 = ErrZHisto2D->GetXaxis()->FindBin(BinEdgesArray[i-1] + 0.0001);
        int bin2 = ErrZHisto2D->GetXaxis()->FindBin(BinEdgesArray[i] - 0.0001);
        slice = ErrZHisto2D->ProjectionY("slice", bin1, bin2);
        if(slice->GetEntries() > 0)
        {
            ZResHisto->SetBinContent(i, slice->GetRMS());
            ZResHisto->SetBinError(i, slice->GetRMSError());
        }
        delete slice;
    }

    if(displayresZvrt) DisplayResiduals2D(ErrZHisto2D, "Residuals vs Z_{true}", "residuals2D_vs_Zvert.png");
    if(displayresZvrt) DisplayResolution(ZResHisto, "Resolution vs Z_{true}", "resolution_vs_Zvert.png");

    inputFile.Close();
    return 0;
}

void DisplayResiduals2D(TH2F* histo, string title, string filename)
{
    TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 1000, 600);
    histo->SetTitle(title.c_str());
    histo->Draw("COLZ");

    gStyle->SetOptStat("e");
    TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
    st->SetY1NDC(0.92);
    st->SetY2NDC(0.98);

    c->SaveAs(("../outputs/plots/" + filename).c_str());
}

void DisplayResiduals(TH1F* histo, string title, string filename)
{
    TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 1000, 600);
    histo->SetMarkerColor(kBlue);
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(0.4);
    histo->SetTitle(title.c_str());
    histo->Draw("E1");
    histo->Fit("gaus", "Q");

    gStyle->SetEndErrorSize(0);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(1100);
    TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
    st->SetY1NDC(0.88);
    st->SetY2NDC(0.98);
    st->SetX1NDC(0.75);
    st->SetX2NDC(0.98);

    c->SaveAs(("../outputs/plots/" + filename).c_str());
}

void DisplayResolution(TH1F* histo, string title, string filename)
{
    TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 1000, 600);

    histo->GetYaxis()->SetRange(histo->GetMinimum()*0.8, histo->GetMaximum()*1.2);
    histo->SetMarkerColor(kBlue);
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(0.4);

    if(title.find("multiplicity") == string::npos)
        histo->SetMarkerSize(0.8);

    histo->SetTitle(title.c_str());
    histo->Draw("E1");

    gStyle->SetEndErrorSize(0);
    gStyle->SetOptStat(0);

    c->SaveAs(("../outputs/plots/" + filename).c_str());
}

void DisplayEfficiency(TEfficiency* histo, string title, string filename)
{
    TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 1000, 600);
    histo->SetMarkerColor(kBlue);
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(0.4);

    if(title.find("multiplicity") == string::npos)
        histo->SetMarkerSize(0.8);

    histo->SetTitle(title.c_str());
    histo->Draw("E1");

    gStyle->SetEndErrorSize(0);
    gStyle->SetOptStat(0);

    c->SaveAs(("../outputs/plots/" + filename).c_str());
}

vector<double> FormBinEdges(TH2F* histo)
{
    vector<double> binEdges;
    TAxis* xAxis = histo->GetXaxis();

    int nBins = xAxis->GetNbins();
    int centralBin = (nBins / 2) + 1;

    binEdges.push_back(xAxis->GetBinLowEdge(centralBin));
    binEdges.push_back(xAxis->GetBinUpEdge(centralBin));

    int edge = centralBin + 1;
    while(edge<=nBins)
    {
        int j=0;
        while(histo->Integral(edge, edge+j, 1, histo->GetNbinsY()) < 1000 && edge+j <= nBins)
            j++;
        binEdges.push_back(xAxis->GetBinUpEdge(edge+j));
        edge += (j+1);
    }

    edge = centralBin-1;
    while(edge >= 1)
    {
        int j=0;
        while(histo->Integral(edge-j, edge, 1, histo->GetNbinsY()) < 1000 && edge-j >= 1)
            j++;
        binEdges.push_back(xAxis->GetBinLowEdge(edge-j));
        edge -= (j+1);
    }

    sort(binEdges.begin(), binEdges.end());

    return binEdges;
}