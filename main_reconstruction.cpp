#include <TClonesArray.h>
#include <TEnv.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <string>
#include <vector>
#include <algorithm>

#include "classes/Config.h"
#include "classes/MyPoint.h"
#include "classes/Tracklet.h"

vector<Tracklet> FormTracklets(const vector<MyPoint>& hitsLayer1, const vector<MyPoint>& hitsLayer2, Config& config);

double ReconstructVertex(const vector<Tracklet>& tracklets, bool& success, Config& config);
double RunningWindow(const vector<double>& zIntersections, double windowSize, bool& success);

double DeltaPhi(double phi1, double phi2);

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
        bool success;
        int mult;} REC; 
    static REC recVertex;

    TFile outputFile("../outputs/reconstruction_output.root","RECREATE");

    TTree *outputTree = new TTree("Tree_RecOut","Reconstruction Tree");
    outputTree->Branch("Vertex",&recVertex.Ztrue,"Ztrue/D:Zrec:success/O:mult/I");


    //================================= Reconstruction ==============================
    vector<MyPoint> hitsL1;
    vector<MyPoint> hitsL2;
    vector<Tracklet> tracklets;

    int Nevents = inputTree->GetEntries();

    TStopwatch timer;
    timer.Start();

    for(int i_event=0; i_event<Nevents;i_event++)
    {
        if(i_event % 10000 == 0)
            cout << "Reconstructing event " << i_event << "/" << Nevents << endl;


        inputTree->GetEvent(i_event);

        recVertex.mult = trueVertex.mult;
        recVertex.Ztrue = trueVertex.Z;

        int nHitsL1 = ptrhits1->GetEntries();
        hitsL1.clear();
        hitsL1.reserve(nHitsL1);

        for(int j=0; j<nHitsL1; j++)
        {
            hitsL1.push_back(*(MyPoint*)ptrhits1->At(j));
            config.MyRandom()->Smear(hitsL1.back(), config.smearZ, config.smearRPhi, config.detectorLength);
        }

        int nHitsL2 = ptrhits2->GetEntries();
        hitsL2.clear();
        hitsL2.reserve(nHitsL2);

        for(int j=0; j<nHitsL2; j++)
        {
            hitsL2.push_back(*(MyPoint*)ptrhits2->At(j));
            config.MyRandom()->Smear(hitsL2.back(), config.smearZ, config.smearRPhi, config.detectorLength);
        }

        //================================= Tracklet construction =================================
        tracklets.clear();  
        tracklets = FormTracklets(hitsL1, hitsL2, config);

        //================================= Vertex reconstruction =================================
        if(tracklets.size()<1)
        {
            recVertex.Zrec = 0.;
            recVertex.success = false;
            outputTree->Fill();
        }
        else if(tracklets.size() == 1)
        {
            recVertex.Zrec = tracklets[0].z_intersection;
            recVertex.success = true;
            outputTree->Fill();
        }
        else
        {
            recVertex.Zrec = ReconstructVertex(tracklets, recVertex.success, config);
            outputTree->Fill();
        }

        ptrhits1->Clear();
        ptrhits2->Clear();
    }

    timer.Stop(); 
    timer.Print();

    outputFile.Write();

    inputFile.Close();
    outputFile.Close();

    return 0;
}

vector<Tracklet> FormTracklets(const vector<MyPoint>& hitsLayer1, const vector<MyPoint>& hitsLayer2, Config& config)
{
    vector<Tracklet> tracklets;
    vector<pair<double, int>> sortedL1, sortedL2;

    sortedL1.reserve(hitsLayer1.size());
    sortedL2.reserve(hitsLayer2.size());
    
    for (int i=0; i<hitsLayer1.size(); i++)
        sortedL1.emplace_back(hitsLayer1[i].GetPhi(), i);
    
    for (int i=0; i<hitsLayer2.size(); i++)
        sortedL2.emplace_back(hitsLayer2[i].GetPhi(), i);

    sort(sortedL1.begin(), sortedL1.end());
    sort(sortedL2.begin(), sortedL2.end());
    
    for (const auto& [phi1, idx1] : sortedL1)
    {
        int j_start = 0;

        while ((j_start < sortedL2.size()) && (DeltaPhi(phi1, sortedL2[j_start].first) > config.deltaPhiCut))
            j_start++;
        
        for (int j = j_start; j < sortedL2.size() && (DeltaPhi(phi1, sortedL2[j].first) < config.deltaPhiCut); j++)
        {
            const auto& [phi2, idx2] = sortedL2[j];    
            Tracklet tracklet(idx1, idx2);
            tracklet.CalculateTrackletIntersection(hitsLayer1[idx1], hitsLayer2[idx2]);
            tracklets.push_back(tracklet);
        }
    }
    return tracklets;
}

double ReconstructVertex(const vector<Tracklet>& tracklets, bool& success, Config& config)
{
    success = false;

    vector<double> zIntersections;
    zIntersections.reserve(tracklets.size());

    for(int i=0; i<tracklets.size(); i++)
        zIntersections.push_back(tracklets[i].z_intersection);
    sort(zIntersections.begin(), zIntersections.end());

    double Zcenter = RunningWindow(zIntersections, config.runningWindowSize, success);

    return Zcenter;
}

double RunningWindow(const vector<double>& zIntersections, double windowSize, bool& success)
{
    int maxCount = 0;
    int maxIndex = 0;
    int right = 0;
    
    for(int i=0; i<zIntersections.size(); i++)
    {
        double zInit = zIntersections[i];
        while((right < zIntersections.size()) && (zIntersections[right]-zInit < windowSize))
            right++;

        if(right-i > maxCount)
        {
            maxCount = right-i;
            maxIndex = i;
        }
    }

    success = (maxCount > 1);
    if (!success) return 0.;

    double Zrec = 0.;
    for(int i=0;i<maxCount;i++)
        Zrec += zIntersections[maxIndex + i];
    
    return Zrec / maxCount;
}

double DeltaPhi(double phi1, double phi2)
{
    double dPhi = fabs(phi1 - phi2);
    
    if (dPhi > TMath::Pi())
        dPhi = 2.0 * TMath::Pi() - dPhi;

    return dPhi;
}