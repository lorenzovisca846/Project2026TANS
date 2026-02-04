#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TEnv.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

#include "../simulation/simclasses/MyPoint.h"

using namespace std;

typedef struct{
    double X, Y, Z;
    int mult;} VTX;

void geometry();
void event();

void eventdisplay()
{
    geometry();
    event();
}

void geometry()
{
    TEnv *config = new TEnv("../config/simDisplay.txt");

    double bpR          = config->GetValue("BeamPipeRadius", 3.0);
    double bpL          = config->GetValue("Length", 27.0);
    double bpW          = config->GetValue("BeamPipeWidth", 0.08);

    double l1R          = config->GetValue("Layer1Radius", 4.0);
    double l1L          = config->GetValue("Length", 27.0);
    double l1W          = config->GetValue("LayerWidth", 0.02);

    double l2R          = config->GetValue("Layer2Radius", 7.0);
    double l2L          = config->GetValue("Length", 27.0);
    double l2W          = config->GetValue("LayerWidth", 0.02);

    delete config;

    TGeoManager *geom = new TGeoManager("DetectorGeom", "Geometria Tracker 3 Layer");
    
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium *medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    
    double world_size = l2L * 1.1; 
    TGeoVolume *top = geom->MakeBox("TOP", medVacuum, world_size, world_size, world_size);
    geom->SetTopVolume(top);
    
    TGeoVolume *beampipe = geom->MakeTube("BEAM_PIPE", medVacuum, bpR, bpR + bpW, bpL *1.4 / 2.0);
    beampipe->SetLineColor(kGray);
    top->AddNode(beampipe, 1);
    
    TGeoVolume *layer1 = geom->MakeTube("LAYER1", medVacuum, l1R, l1R + l1W, l1L / 2.0);
    layer1->SetLineColor(kYellow);
    top->AddNode(layer1, 1);
    
    TGeoVolume *layer2 = geom->MakeTube("LAYER2", medVacuum, l2R, l2R + l2W, l2L / 2.0);
    layer2->SetLineColor(kOrange);
    top->AddNode(layer2, 1);

    geom->CloseGeometry();

    top->Draw("ogl");
}

void event() 
{

    TFile *InputFile = TFile::Open("../sim_display.root","READ");
    TTree *tree = (TTree*)InputFile->Get("Tree_SimOut");

    double VTXx = 0., VTXy = 0., VTXz = 0.;
    TClonesArray* hits_bp = nullptr;
    TClonesArray* hits_l1 = nullptr;
    TClonesArray* hits_l2 = nullptr;
    VTX vertex;

    tree->SetBranchAddress("Hits_L1", &hits_l1);
    tree->SetBranchAddress("Hits_L2", &hits_l2);
    tree->SetBranchAddress("Hits_BP", &hits_bp);
    tree->SetBranchAddress("Vertex", &vertex.X);

    tree->GetEntry(0);

    map<int, vector<MyPoint*>> track_map;
    vector<MyPoint*> noise_points;

    auto sort_hits = [&](TClonesArray* array)
    {
        if (!array) return;
        for (int i = 0; i < array->GetEntriesFast(); ++i)
        {
            MyPoint* p = (MyPoint*)array->At(i);
            int tid = p->GetTrackID();
            if (tid != -1) 
                track_map[tid].push_back(p);
            else
                noise_points.push_back(p);
        }
    };

    sort_hits(hits_bp);
    sort_hits(hits_l1);
    sort_hits(hits_l2);

    for (auto const& [tid, hits] : track_map)
    {
        int n_points = hits.size() + 1;
        TPolyLine3D *line = new TPolyLine3D(n_points);
        TPolyMarker3D *markers = new TPolyMarker3D(hits.size());
        
        line->SetPoint(0, vertex.X, vertex.Y, vertex.Z);
        
        
        for (size_t i = 0; i < hits.size(); ++i)
        {
            line->SetPoint(i + 1, hits[i]->GetX(), hits[i]->GetY(), hits[i]->GetZ());
            if(i>0)
                markers->SetPoint(i, hits[i]->GetX(), hits[i]->GetY(), hits[i]->GetZ());
        }

        markers->SetMarkerStyle(20);
        markers->SetMarkerSize(0.2);
        markers->SetMarkerColor(kSpring);

        if (hits.size() < 3)
            line->SetLineColor(kBlue);
        else
            line->SetLineColor(kSpring);
        line->SetLineWidth(3);
        line->Draw("same");
        markers->Draw("same");
    }

    TPolyMarker3D *noise_markers = new TPolyMarker3D(noise_points.size());
    for (size_t i = 0; i < noise_points.size(); ++i)
        noise_markers->SetPoint(i,noise_points[i]->GetX(), noise_points[i]->GetY(), noise_points[i]->GetZ());

    noise_markers->SetMarkerStyle(20);
    noise_markers->SetMarkerSize(0.2);
    noise_markers->SetMarkerColor(kRed);
    noise_markers->Draw("same");
}