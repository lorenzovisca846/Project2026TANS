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
#include <TCanvas.h>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

#include "classes/MyPoint.h"

using namespace std;

typedef struct{
    double X, Y, Z;
    int mult;} VTX;


void sethits(map<int, vector<MyPoint*>>&, vector<MyPoint*>&, TClonesArray*);

void eventdisplay()
{
    TCanvas *canvas = new TCanvas();

    TEnv *config = new TEnv("config/simDisplay.txt");

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
    
    double world_size = bpL* 1.4; 
    TGeoVolume *top = geom->MakeBox("TOP", medVacuum, world_size, world_size, world_size);
    geom->SetTopVolume(top);
    
    TGeoVolume *beampipe = geom->MakeTube("BEAM_PIPE", medVacuum, bpR, bpR + bpW, bpL * 1.4 / 2.0);
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

    TFile *InputFile = TFile::Open("outputs/sim_display.root","READ");
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

    map<int, vector<MyPoint*>> hits_map;
    vector<MyPoint*> noise_points;

    sethits(hits_map, noise_points, hits_bp);
    sethits(hits_map, noise_points, hits_l1);
    sethits(hits_map, noise_points, hits_l2);

    TPolyMarker3D *noise_markers = new TPolyMarker3D(noise_points.size());
    for(int i=0; i<noise_points.size(); i++)
        noise_markers->SetPoint(i,noise_points[i]->GetX(), noise_points[i]->GetY(), noise_points[i]->GetZ());

    noise_markers->SetMarkerStyle(20);
    noise_markers->SetMarkerSize(0.2);
    noise_markers->SetMarkerColor(kRed);
    noise_markers->Draw("same");


    TPolyMarker3D *markers;
    TPolyLine3D *line;

    double r;

    for(int i=0;i<vertex.mult; i++)
    {
        int Npoints = hits_map[i].size() + 1;
        line = new TPolyLine3D(Npoints);
        markers = new TPolyMarker3D(Npoints);

        line->SetPoint(0, vertex.X, vertex.Y, vertex.Z);

        double Nhits = hits_map[i].size();

        int nLayers = 0;
        for(int j=0;j<Nhits; j++)
        {
            line->SetPoint(j+1, hits_map[i][j]->GetX(), hits_map[i][j]->GetY(), hits_map[i][j]->GetZ());
            
            r = sqrt(pow(hits_map[i][j]->GetX(),2) + pow(hits_map[i][j]->GetY(),2));
            if(r>bpR + 0.001)
            {
                markers->SetPoint(j, hits_map[i][j]->GetX(), hits_map[i][j]->GetY(), hits_map[i][j]->GetZ());
                nLayers++;
            }
            
            markers->SetMarkerStyle(20);
            markers->SetMarkerSize(0.2);
            markers->SetMarkerColor(kSpring);
        }

        markers->Draw("same");
        line->SetLineColor(kSpring);
        if((Nhits < 2) || (Nhits == 2 && nLayers == 1))
            line->SetLineColor(kBlue);
        line->SetLineWidth(3);
        line->Draw("same");
    }
}

void sethits(map<int, vector<MyPoint*>>& hits_map, vector<MyPoint*>& noise_points, TClonesArray* hits)
{
    int Nentries = hits->GetEntries();
    MyPoint* p;
    int trackid;

    for(int i=0; i<Nentries;i++)
    {
        p = (MyPoint*)hits->At(i);
        trackid = p->GetTrackID();
        if(trackid == -1)
            noise_points.push_back(p);
        else
            hits_map[trackid].push_back(p);
    }
}