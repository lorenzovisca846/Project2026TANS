#include "classes/Cylinder.h"
#include "classes/Point.h"
#include "classes/SimRandom.h"

#include <TRandom3.h>

using namespace std;

void simulation(bool msEnabled, double seed=0, double multMax = 80, double multMin = 20)
{
    double vtxXYsigma = 0.01;
    double vtxZsigma = 5.3;

    double bpR = 3.0, bpL = 27.0, bpW = 0.08;
    double l1R = 4.0, l1L = 27.0, l1W = 0.02;
    double l2R = 7.0, l2L = 27.0;

    Cylinder beamPipe(bpR, bpL, bpW, "Be", msEnabled);
    Cylinder Layer1(l1R, l1L, l1W, "Si", msEnabled);
    Cylinder Layer2(l2R, l2L, 0., "Si", false);            // No need to calculate MS for the outer layer

    typedef struct{
    Point VTXcoord;
    int mult;} VTX;

    delete gRandom;
    SimRandom *simrand = new SimRandom(seed);
    gRandom = simrand;

    //============= VERTEX GENERATION =============//

    VTX vertex;

    vertex.mult = simrand->VMult1(multMin, multMax);
    vertex.VTXcoord = simrand->GausPoint(vtxXYsigma, vtxZsigma);

    // WIP
}