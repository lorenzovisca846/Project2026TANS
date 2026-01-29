#include <Riostream.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TStopwatch.h>

#include "classes/Cylinder.h"
#include "classes/Point.h"
#include "classes/SimRandom.h"

//Nota: cambiare questi define con una funzione apposita che linka la funzione chiamata ai puntatori, l'abbiamo fatto a lezione
#define MULT_METHOD 1
#define VTX_METHOD 1

using namespace std;

void simulation(double Nevents = 1000, bool msEnabled, unsigned int seed = 0)
{
    const int multMax = 80;
    const int multMin = 20;

    const int noiseMax = 20;                        // max noise hits per layer

    const double vtxXYsigma = 0.01;
    const double vtxZsigma = 5.3;

    const double bpR = 3.0, bpL = 27.0, bpW = 0.08;
    const double l1R = 4.0, l1L = 27.0, l1W = 0.02;
    const double l2R = 7.0, l2L = 27.0;

    Cylinder beamPipe(bpR, bpL, bpW, "Be");
    Cylinder Layer1(l1R, l1L, l1W, "Si");
    Cylinder Layer2(l2R, l2L, 0., "Si");            // No need to calculate MS for the outer layer

    typedef struct{
        double X, Y, Z;
        int mult;} VTX;
    static VTX vertex;

    /*
    ============================================================
        APRIRE FILE E ISTOGRAMMI PER DISTRIBUZIONI MULT, ETA
        Nota: passare subito gli istogrammi a simRandom, in
        modo da non doverli aprire ad ogni evento

                            WIP
    ============================================================
    */

    TFile *inputFile = new TFile("inputDistributions.root","READ");
    TH1F *multHist= (TH1F*)inputFile->Get("multHist"); 
    TH1F *etaHist= (TH1F*)inputFile->Get("etaHist"); 

    delete gRandom;
    SimRandom *simrand = new SimRandom(seed, multHist, etaHist);
    gRandom = simrand;

    inputFile->Close();
    delete inputFile;

    /*
    ============================================================
        INSERIRE DICHIARAZIONE DEL TREE E DELLE BRANCHES

                            WIP
    ============================================================
    */
    TFile hfile("htree.root","RECREATE");
    TTree *tree = new TTree("Tree","Vertex-Hits TTree");

    int arrdim = Nevents * (multMax + noiseMax + 5);            //5 as safety margin

    TClonesArray *ptrhits1 = new TClonesArray("Point",arrdim);
    TClonesArray &hits1 = *ptrhits1;

    TClonesArray *ptrhits2 = new TClonesArray("Point",arrdim);
    TClonesArray &hits2 = *ptrhits2;

    tree->Branch("Vertex",&vertex.X,"X/D:Y:Z:mult/I");
    tree->Branch("Hits_L1",&ptrhits1);
    tree->Branch("Hits_L2",&ptrhits2);


    for(unsigned int i=0; i<Nevents; i++)
    {
        /*
        =================================
                VERTEX GENERATION

                    WIP
        =================================
        */

        #if MULT_METHOD==1
            vertex.mult = multMax;
        #elif MULT_METHOD==2
            vertex.mult = simrand->VMult1(multMin, multMax);
        #else
            do{vertex.mult = simrand->VMult2(multMin, multMax);} while(vertex.mult < multMin || vertex.mult > multMax);
        #endif

        #if VTX_METHOD==1
            simrand->GausPoint(vertex.X, vertex.Y, vertex.Z, vtxXYsigma, vtxZsigma);
        #elif VTX_METHOD==2
            simrand->OriginPoint(vertex.X, vertex.Y, vertex.Z);
        #else
            simrand->UnifPoint(vertex.X, vertex.Y, vertex.Z, vtxXYsigma, vtxZsigma);
        #endif

        for(unsigned int j=0; j<vertex.mult; j++)
        {
            /*
            ============================
                PARTICLE PROPAGATION

                        WIP
            ============================
            */
        }

        /*
        ========================
            NOISE GENERATION

                WIP
        ========================
        */


        /*
        ======================================================
            FILLING THE TREE WITH VERTEX, HITS L1, HITS L2
    
                            WIP
        ======================================================
        */

        tree->Fill();
        ptrhits1->Clear();
        ptrhits2->Clear();
    }

    /*
    =====================================
        CLEANING UP AND CLOSING FILES

                    WIP
    =====================================
    */


    hfile.Write();
    hfile.Close();
    delete simrand;
    delete ptrhits1;
    delete ptrhits2;

}