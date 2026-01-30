#include <Riostream.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TStopwatch.h>

#include "classes/Cylinder.h"
#include "classes/SimRandom.h"
#include "classes/Particle.h"

//Nota: cambiare questi define con una funzione apposita che linka la funzione chiamata ai puntatori, l'abbiamo fatto a lezione
#define MULT_METHOD 1
#define VTX_METHOD 1

using namespace std;

void Transport(Particle* part, Cylinder& layer, TClonesArray& hits, bool detector, bool msEnabled);

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
                  INPUT FILE WITH DISTRIBUTIONS

                  WIP (file does not exist yet)
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
    ====================================================
            DECLARATION OF OUTPUT FILE AND TTREE
    ====================================================
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

    TStopwatch timer;
    timer.Start();

    Particle *ptrPart = new Particle(simrand);

    for(unsigned int i=0; i<Nevents; i++)
    {
        /*
        =================================
                VERTEX GENERATION
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

            ptrPart->Init(vertex.X, vertex.Y, vertex.Z, 1.0, 0.8, 1); // beta=1, p=0.8GeV/c, Q=1e
            Transport(ptrPart, beamPipe, hits1, false, msEnabled);
            Transport(ptrPart, Layer1, hits1, true, msEnabled);
            Transport(ptrPart, Layer2, hits2, true, false);
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

    timer.Stop(); 
    timer.Print();

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

void Transport(Particle* part, Cylinder& layer, TClonesArray& hits, bool detector, bool msEnabled)
{
    // WIP
    // 1) Chiamare la funzione di propagazione della particella fino al raggio interno del layer (da scrivere)
    // 2) Verificare se la particella interseca il layer
    // 3) Se è in accettanza e il layer è un rivelatore, salvare la hit nel TClonesArray. Capire se in cartesiane o cilindriche
    // 4) Se è in accettanza e il MS è acceso, chiamare la funzione di ms della particella (da scrivere)
}