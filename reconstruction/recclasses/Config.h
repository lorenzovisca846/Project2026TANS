#ifndef CONFIG_H
#define CONFIG_H

struct config
{
    // Geometria del detector
    double beamPipeRadius = 3.0;           // cm (raggio interno del beam pipe)
    double beamPipeThickness = 0.08;       // cm (spessore del beam pipe)
    
    double layer1Radius = 4.0;             // cm (raggio del primo layer sensibile)
    double layer1Thickness = 0.02;         // cm (spessore del primo layer)
    
    double layer2Radius = 7.0;             // cm (raggio del secondo layer sensibile)
    double layer2Thickness = 0.02;         // cm (spessore del secondo layer)
    
    double detectorLength = 27.0;          // cm (lunghezza totale del detector lungo z)

    // Parametri del vertice
    double vertexZSigma = 5.3;             // cm (deviazione standard della distribuzione del vertice in z)
    double vertexXYSigma = 0.01;           // cm (deviazione standard della distribuzione del vertice in xy)

    // Parametri di smearing (risoluzione del detector)
    double smearZ = 0.0120;                // cm (120 μm, risoluzione in z)
    double smearRPhi = 0.0030;             // cm (30 μm, risoluzione in rφ)

    // Proprietà dei materiali
    double X0_Be = 35.28;                  // cm (lunghezza di radiazione del Berillio)
    double X0_Si = 9.37;                   // cm (lunghezza di radiazione del Silicio)
    double mspConstant = 0.0136;           // MeV/c (costante per scattering multiplo)

    // Parametri di simulazione
    int nEvents = 1000;                    // Numero di eventi da simulare
    int multiplicityMin = 20;              // Molteplicità minima per evento
    int multiplicityMax = 80;              // Molteplicità massima per evento
    double noiseRateLayer = 5.0;           // % (tasso di rumore per layer)

    // Parametri per l'algoritmo Metropolis
    double metropolisStepSize = 0.01;      // cm (dimensione del passo per Metropolis)
    int metropolisNSteps = 10000;          // Numero di passi per Metropolis

    // Taglio per formare i tracklets
    double deltaPhiCut = 0.1;              // rad (massima differenza in φ per matching hits)

    // Metodo per stampare la configurazione
    void Print()
    {
        cout << "\n=== CONFIGURAZIONE DEL DETECTOR ===" << endl;
        cout << "Geometria:" << endl;
        cout << "  - BeamPipe: R=" << beamPipeRadius << " cm, spessore=" << beamPipeThickness << " cm" << endl;
        cout << "  - Layer1: R=" << layer1Radius << " cm, spessore=" << layer1Thickness << " cm" << endl;
        cout << "  - Layer2: R=" << layer2Radius << " cm, spessore=" << layer2Thickness << " cm" << endl;
        cout << "  - Lunghezza totale: " << detectorLength << " cm" << endl;
        cout << "\nParametri vertice:" << endl;
        cout << "  - σ_z = " << vertexZSigma << " cm" << endl;
        cout << "  - σ_xy = " << vertexXYSigma << " cm" << endl;
        cout << "\nRisoluzione detector:" << endl;
        cout << "  - σ_z = " << smearZ*10000 << " μm" << endl;
        cout << "  - σ_rφ = " << smearRPhi*10000 << " μm" << endl;
        cout << "\nSimulazione:" << endl;
        cout << "  - Eventi: " << nEvents << endl;
        cout << "  - Molteplicità: [" << multiplicityMin << ", " << multiplicityMax << "]" << endl;
        cout << "  - Taglio Δφ: " << deltaPhiCut << " rad" << endl;
        cout << "====================================\n" << endl;
    }
};



#endif // CONFIG_H