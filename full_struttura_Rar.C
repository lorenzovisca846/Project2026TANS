#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <TEllipse.h>
#include <TLine.h>

using namespace std;

// ============================================================
// STRUTTURA DI CONFIGURAZIONE
// ============================================================ 

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

// ============================================================
// CLASSI PER I DATI DELLA SIMULAZIONE
// ============================================================

/**
 * @class Particle
 * @brief Rappresenta una particella generata nell'evento
 * 
 * Contiene tutte le informazioni cinematiche e di identificazione
 * di una singola particella.
 */
class Particle
{
public: 
    int trackID;                // ID univoco della traccia (>=0 per segnale, -1 per rumore)
    TVector3 position;          // Posizione 3D (x, y, z)
    TVector3 direction;         // Direzione 3D normalizzata (versore)
    double phi;                 // Angolo azimutale φ [rad]
    double theta;               // Angolo polare θ [rad]
    double p;                   // Impulso totale [GeV/c]
    double pt;                  // Impulso trasverso [GeV/c]
    int charge;                 // Carica elettrica (±1)
    double beta;                // Velocità relativistica β = v/c

    // Costruttore di default
    Particle() : trackID(-1), theta(0), phi(0), p(0), pt(0), charge(0), beta(1.0) {}
};

/**
 * @class Hit
 * @brief Rappresenta un hit (segnale) nel detector
 * 
 * Contiene le informazioni di un hit misurato, sia segnale che rumore.
 */
class Hit
{
public:
    int layer;                  // Layer del detector (1 o 2)
    double x;                   // Coordinata x misurata [cm]
    double y;                   // Coordinata y misurata [cm]
    double z;                   // Coordinata z misurata [cm]
    double r;                   // Raggio r = √(x²+y²) [cm]
    double phi;                 // Angolo azimutale φ = atan2(y,x) [rad]
    double z_smeared;           // Coordinata z con smearing [cm]
    double phi_smeared;         // Angolo φ con smearing [rad]
    bool isSignal;              // true = hit di segnale, false = hit di rumore
    int trackID;                // ID della traccia associata (-1 se rumore)

    // Costruttore di default
    Hit() : layer(0), x(0), y(0), z(0), r(0), phi(0), 
            z_smeared(0), phi_smeared(0), isSignal(true), trackID(-1) {}
};

/**
 * @class Tracklet
 * @brief Coppia di hit matched tra due layer
 * 
 * Rappresenta un segmento di traccia formato da un hit nel layer1
 * e un hit nel layer2 che probabilmente appartengono alla stessa particella.
 */
class Tracklet
{
public:
    int hit1_idx;               // Indice dell'hit nel vettore hitsLayer1
    int hit2_idx;               // Indice dell'hit nel vettore hitsLayer2
    double z_intersection;      // Coordinata z dell'intersezione a r=0 [cm]
    double slope;               // Pendenza dz/dr del tracklet
    double chi2;                // χ² del fit lineare
    double sigma_z;             // Errore stimato su z_intersection [cm]
    bool valid;                 // true se il tracklet è valido

    // Costruttore di default
    Tracklet() : hit1_idx(-1), hit2_idx(-1), z_intersection(0), 
                 slope(0), chi2(0), sigma_z(0.01), valid(false) {}
    
    /**
     * @brief Calcola il χ² del tracklet
     * @return Valore del χ²
     * 
     * Il χ² misura la bontà del fit lineare tra i due hit.
     * Valori bassi indicano un buon allineamento.
     */
    double CalculateChi2(double smearZ) const
    {
        // Per un tracklet con 2 punti e 2 parametri (z0, slope),
        // il χ² ha 0 gradi di libertà. Questo calcolo è principalmente
        // per scopi di monitoraggio e selezione.
        if (!valid) return 1e6;
        
        // Residuo normalizzato: (z_mis - z_pred)/σ
        // In questo caso semplice, assumiamo che il fit sia perfetto
        // quindi il residuo è zero. Un'implementazione reale userebbe
        // i residui rispetto alla retta fittata.
        return 0.0;
    }
};

/**
 * @class Event
 * @brief Contiene tutte le informazioni di un evento simulato
 */
class Event
{
public:
    int eventID;                        // ID univoco dell'evento
    int multiplicity;                   // Numero di particelle generate
    double trueVertexZ;                 // Coordinata z vera del vertice [cm]
    double recoVertexZ;                 // Coordinata z ricostruita [cm]
    double recoVertexZMetropolis;       // Coordinata z ricostruita con Metropolis [cm]
    bool vertexFound;                   // true se il vertice è stato ricostruito
    int nTracklets;                     // Numero di tracklets formati

    // Costruttore di default
    Event() : eventID(0), multiplicity(0), trueVertexZ(0), 
              recoVertexZ(0), recoVertexZMetropolis(0), 
              vertexFound(false), nTracklets(0) {}
};

// ============================================================
// CLASSE DETECTOR
// ============================================================

/**
 * @class Detector
 * @brief Gestisce la geometria e l'accettanza del detector
 * 
 * Contiene tutti i parametri geometrici e fornisce metodi
 * per verificare se un punto è all'interno dell'accettanza.
 */
class Detector
{
private:
    config fConfig;  // Configurazione del detector

public:
    /**
     * @brief Costruttore della classe Detector
     * @param config Configurazione del detector
     */
    Detector(const config& cfg) : fConfig(cfg) {} 

    /**
     * @brief Verifica se un punto è nell'accettanza del layer
     * @param punto Posizione 3D da verificare
     * @param layer Numero del layer (1 o 2)
     * @return true se il punto è nell'accettanza, false altrimenti
     * 
     * Un punto è accettato se:
     * 1. La sua distanza radiale è vicina al raggio del layer (entro lo spessore)
     * 2. La sua coordinata z è entro la metà della lunghezza del detector
     */
    bool Accettanza(const TVector3& punto, int layer) const
    {
        double raggio_target;
        double spessore;
        
        // Seleziona il raggio e lo spessore in base al layer
        if (layer == 1) {
            raggio_target = fConfig.layer1Radius;
            spessore = fConfig.layer1Thickness;
        } else if (layer == 2) {
            raggio_target = fConfig.layer2Radius;
            spessore = fConfig.layer2Thickness;
        } else {
            // Layer non valido
            return false;
        }
        
        // Calcola la distanza radiale dal centro
        double r = punto.Perp();
        
        // Verifica se il punto è vicino alla superficie del layer
        bool dentro_spessore = fabs(r - raggio_target) < spessore/2.0;
        
        // Verifica se il punto è entro la lunghezza del detector
        bool dentro_lunghezza = fabs(punto.Z()) < fConfig.detectorLength/2.0;
        
        return dentro_spessore && dentro_lunghezza;
    }

    /**
     * @brief Disegna la geometria del detector
     * @param pad Pannello di ROOT su cui disegnare
     * 
     * Disegna la geometria del detector in vista trasversa (xy).
     * I layer sono rappresentati come cerchi concentrici.
     */
    void Geometria(TVirtualPad* pad) const
    {
        pad->cd();  // Seleziona il pannello per il disegno
        
        // ----- BEAM PIPE -----
        // Cerchio interno del beam pipe
        TEllipse* pipeInner = new TEllipse(0, 0, fConfig.beamPipeRadius);
        pipeInner->SetFillStyle(0);      // Riempimento vuoto
        pipeInner->SetLineColor(kBlue);  // Colore blu
        pipeInner->SetLineWidth(2);      // Spessore linea
        pipeInner->Draw();               // Disegna
        
        // Cerchio esterno del beam pipe
        TEllipse* pipeOuter = new TEllipse(0, 0, 
                                           fConfig.beamPipeRadius + fConfig.beamPipeThickness);
        pipeOuter->SetFillStyle(0);      // Riempimento vuoto
        pipeOuter->SetLineColor(kBlue);  // Colore blu
        pipeOuter->SetLineWidth(2);      // Spessore linea
        pipeOuter->Draw("same");         // Disegna sullo stesso pannello
        
        // ----- LAYER 1 -----
        // Cerchio interno del layer 1
        TEllipse* layer1Inner = new TEllipse(0, 0, fConfig.layer1Radius);
        layer1Inner->SetFillStyle(0);    // Riempimento vuoto
        layer1Inner->SetLineColor(kRed); // Colore rosso
        layer1Inner->SetLineWidth(2);    // Spessore linea
        layer1Inner->Draw("same");       // Disegna sullo stesso pannello
        
        // Cerchio esterno del layer 1
        TEllipse* layer1Outer = new TEllipse(0, 0, 
                                             fConfig.layer1Radius + fConfig.layer1Thickness);
        layer1Outer->SetFillStyle(0);    // Riempimento vuoto
        layer1Outer->SetLineColor(kRed); // Colore rosso
        layer1Outer->SetLineWidth(2);    // Spessore linea
        layer1Outer->Draw("same");       // Disegna sullo stesso pannello
        
        // ----- LAYER 2 -----
        // Cerchio interno del layer 2
        TEllipse* layer2Inner = new TEllipse(0, 0, fConfig.layer2Radius);
        layer2Inner->SetFillStyle(0);      // Riempimento vuoto
        layer2Inner->SetLineColor(kGreen); // Colore verde
        layer2Inner->SetLineWidth(2);      // Spessore linea
        layer2Inner->Draw("same");         // Disegna sullo stesso pannello
        
        // Cerchio esterno del layer 2
        TEllipse* layer2Outer = new TEllipse(0, 0, 
                                             fConfig.layer2Radius + fConfig.layer2Thickness);
        layer2Outer->SetFillStyle(0);      // Riempimento vuoto
        layer2Outer->SetLineColor(kGreen); // Colore verde
        layer2Outer->SetLineWidth(2);      // Spessore linea
        layer2Outer->Draw("same");         // Disegna sullo stesso pannello
        
        // ----- ASSE Z (per riferimento) -----
        TLine* zAxis = new TLine(-fConfig.detectorLength/2, 0, 
                                  fConfig.detectorLength/2, 0);
        zAxis->SetLineColor(kBlack);  // Colore nero
        zAxis->SetLineStyle(2);       // Linea tratteggiata
        zAxis->Draw("same");          // Disegna sullo stesso pannello
    }
};

// ============================================================
// CLASSE MULTIPLE SCATTERING
// ============================================================

/**
 * @class MultipleScattering
 * @brief Gestisce lo scattering multiplo delle particelle
 * 
 * Simula l'effetto dello scattering multiplo quando una particella
 * attraversa un materiale. Lo scattering è modellato come una deviazione
 * gaussiana angolare.
 */
class MultipleScattering
{
private:
    config fConfig;    // Configurazione
    TRandom3 fRandom;  // Generatore di numeri casuali

public:
    /**
     * @brief Costruttore della classe MultipleScattering
     * @param config Configurazione del detector
     */
    MultipleScattering(const config& cfg) : fConfig(cfg), fRandom(0) {}

    /**
     * @brief Calcola l'angolo di scattering RMS
     * @param p Impulso della particella [GeV/c]
     * @param beta Velocità relativistica β = v/c
     * @param thickness Spessore del materiale attraversato [cm]
     * @param charge Carica della particella [e]
     * @param material Nome del materiale ("Be" o "Si")
     * @return Angolo RMS di scattering [rad]
     * 
     * Usa la formula di Highland per lo scattering multiplo:
     * θ₀ = (0.0136 GeV/c * |Z|) / (βp) * √(t/X₀) * [1 + 0.038 ln(t/X₀)]
     * dove t è lo spessore e X₀ la lunghezza di radiazione.
     */
    double CalcoloTheta0(double p, double beta, double thickness, 
                         int charge, const string& material) const
    {
        // Seleziona la lunghezza di radiazione in base al materiale
        double X0;
        if (material == "Be") {
            X0 = fConfig.X0_Be;
        } else if (material == "Si") {
            X0 = fConfig.X0_Si;
        } else {
            // Materiale sconosciuto, usa un valore di default
            X0 = 10.0;
        }
        
        // Frazione di lunghezza di radiazione attraversata
        double t_over_X0 = thickness / X0;
        
        // Formula di Highland per l'angolo RMS di scattering
        double theta0 = (fConfig.mspConstant * fabs(charge) / (beta * p)) 
                       * sqrt(t_over_X0) 
                       * (1.0 + 0.038 * log(t_over_X0));
        
        return theta0;
    }

    /**
     * @brief Applica lo scattering multiplo a una particella
     * @param particle Particella da cui applicare lo scattering
     * @param thickness Spessore del materiale attraversato [cm]
     * @param material Nome del materiale ("Be" o "Si")
     * 
     * Genera deviazioni angolari gaussiane in x e y e ruota
     * la direzione della particella di conseguenza.
     */
    void ApplicaScattering(Particle& particle, double thickness, const string& material)
    {
        // Calcola l'angolo RMS di scattering
        double theta0 = CalcoloTheta0(particle.p, particle.beta, thickness, 
                                      particle.charge, material);
        
        // Genera deviazioni angolari gaussiane in x e y
        // Il fattore √2 viene dalla proiezione 2D
        double theta_x = fRandom.Gaus(0, theta0 * sqrt(2.0));
        double theta_y = fRandom.Gaus(0, theta0 * sqrt(2.0));
        
        // Applica le deviazioni angolari alla direzione della particella
        TVector3 newDir = particle.direction;
        
        // Ruota attorno all'asse x (deviazione in yz)
        newDir.RotateX(theta_x);
        
        // Ruota attorno all'asse y (deviazione in xz)
        newDir.RotateY(theta_y);
        
        // Rinormalizza il versore di direzione
        newDir = newDir.Unit();
        
        // Aggiorna la direzione e gli angoli della particella
        particle.direction = newDir;
        particle.theta = newDir.Theta();  // Angolo polare
        particle.phi = newDir.Phi();      // Angolo azimutale
    }
};

// ============================================================
// CLASSE TRACK INTERSECTOR
// ============================================================

/**
 * @class TrackIntersector
 * @brief Calcola le intersezioni delle tracce con i layer del detector
 * 
 * Data una particella con posizione e direzione iniziali, calcola
 * dove interseca un layer cilindrico di dato raggio.
 */
class TrackIntersector
{
private:
    config fConfig;  // Configurazione del detector

public:
    /**
     * @brief Costruttore della classe TrackIntersector
     * @param config Configurazione del detector
     */
    TrackIntersector(const config& cfg) : fConfig(cfg) {}

    /**
     * @brief Calcola l'intersezione di una traccia con un cilindro
     * @param pos Posizione iniziale della particella [cm]
     * @param dir Direzione della particella (versore)
     * @param raggio Raggio del cilindro [cm]
     * @param intersezione Punto di intersezione (output) [cm]
     * @param t Parametro di percorso (output) [cm]
     * @return true se c'è intersezione, false altrimenti
     * 
     * Risolve l'equazione parametrica della retta con l'equazione del cilindro:
     * (x₀ + c₁t)² + (y₀ + c₂t)² = R²
     * At² + Bt + C = 0
     */
    bool CalcoloIntersezione(const TVector3& pos, const TVector3& dir, 
                             double raggio, TVector3& intersezione, double& t) const
    {
        double x0 = pos.X();  // Coordinata x iniziale
        double y0 = pos.Y();  // Coordinata y iniziale
        double z0 = pos.Z();  // Coordinata z iniziale
        
        double c1 = dir.X();  // Componente x della direzione
        double c2 = dir.Y();  // Componente y della direzione
        double c3 = dir.Z();  // Componente z della direzione
        
        // Coefficienti dell'equazione quadratica
        double A = c1*c1 + c2*c2;                    // a = c₁² + c₂²
        double B = 2.0 * (x0*c1 + y0*c2);            // b = 2(x₀c₁ + y₀c₂)
        double C = x0*x0 + y0*y0 - raggio*raggio;    // c = x₀² + y₀² - R²
        
        // Discriminante
        double delta = B*B - 4.0*A*C;
        
        // Se delta < 0, la retta non interseca il cilindro
        if (delta < 0) return false;
        
        // Soluzioni dell'equazione quadratica
        double sqrt_delta = sqrt(delta);
        double t1 = (-B + sqrt_delta) / (2.0*A);
        double t2 = (-B - sqrt_delta) / (2.0*A);
        
        // Scegli la soluzione positiva più piccola (prima intersezione)
        t = 1e9;
        if (t1 > 0 && t1 < t) t = t1;
        if (t2 > 0 && t2 < t) t = t2;
        
        // Se nessuna soluzione è positiva, non c'è intersezione in avanti
        if (t > 1e8) return false;
        
        // Calcola il punto di intersezione
        intersezione = pos + t * dir;
        
        // Verifica che l'intersezione sia entro la lunghezza del detector
        if (fabs(intersezione.Z()) > fConfig.detectorLength/2.0) {
            return false;
        }
        
        return true;
    }
};

// ============================================================
// CLASSE METROPOLIS VERTICE RECONSTRUCTOR
// ============================================================

/**
 * @class MetropolisVertexReconstructor
 * @brief Ricostruisce il vertice usando l'algoritmo Metropolis-Hastings
 * 
 * Implementa un algoritmo Markov Chain Monte Carlo (MCMC) per trovare
 * la posizione del vertice che massimizza la likelihood dei tracklets.
 */
class MetropolisVertexReconstructor
{
private:
    config fConfig;           // Configurazione del detector
    TRandom3 fRandom;         // Generatore di numeri casuali
    double fStepSize;         // Dimensione del passo per Metropolis [cm]
    int fIterations;          // Numero di iterazioni
    double fBestZ;            // Migliore stima di z trovata [cm]
    double fBestLLH;          // Log-Likelihood della migliore stima

public:
    /**
     * @brief Costruttore della classe MetropolisVertexReconstructor
     * @param config Configurazione del detector
     * @param stepSize Dimensione del passo per Metropolis [cm]
     * @param iterations Numero di iterazioni
     */
    MetropolisVertexReconstructor(const config& cfg, 
                                  double stepSize = 0.01, 
                                  int iterations = 5000) 
        : fConfig(cfg), fRandom(0), fStepSize(stepSize), 
          fIterations(iterations), fBestZ(0), fBestLLH(-1e9) {}

    /**
     * @brief Stima la posizione del vertice usando Metropolis-Hastings
     * @param tracklets Vettore di tracklets validi
     * @param initialZ Stima iniziale di z [cm]
     * @return Stima ottimale di z del vertice [cm]
     * 
     * L'algoritmo Metropolis-Hastings è una catena di Markov che
     * cammina nello spazio dei parametri (z) accettando passi che
     * migliorano la likelihood o occasionalmente passi che la peggiorano.
     * Questo permette di evitare massimi locali.
     */
    double FitVertex(const vector<Tracklet>& tracklets, double initialZ)
    {
        // Inizializza con la stima iniziale
        double currentZ = initialZ;
        double currentLLH = CalculateLLH(tracklets, currentZ);
        
        // Inizializza la migliore stima
        fBestZ = currentZ;
        fBestLLH = currentLLH;
        
        // Iterazioni dell'algoritmo Metropolis-Hastings
        for (int i = 0; i < fIterations; i++)
        {
            // Proponi un nuovo passo
            double proposedZ = currentZ + fRandom.Gaus(0.0, fStepSize);
            double proposedLLH = CalculateLLH(tracklets, proposedZ);
            
            // Calcola il rapporto di accettazione
            double deltaLLH = proposedLLH - currentLLH;
            double acceptanceRatio = min(1.0, exp(deltaLLH));
            
            // Accetta o rifiuta il passo
            if (fRandom.Uniform() < acceptanceRatio)
            {
                // Accetta il passo proposto
                currentZ = proposedZ;
                currentLLH = proposedLLH;
                
                // Aggiorna la migliore stima trovata
                if (currentLLH > fBestLLH)
                {
                    fBestZ = currentZ;
                    fBestLLH = currentLLH;
                }
            }
            // Se il passo è rifiutato, rimani nella posizione corrente
        }
        
        // Restituisce la migliore stima trovata
        return fBestZ;
    }

private:
    /**
     * @brief Calcola la Log-Likelihood per una data posizione del vertice
     * @param tracklets Vettore di tracklets validi
     * @param zVert Posizione z del vertice da testare [cm]
     * @return Valore della Log-Likelihood
     * 
     * La likelihood è basata sull'assunzione che i residui (differenze
     * tra l'intersezione del tracklet e la posizione del vertice)
     * siano distribuiti normalmente con σ dato dall'errore del tracklet.
     */
    double CalculateLLH(const vector<Tracklet>& tracklets, double zVert)
    {
        double llh = 0.0;
        
        // Somma i contributi di tutti i tracklets validi
        for (const auto& tracklet : tracklets)
        {
            if (!tracklet.valid) continue;
            
            // Residuo: differenza tra intersezione e vertice
            double residual = tracklet.z_intersection - zVert;
            
            // Errore del tracklet (può essere diverso per ogni tracklet)
            double sigma = tracklet.sigma_z;
            
            // Contributo alla Log-Likelihood: -½(residual/σ)²
            // (trascurando il termine costante che non influisce sul massimo)
            llh += -0.5 * residual * residual / (sigma * sigma);
        }
        
        return llh;
    }
};

// ============================================================
// CLASSE VERTEX RECONSTRUCTOR (PRINCIPALE)
// ============================================================

/**
 * @class VertexReconstructor
 * @brief Ricostruisce il vertice di collisione dai tracklets
 * 
 * Questa è la classe principale per la ricostruzione del vertice.
 * Implementa:
 * 1. Formazione dei tracklets da hit nei layer 1 e 2
 * 2. Stima iniziale del vertice (mediana)
 * 3. Ricostruzione iterativa robusta
 * 4. Gestione degli errori e qualità del fit
 */
class VertexReconstructor
{
private:
    config fConfig;      // Configurazione del detector
    TRandom3 fRandom;    // Generatore di numeri casuali
    
public:
    /**
     * @brief Costruttore della classe VertexReconstructor
     * @param config Configurazione del detector
     */
    VertexReconstructor(const config& cfg) : fConfig(cfg), fRandom(0) {}
    
    /**
     * @brief Forma i tracklets abbinando hit tra layer 1 e 2
     * @param hitsLayer1 Vettore di hit nel layer 1
     * @param hitsLayer2 Vettore di hit nel layer 2
     * @return Vettore di tracklets validi
     * 
     * I tracklets sono formati abbinando hit nei due layer che:
     * 1. Hanno Δφ < deltaPhiCut
     * 2. Hanno Δz ragionevole (entro 10σ)
     * 3. Producono un'intersezione plausibile
     * 
     * Algoritmo ottimizzato: O(n log n) invece di O(n²) tramite:
     * - Ordinamento degli hit per φ
     * - Finestra scorrevole per il matching
     */
    vector<Tracklet> FormTracklets(const vector<Hit>& hitsLayer1, 
                                   const vector<Hit>& hitsLayer2) const
    {
        vector<Tracklet> tracklets;
        
        // ----- PREPARAZIONE DEI DATI -----
        // Creiamo copie degli hit con i loro indici, ordinati per φ
        vector<pair<double, int>> sortedL1, sortedL2;
        
        // Layer 1: coppia (φ, indice_originale)
        for (size_t i = 0; i < hitsLayer1.size(); i++)
        {
            sortedL1.emplace_back(hitsLayer1[i].phi_smeared, i);
        }
        
        // Layer 2: coppia (φ, indice_originale)
        for (size_t i = 0; i < hitsLayer2.size(); i++)
        {
            sortedL2.emplace_back(hitsLayer2[i].phi_smeared, i);
        }
        
        // Ordina per φ (crescente)
        sort(sortedL1.begin(), sortedL1.end());
        sort(sortedL2.begin(), sortedL2.end());
        
        // ----- MATCHING CON FINESTRA SCORREVOLE -----
        // Inizializza l'indice di partenza per la finestra in L2
        size_t j_start = 0;
        
        // Itera su tutti gli hit del layer 1
        for (const auto& [phi1, idx1] : sortedL1)
        {
            // Avanza j_start fino a che φ_L2 >= φ_L1 - Δφ - margine
            while (j_start < sortedL2.size() && 
                   sortedL2[j_start].first < phi1 - fConfig.deltaPhiCut - 0.05)
            {
                j_start++;
            }
            
            // Itera sulla finestra in L2 (da φ1-Δφ a φ1+Δφ)
            for (size_t j = j_start; j < sortedL2.size(); j++)
            {
                const auto& [phi2, idx2] = sortedL2[j];
                
                // Esci se abbiamo superato φ1+Δφ+margine
                if (phi2 > phi1 + fConfig.deltaPhiCut + 0.05) break;
                
                // Calcola Δφ, gestendo la periodicità 2π
                double deltaPhi = fabs(phi1 - phi2);
                if (deltaPhi > TMath::Pi())
                {
                    deltaPhi = 2.0 * TMath::Pi() - deltaPhi;
                }
                
                // ----- TAGLIO IN Δφ -----
                if (deltaPhi < fConfig.deltaPhiCut)
                {
                    const Hit& hit1 = hitsLayer1[idx1];
                    const Hit& hit2 = hitsLayer2[idx2];
                    
                    // ----- TAGLIO IN Δz (opzionale) -----
                    // Taglio in differenza di z per ridurre falsi match
                    double deltaZ = fabs(hit1.z_smeared - hit2.z_smeared);
                    double deltaZCut = 10.0 * fConfig.smearZ;  // 10σ
                    if (deltaZ > deltaZCut) continue;
                    
                    // ----- CREA E VALIDA IL TRACKLET -----
                    Tracklet tracklet;
                    tracklet.hit1_idx = idx1;
                    tracklet.hit2_idx = idx2;
                    
                    if (CalculateTrackletIntersection(tracklet, hit1, hit2))
                    {
                        // Stima l'errore sul tracklet
                        tracklet.sigma_z = EstimateTrackletError(tracklet, hit1, hit2);
                        
                        // Calcola il χ² del fit
                        tracklet.chi2 = CalculateChi2(tracklet, hit1, hit2);
                        
                        // ----- TAGLIO IN χ² PER QUALITÀ -----
                        if (tracklet.chi2 < 10.0)  // Valore empirico
                        {
                            tracklets.push_back(tracklet);
                        }
                    }
                }
            }
        }
        
        return tracklets;
    }
    
    /**
     * @brief Ricostruzione semplice del vertice (mediana pesata)
     * @param tracklets Vettore di tracklets validi
     * @return Stima di z del vertice [cm]
     * 
     * Calcola la mediana delle intersezioni dei tracklets,
     * pesata con l'inverso della varianza (1/σ²).
     * La mediana è robusta contro outliers.
     */
    double ReconstructVertexSimple(const vector<Tracklet>& tracklets)
    {
        if (tracklets.empty()) return 0.0;
        
        // Vettore di coppie (z_intersection, peso)
        // Peso = 1/σ² (inverso della varianza)
        vector<pair<double, double>> weightedIntersections;
        
        // Raccogli tutte le intersezioni valide con i loro pesi
        for (const auto& tracklet : tracklets)
        {
            if (!tracklet.valid) continue;
            
            double weight = 1.0 / (tracklet.sigma_z * tracklet.sigma_z);
            weightedIntersections.emplace_back(tracklet.z_intersection, weight);
        }
        
        if (weightedIntersections.empty()) return 0.0;
        
        // Ordina per z (crescente)
        sort(weightedIntersections.begin(), weightedIntersections.end(),
             [](const pair<double,double>& a, const pair<double,double>& b) {
                 return a.first < b.first;
             });
        
        // Calcola il peso totale
        double totalWeight = 0.0;
        for (const auto& [z, w] : weightedIntersections)
        {
            totalWeight += w;
        }
        
        // Trova la mediana pesata
        double halfWeight = totalWeight / 2.0;
        double accumulated = 0.0;
        
        for (const auto& [z, w] : weightedIntersections)
        {
            accumulated += w;
            if (accumulated >= halfWeight)
            {
                return z;  // Mediana pesata
            }
        }
        
        // Fallback: mediana semplice (non pesata)
        size_t n = weightedIntersections.size();
        return weightedIntersections[n/2].first;
    }
    
    /**
     * @brief Ricostruzione iterativa robusta del vertice
     * @param tracklets Vettore di tracklets validi
     * @param maxIterations Numero massimo di iterazioni (default: 10)
     * @param nSigmaCut Taglio in σ per escludere outliers (default: 3σ)
     * @return Stima robusta di z del vertice [cm]
     * 
     * Algoritmo iterativo che:
     * 1. Stima iniziale con mediana
     * 2. Calcola residui normalizzati
     * 3. Pesa i tracklets con funzione robusta (Tukey-like)
     * 4. Ristima il vertice
     * 5. Ripete fino a convergenza
     * 
     * La funzione di peso robusta riduce l'influenza dei punti
     * lontani dal vertice (outliers).
     */
    double ReconstructVertexIterative(const vector<Tracklet>& tracklets, 
                                      int maxIterations = 10, 
                                      double nSigmaCut = 3.0)
    {
        // Almeno 3 tracklets per una stima robusta
        if (tracklets.size() < 3) return 0.0;
        
        // Stima iniziale (mediana semplice)
        double currentZ = ReconstructVertexSimple(tracklets);
        
        // Ciclo iterativo
        for (int iter = 0; iter < maxIterations; iter++)
        {
            double sumZ = 0.0;    // Somma pesata delle z
            double sumW = 0.0;    // Somma dei pesi
            int nUsed = 0;        // Numero di tracklets usati
            
            // Itera su tutti i tracklets
            for (const auto& tracklet : tracklets)
            {
                if (!tracklet.valid) continue;
                
                // Residuo normalizzato
                double residual = tracklet.z_intersection - currentZ;
                double normalizedResidual = fabs(residual) / tracklet.sigma_z;
                
                // ----- TAGLIO ITERATIVO SUGLI OUTLIERS -----
                if (normalizedResidual < nSigmaCut)
                {
                    // Peso base: inverso della varianza
                    double baseWeight = 1.0 / (tracklet.sigma_z * tracklet.sigma_z);
                    
                    // ----- FUNZIONE DI PESO ROBUSTA (TUKEY-LIKE) -----
                    // Il peso diminuisce quadraticamente con il residuo normalizzato
                    // Peso totale = baseWeight * [1 - (residuo/nSigmaCut)²]²
                    // Questo riduce l'influenza dei punti vicini al taglio
                    double robustFactor = 1.0 - (normalizedResidual/nSigmaCut) * 
                                                (normalizedResidual/nSigmaCut);
                    double robustWeight = baseWeight * robustFactor * robustFactor;
                    
                    // Aggiorna le somme
                    sumZ += tracklet.z_intersection * robustWeight;
                    sumW += robustWeight;
                    nUsed++;
                }
            }
            
            // Verifica di avere abbastanza punti
            if (nUsed < 3 || sumW < 1e-9) break;
            
            // Nuova stima (media pesata)
            double newZ = sumZ / sumW;
            
            // ----- CRITERIO DI CONVERGENZA -----
            // Convergenza se la variazione è < 10 μm
            if (fabs(newZ - currentZ) < 0.001)  // 10 μm
            {
                return newZ;
            }
            
            // Prepara la prossima iterazione
            currentZ = newZ;
        }
        
        // Restituisce l'ultima stima (non converguta entro maxIterations)
        return currentZ;
    }
    
private:
    /**
     * @brief Calcola l'intersezione di un tracklet con r=0
     * @param tracklet Tracklet da calcolare (input/output)
     * @param hit1 Primo hit (layer 1)
     * @param hit2 Secondo hit (layer 2)
     * @return true se il tracklet è valido, false altrimenti
     * 
     * Esegue un fit lineare in 2D (r→z) e estrapola a r=0.
     * I parametri del fit sono:
     * - slope = (z2 - z1) / (r2 - r1)
     * - z_intersection = z1 - slope * r1
     */
    bool CalculateTrackletIntersection(Tracklet& tracklet, 
                                       const Hit& hit1, 
                                       const Hit& hit2) const
    {
        double r1 = hit1.r;            // Raggio del primo hit
        double z1 = hit1.z_smeared;    // z del primo hit (con smearing)
        double r2 = hit2.r;            // Raggio del secondo hit
        double z2 = hit2.z_smeared;    // z del secondo hit (con smearing)
        
        // Evita divisione per zero (hit sullo stesso raggio)
        if (fabs(r2 - r1) < 1e-6)
        {
            tracklet.valid = false;
            return false;
        }
        
        // Calcola la pendenza (dz/dr)
        tracklet.slope = (z2 - z1) / (r2 - r1);
        
        // Estrapola a r=0: z = z1 - slope * r1
        tracklet.z_intersection = z1 - tracklet.slope * r1;
        
        // ----- CRITERI DI VALIDITÀ -----
        // 1. Intersezione entro ±30 cm (regione fisica plausibile)
        // 2. Pendenza ragionevole (<5, corrisponde a θ < ~78°)
        tracklet.valid = (fabs(tracklet.z_intersection) < 30.0) && 
                         (fabs(tracklet.slope) < 5.0);
        
        return tracklet.valid;
    }
    
    /**
     * @brief Stima l'errore sull'intersezione di un tracklet
     * @param tracklet Tracklet di cui stimare l'errore
     * @param hit1 Primo hit (layer 1)
     * @param hit2 Secondo hit (layer 2)
     * @return Errore stimato σ_z [cm]
     * 
     * L'errore totale ha tre componenti:
     * 1. Risoluzione intrinseca del detector (σ_z_hit)
     * 2. Propagazione geometrica (fattore r/Δr)
     * 3. Scattering multiplo (dipendente dall'angolo)
     */
    double EstimateTrackletError(const Tracklet& tracklet, 
                                 const Hit& hit1, 
                                 const Hit& hit2) const
    {
        // ----- 1. ERRORE BASE DAL DETECTOR -----
        // Risoluzione intrinseca in z (120 μm)
        double sigma_z_hit = fConfig.smearZ;
        
        // ----- 2. PROPAGAZIONE GEOMETRICA -----
        // Fattore geometrico: errore cresce quando r1/Δr è grande
        // (estrapolazione lunga → maggiore incertezza)
        double r1 = hit1.r;
        double r2 = hit2.r;
        double deltaR = r2 - r1;
        
        // Fattore di propagazione: √[1 + (r1/Δr)²]
        double geometric_factor = sqrt(1.0 + (r1/deltaR)*(r1/deltaR));
        
        // ----- 3. CONTRIBUTO DELLO SCATTERING MULTIPLO -----
        // Calcola l'angolo del tracklet
        double theta = 0.0;
        if (fabs(hit2.z_smeared - hit1.z_smeared) > 1e-6)
        {
            theta = atan2(deltaR, fabs(hit2.z_smeared - hit1.z_smeared));
        }
        
        // Errore da scattering: proporzionale a r1/sin(θ)
        // Approssimazione: 0.1 mrad per X0
        double scattering_error = 0.0;
        if (theta > 0.01)  // Evita divisione per numeri piccoli
        {
            scattering_error = 0.0001 * r1 / sin(theta);  // 0.1 mrad in rad
        }
        
        // ----- ERRORE TOTALE -----
        // Combinazione in quadratura dei contributi
        double sigma_total = sqrt(pow(sigma_z_hit * geometric_factor, 2) + 
                                  pow(scattering_error, 2));
        
        // Limite minimo di risoluzione: 5 μm
        // (al di sotto di questo, i miglioramenti non sono significativi)
        return max(sigma_total, 0.0005);  // Minimo 5 μm
    }
    
    /**
     * @brief Calcola il χ² di un tracklet
     * @param tracklet Tracklet da valutare
     * @param hit1 Primo hit (layer 1)
     * @param hit2 Secondo hit (layer 2)
     * @return Valore del χ²
     * 
     * Il χ² misura la compatibilità dei due hit con una retta.
     * Per due punti e due parametri (z0, slope), i gradi di libertà sono 0.
     * Questo calcolo è principalmente per monitoraggio e selezione.
     */
    double CalculateChi2(const Tracklet& tracklet, 
                         const Hit& hit1, 
                         const Hit& hit2) const
    {
        if (!tracklet.valid) return 1e6;
        
        double r1 = hit1.r;
        double z1 = hit1.z_smeared;
        double r2 = hit2.r;
        double z2 = hit2.z_smeared;
        
        // ----- PREDIZIONI DAL FIT -----
        // Il fit lineare passa esattamente per i due punti,
        // quindi le predizioni coincidono con le misure.
        // In un'implementazione reale con errori diversi sui due punti,
        // il fit potrebbe non passare esattamente.
        double z1_pred = tracklet.z_intersection + tracklet.slope * r1;
        double z2_pred = tracklet.z_intersection + tracklet.slope * r2;
        
        // ----- ERRORI SUGLI HIT -----
        // Assumiamo lo stesso errore per tutti gli hit
        double sigma_z = fConfig.smearZ;
        
        // ----- RESIDUI NORMALIZZATI -----
        double resid1 = (z1 - z1_pred) / sigma_z;
        double resid2 = (z2 - z2_pred) / sigma_z;
        
        // ----- χ² TOTALE -----
        // Somma dei quadrati dei residui normalizzati
        return resid1*resid1 + resid2*resid2;
    }
    
    /**
     * @brief Calcola la mediana di un vettore
     * @tparam T Tipo degli elementi (deve supportare ordinamento)
     * @param v Vettore di elementi
     * @return Mediana degli elementi
     * 
     * La mediana è robusta contro outliers rispetto alla media.
     * Per un numero pari di elementi, restituisce la media dei due centrali.
     */
    template<typename T>
    T Median(vector<T>& v) const
    {
        if (v.empty()) return T(0);
        
        // Ordina il vettore
        sort(v.begin(), v.end());
        
        size_t n = v.size();
        if (n % 2 == 0)  // Numero pari di elementi
        {
            // Media dei due elementi centrali
            return (v[n/2 - 1] + v[n/2]) / 2.0;
        }
        else  // Numero dispari di elementi
        {
            // Elemento centrale
            return v[n/2];
        }
    }
};
