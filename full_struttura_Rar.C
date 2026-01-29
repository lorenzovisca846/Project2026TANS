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
// configurazione
// ============================================================ 

struct config
{
    // geometria
    double beamPipeRadius = 3.0;           // cm // berillio    
    double beamPipeThickness = 0.08;       // cm

    double layer1Radius = 4.0;             // cm
    double layer1Thickness = 0.02;         // cm

    double layer2Radius = 7.0;             // cm
    double layer2Thickness = 0.02;         // cm

    double detectorLength = 27.0;          // cm

    // parametri vertice
    double vertexZSigma = 5.3;             // cm
    double vertexXYSigma = 0.01;           // cm

    // smearing
    double smearZ = 0.0120;                // cm
    double smearRPhi = 0.0030;             // cm
    
    // Materiali
    double X0_Be = 35.28;                  // cm
    double X0_Si = 9.37;                   // cm
    double mspConstant = 0.0136;           // MeV


    // parametri simulazione
    int nEvents = 1000;
    int multiplicityMin = 20;
    int multiplicityMax = 80;
    double noiseRateLayer = 5.0;           // %  

    // Metropolis
    double metropolisStepSize = 0.01;
    int metropolisNSteps = 10000;

    // Taglio in phi per formare tracklets
    double deltaPhiCut = 0.1;              // rad

    void Print()
    {
        cout << "=== CONFIGURAZIONE ===" << endl;
        cout << "Geometria: BeamPipe R=" << beamPipeRadius
             << "cm, Layer1 R=" << layer1Radius
             << "cm, Layer2 R=" << layer2Radius << "cm" << endl;
        cout << "Vertice: σ_z=" << vertexZSigma
             << "cm, σ_xy=" << vertexXYSigma << "cm" << endl;
        cout << "Delta Phi Cut: " << deltaPhiCut << " rad" << endl;
        cout << "Eventi: " << nEvents << endl;
    }
};

// ============================================================
// CLASSi dei dati
// ============================================================

class Particle
{
public: 
    int trackID;
    TVector3 position;
    TVector3 direction;
    double phi;
    double theta;
    double p;
    double pt;
    int charge;
    double beta;    

    // Costruttore di Particle
    Particle() : trackID(-1), theta(0), phi(0), p(0), pt(0), charge(0), beta(1.0) {}
};

class Hit
{
public:
    int layer;
    double x;
    double y;
    double z;
    double r;
    double phi;
    double z_smeared;
    double phi_smeared;
    bool isSignal;
    int trackID;   

    // Costruttore di Hit
    Hit() : layer(0), x(0), y(0), z(0), r(0), phi(0), z_smeared(0), phi_smeared(0), isSignal(true), trackID(-1) {}
};

class Tracklet
{
public:
    Hit hit1;
    Hit hit2;
    double z_intersection;
    double slope;
    double chi2;
    bool valid;

    // Costruttore di Tracklet
    Tracklet() : z_intersection(0), slope(0), chi2(0), valid(false) {}
    double CalculateChi2()
    {
        return 1.0; // Semplificato
    }
};

class Event
{
public:
    int eventID;
    int multiplicity;
    double trueVertexZ;
    double recoVertexZ;
    double recoVertexZMetropolis;
    bool vertexFound;
    int nTracklets;

    // Costruttore di Event
    Event() : eventID(0), multiplicity(0), trueVertexZ(0), recoVertexZ(0), recoVertexZMetropolis(0), vertexFound(false), nTracklets(0) {}

};

// ============================================================
// CLASSI DETECTOR
// ============================================================

class Detector
{
private:
    config fConfig;

public:
    // Costruttore di Detector, usa una lista di inizializzazione
    Detector(const config& config) : fConfig(config) {} 

    bool Accettanza(const TVector3& punto , int layer) const
    {
        double raggio = (layer ==1) ? fConfig.layer1Radius : fConfig.layer2Radius;
        double r = punto.Perp(); //.Perp() calcola la distanza radiale dal centro
        return (fabs(r-raggio) < fConfig.layerThickness && fabs(punto.Z()) < fConfig.detectorLength/2);
    }

    void Geometria (TVirtualPad* pad) const //TVirtualPad e' una classe di ROOT che rappresenta un'area di disegno
    {
        pad->cd(); 
        //cd=change directory, va a disegna in sto pad

        // crea un'ellisse centrata in (0,0) con raggio beamPipeRadius, perchè un'ellisse e non un cerchio? boh
        TEllipse* pipeInner = new TEllipse(0,0,fConfig.beamPipeRadius); 
        pipeInner->SetFillStyle(0); 
        pipeInner->SetLineColor(kBlue);
        pipeInner->SetLineWidth(2);
        pipeInner->Draw(); 
        
        TEllipse* pipeOuter = new TEllipse(0,0,fConfig.beamPipeRadius+fConfig.beamPipeThickness);
        pipeOuter->SetFillStyle(0);    
        pipeOuter->SetFillColor(kBlue);
        pipeOuter->SetLineWidth(2);
        pipeOuter->Draw("same"); // same per disegnare sullo stesso pad senza cancellare quello che c'era prima

        
        
        TEllipse* layer1Inner = new TEllipse(0,0,fConfig.layer1Radius);
        layer1Inner->SetFillStyle(0);    
        layer1Inner->SetLineColor(kRed);
        layer1Inner->SetLineWidth(2);
        layer1Inner->Draw("same"); // nartra vorta
        
        TEllipse* layer1Outer = new TEllipse(0,0,fConfig.layer1Radius+fConfig.layer1Thickness);
        layer1Outer->SetFillStyle(0);    
        layer1Outer->SetLineColor(kRed);
        layer1Outer->SetLineWidth(2);
        layer1Outer->Draw("same"); 
        
        

        TEllipse* layer2Inner = new TEllipse(0,0,fConfig.layer2Radius);
        layer2Inner->SetFillStyle(0);   
        layer2Inner->SetLineColor(kGreen);
        layer2Inner->SetLineWidth(2);
        layer2Inner->Draw("same"); 

        TEllipse* layer2Outer = new TEllipse(0,0,fConfig.layer2Radius+fConfig.layer2Thickness);
        layer2Outer->SetFillStyle(0);    
        layer2Outer->SetLineColor(kGreen);
        layer2Outer->SetLineWidth(2);
        layer2Outer->Draw("same"); 
        


        TLine* zAxis = new TLine(-fConfig.detectorLength/2,0,fConfig.detectorLength/2,0);
        zAxis->SetLineColor(kBlack);
        zAxis->SetLineStyle(2);
        zAxis->Draw("same");
    }
};

// ============================================================
// CLASSE SCATTERING
// ============================================================

class MultipleScattering {
private:
    Config fConfig;
    TRandom3 fRandom;  

public:
    // Costruttore di MultipleScattering, usa una lista di inizializzazione
    MultipleScattering(const Config& config) : fConfig(config), fRandom(0)

    double CalcoloTheta0(double p, double beta, double thickness, double charge , const string& material) const
    {
        double X0 = (material == "Be") ? fConfig.X0_Be : fConfig.X0_Si;
        
        double theta0 = (fConfig.mspConstant * fabs(charge) / (beta * p)) * sqrt(thickness / X0) * (1.0 + 0.038 * log(thickness / X0));
        return theta0;    
    }

    void ApplicaScattering(Particle& particle, double thickness, const string& material)
    {
        double theta0= CalcoloTheta0(particle.p, particle.beta, thickness, particle.charge, material);
        
        double theta_x = fRandom.Gaus(0, theta0*sqrt(2)); 
        double theta_y = fRandom.Gaus(0, theta0*sqrt(2));

        //applica la deviazione angolare alla direzione della particella
        TVector3 newDir=particle.direction;
        newDir.RotateX(theta_x);
        newDir.RotateY(theta_y);
        newDir.Normalize();
     
        particle.direction = newDir;
        particle.theta = newDir.Theta();
        particle.phi = newDir.Phi();
    }
    
};

// ============================================================
// E CLASSI PAA TRACCIA
// ============================================================

class TrackIntersector {
private:
    Config fConfig; 

public:
    // Costruttore di TrackIntersector, usa una lista di inizializzazione   
    TrackIntersector(const Config& config) : fConfig(config) {}


    bool CalcoloIntersezione(const TVector3& pos, const TVector3& dir, double raggio, TVector3& intersezione, double& t) const
    {
        double x0 = pos.X();
        double y0 = pos.Y(); 
        double z0 = pos.Z();

        double c1 = dir.X();
        double c2 = dir.Y();
        double c3 = dir.Z();

        double a = c1*c1 + c2*c2;
        double b = 2.0 * (x0*c1 + y0*c2);
        double c = x0*x0 + y0*y0 - raggio*raggio;

        double delta = b*b - 4.0*a*c;
        if (delta < 0) return false;        

        double t1 = (-b + sqrt(delta)) / (2.0*a);
        double t2 = (-b - sqrt(delta)) / (2.0*a);
        //ma guarda te se devo impostare una cazzo di quadratica 

        t = (t1>0 && (t1<t2 || t2<0)) ? t1 : t2;
        if (t <= 0) return false;

        intersezione = pos + t*dir;

        if (fabs(intersezione.Z()) > fConfig.detectorLength/2) 
        {
            return false;
        }
        return true;
    }
};

// ============================================================
// CLASSE METROPOLIS
// ============================================================

class MetropolisVertexReconstructor {
private:
    Config fConfig;
    TRandom3 fRandom;
    double fVertexZ;
    double fCurrentLLH; // log-likelihood corrente questo mea suggerito copailot
    double fStepSize;
    int fIterations;
    double fBestZ; 
    double fBestLLH;

public:
    MetropolisVertexReconstructor(double stepSize = 0.01, int iterations = 5000) : fRandom(0), fStepSize(stepSize), fIterations(iterations),  fBestZ(0), fBestLLH(-1e9) {}

    double FitVertex(const vector<Tracklet>& tracklets, double initialZ)
    {
        double currentZ = initialZ;
        double currentLLH = CalculateLLH(tracklets, currentZ);

        fBestZ = currentZ;
        fBestLLH = currentLLH;

        for (int i = 0; i < fIterations; i++)
        {
            double proposedZ = currentZ + fRandom.Gaus(0.0, fStepSize);
            double proposedLLH = CalculateLLH(tracklets, proposedZ);
            double deltaLLH = proposedLLH - currentLLH;
            double acceptanceRatio = min(1.0, exp(deltaLLH));

            if (fRandom.Uniform() < acceptanceRatio)
            {
                currentZ = proposedZ;
                currentLLH = proposedLLH;

                if (currentLLH > fBestLLH)
                {
                    fBestZ = currentZ;
                    fBestLLH = currentLLH;
                }
            }
        return fBestZ;

        // in pratica è un random walk nello spazio dei parametri Z,
        // accettando passi che migliorano la likelihood e occasionalmente passi che la peggiorano, 
        // permettendo di trovare il massimo globale della funzione di likelihood anche in presenza di massimi locali.
            


    }

    double CalculateLLH(const vector<Tracklet>& tracklets, double zVert)
    {
        double llh = 0.0;
        double sigma = 0.012; // smearing in z delle hit
        
        for (const auto& tracklet : tracklets)
        {
            if (!tracklet.valid) continue;
            double residual = tracklet.z_intersection - zVert;
            llh += -0.5 * residual * residual / (sigma * sigma);
        }
        return llh;
    }
};
}

// ============================================================
// CLASSE RICOSTRUTTORE PUTTANA LA VELOCITA'
// ============================================================

class VertexReconstructor {
private:
    Config fConfig;
    TRandom3 fRandom;
    
    // usa indici invece di copie
    struct Tracklet {
        int hit1_idx;           // Indice nel vettore hitsLayer1
        int hit2_idx;           // Indice nel vettore hitsLayer2
        double z_intersection;
        double slope;
        double chi2;
        double sigma_z;         // Errore stimato su z_intersection
        bool valid;
        
        Tracklet() : hit1_idx(-1), hit2_idx(-1), z_intersection(0), slope(0), chi2(0), sigma_z(0.01), valid(false) {}
    };
    
public:
    VertexReconstructor(const Config& config) : fConfig(config), fRandom(0) {}
    
    vector<Tracklet> FormTracklets(const vector<Hit>& hitsLayer1, const vector<Hit>& hitsLayer2) const {
        vector<Tracklet> tracklets;
        
        // copie ordinate per phi
        vector<pair<double, int>> sortedL1, sortedL2;
        
        // Layer 1 con indici
        for (size_t i = 0; i < hitsLayer1.size(); i++) {
            sortedL1.emplace_back(hitsLayer1[i].phi_smeared, i);
        }
        //emplace_back costruisce l'oggetto direttamente nel vettore senza fare copie inutili
        // mentre push_back richiede la copia dell'oggetto

        // Layer 2 con indici  
        for (size_t i = 0; i < hitsLayer2.size(); i++) {
            sortedL2.emplace_back(hitsLayer2[i].phi_smeared, i);
        }
        
        // Ordina per φ
        sort(sortedL1.begin(), sortedL1.end());
        sort(sortedL2.begin(), sortedL2.end());
        //sort ordnia in ordine crescente
        
        // finestra scorrevole per Δφ
        size_t j_start = 0;
        
        for (const auto& [phi1, idx1] : sortedL1) {
            // Trova inizio finestra in L2 (φ1 - Δφ)
            while (j_start < sortedL2.size() && sortedL2[j_start].first < phi1 - fConfig.deltaPhiCut - 0.05) //0.05 margine di sicurezza
            {
                j_start++;
            }
            
            // Scorri finestra in L2 (da φ1-Δφ a φ1+Δφ)
            for (size_t j = j_start; j < sortedL2.size(); j++) {
                const auto& [phi2, idx2] = sortedL2[j];
                
                // Esci se superi φ1+Δφ
                if (phi2 > phi1 + fConfig.deltaPhiCut + 0.05) break;
                
                double deltaPhi = fabs(phi1 - phi2);
                if (deltaPhi > TMath::Pi()) {
                    deltaPhi = 2*TMath::Pi() - deltaPhi;
                }
                
                // taglio in Δφ
                if (deltaPhi < fConfig.deltaPhiCut) {
                    const Hit& hit1 = hitsLayer1[idx1];
                    const Hit& hit2 = hitsLayer2[idx2];
                    
                    // taglio in Δz 
                    double deltaZ = fabs(hit1.z_smeared - hit2.z_smeared);
                    double deltaZCut = 10.0 * fConfig.smearZ;  // 10σ
                    if (deltaZ > deltaZCut) continue;
                    
                    Tracklet tracklet;
                    tracklet.hit1_idx = idx1;
                    tracklet.hit2_idx = idx2;
                    
                    if (CalculateTrackletIntersection(tracklet, hit1, hit2)) {
                        // calcolo di sigma_z in base all'angolo
                        tracklet.sigma_z = EstimateTrackletError(tracklet, hit1, hit2);
                        tracklet.chi2 = CalculateChi2(tracklet, hit1, hit2);
                        
                        // taglio in chi quadro per qualità del tracklet
                        if (tracklet.chi2 < 10.0) {  // valore scelto un po' a caso
                            tracklets.push_back(tracklet);
                        }
                    }
                }
            }
        }
        
        return tracklets;
    }
    
    double ReconstructVertexSimple(const vector<Tracklet>& tracklets) {
        if (tracklets.empty()) return 0.0;
        
        // Mediana pesata con errori
        vector<pair<double, double>> weightedIntersections; // (z, weight=1/σ²)
        
        for (const auto& tracklet : tracklets) {
            if (!tracklet.valid) continue;
            double weight = 1.0 / (tracklet.sigma_z * tracklet.sigma_z);
            weightedIntersections.emplace_back(tracklet.z_intersection, weight);
        }
        
        if (weightedIntersections.empty()) return 0.0;
        
        // Ordina per z
        sort(weightedIntersections.begin(), weightedIntersections.end(),
             [](const pair<double,double>& a, const pair<double,double>& b) {
                 return a.first < b.first;
             });
        
        // Trova mediana 
        double totalWeight = 0.0;
        for (const auto& [z, w] : weightedIntersections) {
            totalWeight += w;
        }
        
        double halfWeight = totalWeight / 2.0;
        double accumulated = 0.0;
        
        for (const auto& [z, w] : weightedIntersections) {
            accumulated += w;
            if (accumulated >= halfWeight) {
                return z;
            }
        }
        
        // Fallback: mediana semplice
        return weightedIntersections[weightedIntersections.size()/2].first;
    }
    
    // Ricostruzione iterativa con taglio sugli outlier
    double ReconstructVertexIterative(const vector<Tracklet>& tracklets, int maxIterations = 10, double nSigmaCut = 3.0) {
        if (tracklets.size() < 3) return 0.0;
        
        double currentZ = ReconstructVertexSimple(tracklets);
        
        for (int iter = 0; iter < maxIterations; iter++) {
            double sumZ = 0.0;
            double sumW = 0.0;
            int nUsed = 0;
            
            for (const auto& tracklet : tracklets) {
                if (!tracklet.valid) continue;
                
                double residual = tracklet.z_intersection - currentZ;
                double normalizedResidual = fabs(residual) / tracklet.sigma_z;
                
                // Taglio iterativo sugli outlier
                if (normalizedResidual < nSigmaCut) {
                    double weight = 1.0 / (tracklet.sigma_z * tracklet.sigma_z);
                    
                    // Peso robusto per ridurre l'influenza dei punti vicini al taglio
                    double heavyWeight = weight * 
                        (1.0 - (normalizedResidual/nSigmaCut)*(normalizedResidual/nSigmaCut));
                    
                    sumZ += tracklet.z_intersection * heavyWeight;
                    sumW += heavyWeight;
                    nUsed++;
                }
            }
            
            if (nUsed < 3 || sumW < 1e-9) break;
            
            double newZ = sumZ / sumW;
            
            // Convergenza
            if (fabs(newZ - currentZ) < 0.001) {  // 10 μm
                return newZ;
            }
            
            currentZ = newZ;
        }
        
        return currentZ;
    }
    
private:
    bool CalculateTrackletIntersection(Tracklet& tracklet, const Hit& hit1, const Hit& hit2) const {
        double r1 = hit1.r;
        double z1 = hit1.z_smeared;
        double r2 = hit2.r;
        double z2 = hit2.z_smeared;
        
        if (fabs(r2 - r1) < 1e-6) {
            tracklet.valid = false;
            return false;
        }
        
        tracklet.slope = (z2 - z1) / (r2 - r1);
        tracklet.z_intersection = z1 - tracklet.slope * r1;
        
        // Validità del tracklet
        tracklet.valid = (fabs(tracklet.z_intersection) < 30.0) &&  // entro 30 cm
                         (fabs(tracklet.slope) < 5.0);              // slope arbitrario
        
        return tracklet.valid;
    }
    
    double EstimateTrackletError(const Tracklet& tracklet, const Hit& hit1, const Hit& hit2) const {
        // Errore base da risoluzione del detector
        double sigma_z_hit = fConfig.smearZ;  // 120 μm
        
        // Propagazione errore dal fit lineare
        double r1 = hit1.r;
        double r2 = hit2.r;
        double deltaR = r2 - r1;
        
        // Fattore geometrico: errore cresce con angolo piccolo
        double geometric_factor = sqrt(1.0 + (r1/deltaR)*(r1/deltaR));
        
        // Errore angolare (da scattering multiplo)
        double theta = atan2(deltaR, fabs(hit2.z_smeared - hit1.z_smeared));
        double scattering_error = 0.0;
        if (theta > 0.01) {
            // Approssimazione scattering: 0.1 mrad per X0
            scattering_error = 0.0001 * r1 / sin(theta);  // 0.1 mrad in rad
        }
        
        // Errore totale
        double sigma_total = sigma_z_hit * geometric_factor + scattering_error;
        
        // Non scendere sotto il limite intrinseco
        return max(sigma_total, 0.0005);  // Minimo 5 μm
    }
    
    double CalculateChi2(const Tracklet& tracklet, const Hit& hit1, const Hit& hit2) const {
        // Fit lineare
        double r1 = hit1.r;
        double z1 = hit1.z_smeared;
        double r2 = hit2.r;
        double z2 = hit2.z_smeared;
        
        // Predizioni dal fit
        double z1_pred = tracklet.z_intersection + tracklet.slope * r1;
        double z2_pred = tracklet.z_intersection + tracklet.slope * r2;
        
        // Errori sugli hit
        double sigma_z = fConfig.smearZ;
        
        // Residui normalizzati
        double resid1 = (z1 - z1_pred) / sigma_z;
        double resid2 = (z2 - z2_pred) / sigma_z;
        
        // χ² per 2 punti - 2 parametri = 0 gradi di libertà
        return resid1*resid1 + resid2*resid2;
    }
    
    // Funzione per trovare la mediana di un vettore
    template<typename T>
    T Median(vector<T>& v) const {
        if (v.empty()) return T(0);
        sort(v.begin(), v.end());
        size_t n = v.size();
        if (n % 2 == 0) {
            return (v[n/2 - 1] + v[n/2]) / 2.0;
        } else {
            return v[n/2];
        }
    }
};