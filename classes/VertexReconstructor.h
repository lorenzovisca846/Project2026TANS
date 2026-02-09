#ifndef VERTEXRECONSTRUCTOR_H
#define VERTEXRECONSTRUCTOR_H
#include "Config.h"
#include "Hit.h"
#include "Tracklet.h"
#include <TRandom3.h>
#include <TMath.h>

using namespace std;

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
    Config fConfig;      // Configurazione del detector
    TRandom3 fRandom;    // Generatore di numeri casuali
    
public:
    /**
     * @brief Costruttore della classe VertexReconstructor
     * @param config Configurazione del detector
     */
    VertexReconstructor(const Config& cfg) : fConfig(cfg), fRandom(0) {}
    
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


#endif // VERTEXRECONSTRUCTOR_H