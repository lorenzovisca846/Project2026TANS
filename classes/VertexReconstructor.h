#ifndef VERTEXRECONSTRUCTOR_H
#define VERTEXRECONSTRUCTOR_H
#include <TRandom3.h>
#include <TMath.h>

#include "Config.h"
#include "Hit.h"
#include "Tracklet.h"
#include "SimRandom.h"

using namespace std;


class VertexReconstructor
{
    public:
        VertexReconstructor(const Config& cfg, SimRandom* srnd):fConfig(cfg),fRandom(srnd){}

        

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
                                       const vector<Hit>& hitsLayer2) const;



        /**
         * @brief Ricostruzione semplice del vertice (mediana pesata)
         * @param tracklets Vettore di tracklets validi
         * @return Stima di z del vertice [cm]
         * 
         * Calcola la mediana delle intersezioni dei tracklets,
         * pesata con l'inverso della varianza (1/σ²).
         * La mediana è robusta contro outliers.
         */
        double ReconstructVertexSimple(const vector<Tracklet>& tracklets);



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
                                         double nSigmaCut = 3.0);

    private:
        Config fConfig;      // Configurazione del detector
        SimRandom* fRandom;    // Generatore di numeri casuali

        bool CalculateTrackletIntersection(Tracklet& tracklet, const Hit& hit1, const Hit& hit2) const;

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
        double EstimateTrackletError(const Tracklet& tracklet, const Hit& hit1, const Hit& hit2) const;

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
        double CalculateChi2(const Tracklet& tracklet, const Hit& hit1, const Hit& hit2) const;

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