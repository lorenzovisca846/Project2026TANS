#ifndef METROPOLISVERTEXRECONSTRUCTOR_H
#define METROPOLISVERTEXRECONSTRUCTOR_H

#include "Config.h"
#include "Hit.h"
#include "Tracklet.h"
#include <vector>
#include <TRandom3.h>

using namespace std;

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

#endif // METROPOLISVERTEXRECONSTRUCTOR_H