
#ifndef TRACKLET_H
#define TRACKLET_H

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


#endif // TRACKLET_H