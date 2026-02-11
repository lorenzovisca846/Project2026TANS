#include "VertexReconstructor.h"


vector<Tracklet> VertexReconstructor::FormTracklets(const vector<MyPoint>& hitsLayer1, 
                                       const vector<MyPoint>& hitsLayer2) const
{
    vector<Tracklet> tracklets;
        
    // ----- PREPARAZIONE DEI DATI -----
    // Creiamo copie degli hit con i loro indici, ordinati per φ
    vector<pair<double, int>> sortedL1, sortedL2;
    
    // Layer 1: coppia (φ, indice_originale)
    for (size_t i = 0; i < hitsLayer1.size(); i++)
    {
        sortedL1.emplace_back(hitsLayer1[i].GetPhi(), i);
    }
    
    // Layer 2: coppia (φ, indice_originale)
    for (size_t i = 0; i < hitsLayer2.size(); i++)
    {
        sortedL2.emplace_back(hitsLayer2[i].GetPhi(), i);
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
                const MyPoint& hit1 = hitsLayer1[idx1];
                const MyPoint& hit2 = hitsLayer2[idx2];
                
                // ----- TAGLIO IN Δz (opzionale) -----
                // Taglio in differenza di z per ridurre falsi match
                double deltaZ = fabs(hit1.GetZ() - hit2.GetZ());
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




double VertexReconstructor::ReconstructVertexSimple(const vector<Tracklet>& tracklets)
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





double VertexReconstructor::ReconstructVertexIterative(const vector<Tracklet>& tracklets, 
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



bool VertexReconstructor::CalculateTrackletIntersection(Tracklet& tracklet, const MyPoint& hit1, const MyPoint& hit2) const
{
    double r1 = hit1.GetR();            // Raggio del primo hit
    double z1 = hit1.GetZ();    // z del primo hit (con smearing)
    double r2 = hit2.GetR();            // Raggio del secondo hit
    double z2 = hit2.GetZ();    // z del secondo hit (con smearing)
    
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


  double VertexReconstructor::EstimateTrackletError(const Tracklet& tracklet, 
                                 const MyPoint& hit1, 
                                 const MyPoint& hit2) const
    {
        // ----- 1. ERRORE BASE DAL DETECTOR -----
        // Risoluzione intrinseca in z (120 μm)
        double sigma_z_hit = fConfig.smearZ;
        
        // ----- 2. PROPAGAZIONE GEOMETRICA -----
        // Fattore geometrico: errore cresce quando r1/Δr è grande
        // (estrapolazione lunga → maggiore incertezza)
        double r1 = hit1.GetR();
        double r2 = hit2.GetR();
        double deltaR = r2 - r1;
        
        // Fattore di propagazione: √[1 + (r1/Δr)²]
        double geometric_factor = sqrt(1.0 + (r1/deltaR)*(r1/deltaR));
        
        // ----- 3. CONTRIBUTO DELLO SCATTERING MULTIPLO -----
        // Calcola l'angolo del tracklet
        double theta = 0.0;
        if (fabs(hit2.GetZ() - hit1.GetZ()) > 1e-6)
        {
            theta = atan2(deltaR, fabs(hit2.GetZ() - hit1.GetZ()));
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



    double VertexReconstructor::CalculateChi2(const Tracklet& tracklet, 
                         const MyPoint& hit1, 
                         const MyPoint& hit2) const
    {
        if (!tracklet.valid) return 1e6;
        
        double r1 = hit1.GetR();
        double z1 = hit1.GetZ();
        double r2 = hit2.GetR();
        double z2 = hit2.GetZ();
        
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