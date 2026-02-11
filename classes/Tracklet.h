
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
    bool valid;                 // true se il tracklet è valido

    // Costruttore di default
    Tracklet() : hit1_idx(-1), hit2_idx(-1), z_intersection(0), 
                 slope(0), valid(false) {}

    Tracklet(int idx1, int idx2) : hit1_idx(idx1), hit2_idx(idx2), z_intersection(0), 
                                  slope(0), valid(false) {}
};


#endif // TRACKLET_H