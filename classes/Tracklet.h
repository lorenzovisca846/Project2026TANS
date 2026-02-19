
#ifndef TRACKLET_H
#define TRACKLET_H

#include "MyPoint.h"

class Tracklet
{
public:
    Tracklet() : hit1_idx(-1), hit2_idx(-1), z_intersection(0), 
                 slope(0){}
    Tracklet(int idx1, int idx2) : hit1_idx(idx1), hit2_idx(idx2), z_intersection(0), 
                                  slope(0){}
    
    void CalculateTrackletIntersection( const MyPoint& hit1, const MyPoint& hit2);

    int hit1_idx;               // Indice dell'hit nel vettore hitsLayer1
    int hit2_idx;               // Indice dell'hit nel vettore hitsLayer2
    double z_intersection;      // Coordinata z dell'intersezione a r=0 [cm]
    double slope;               // Pendenza dz/dr della tracklet
};


#endif // TRACKLET_H