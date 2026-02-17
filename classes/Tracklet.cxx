#include "Tracklet.h"


void Tracklet::CalculateTrackletIntersection( const MyPoint& hit1, const MyPoint& hit2)
{
    double r1 = hit1.GetR();    // Raggio del primo hit
    double z1 = hit1.GetZ();    // z del primo hit (con smearing)
    double r2 = hit2.GetR();    // Raggio del secondo hit
    double z2 = hit2.GetZ();    // z del secondo hit (con smearing)

    slope = (z2 - z1) / (r2 - r1);
    
    // Estrapola a r=0: z = z1 - slope * r1
    z_intersection = z1 - slope * r1;
}