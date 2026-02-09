#ifndef HIT_H
#define HIT_H

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

#endif // HIT_H