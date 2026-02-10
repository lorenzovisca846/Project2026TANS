#include "Tracklet.h"

double Tracklet::CalculateChi2(double smearZ) const
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