#include "Cylinder.h"
#include "TObject.h"

ClassImp(Cylinder)

Cylinder::Cylinder():TObject(),
    fR(0.0),
    fL(0.0),
    fW(0.0),
    fMsCoefficient(0.0),
    fMsEnabled(false)
{
    // Default constructor
}

Cylinder::Cylinder(double R, double L, double W, string Material, bool msEnabled):TObject(),
    fR(R),
    fL(L),
    fW(W),
    fMsEnabled(msEnabled)
{
    if(fMsEnabled)
    {
        double X0 = 0.0;
        if(Material=="Be")
            X0 = 35.28;
        else if(Material=="Si")
            X0 = 9.36;

        fMsCoefficient = 0.0136 * sqrt(W/X0) * (1 + 0.038 * log(W/X0));
    }
    else
        fMsCoefficient = 0.0;
}

Cylinder::~Cylinder()
{
    // Destructor
}