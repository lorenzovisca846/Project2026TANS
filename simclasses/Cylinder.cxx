#include "Cylinder.h"
#include "TObject.h"

ClassImp(Cylinder)

Cylinder::Cylinder():TObject(),
    fR(0.0),
    fL(0.0),
    fW(0.0),
    fX0(1000.0)
{
    // Default constructor
}

Cylinder::Cylinder(const double& R, const double& L, const double& W, string Material):TObject(),
    fR(R),
    fL(L),
    fW(W),
    fX0(1000.0)
{
    if(Material=="Be")
        fX0 = 35.28;
    else if(Material=="Si")
        fX0 = 9.36;

    //MsCoefficient = 0.0136 * sqrt(W/X0) * (1 + 0.038 * log(W/X0));
}