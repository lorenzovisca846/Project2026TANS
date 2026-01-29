#include "Cylinder.h"
#include "TObject.h"

ClassImp(Cylinder)

Cylinder::Cylinder():TObject(),
    fR(0.0),
    fL(0.0),
    fW(0.0),
    fMsCoefficient(0.0)
{
    // Default constructor
}

Cylinder::Cylinder(const double& R, const double& L, const double& W, string Material):TObject(),
    fR(R),
    fL(L),
    fW(W)
{
    double X0 = 100.0;
    if(Material=="Be")
        X0 = 35.28;
    else if(Material=="Si")
        X0 = 9.36;

    fMsCoefficient = 0.0136 * sqrt(W/X0) * (1 + 0.038 * log(W/X0));
}

Cylinder::~Cylinder()
{
    // Destructor
}