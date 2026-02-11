#include "Cylinder.h"
#include "TObject.h"

ClassImp(Cylinder)

Cylinder::Cylinder(double R, double L, double W, const string& Material, int layer):
    TObject(),
    fR(R),fL(L),fW(W),
    fX0(1000.0),
    fLayer(layer)
{
    if(Material=="Be")
        fX0 = 35.28;
    else if(Material=="Si")
        fX0 = 9.36;
}