#include "Point.h"
#include "TObject.h"

ClassImp(Point)

Point::Point():TObject(),
    fX(0.0),
    fY(0.0),
    fZ(0.0),
    fR(0.0),
    fPhi(0.0)
{
    // Default constructor
}

Point::Point(const double& X,const double& Y,const double& Z):TObject(),
    fX(X),
    fY(Y),
    fZ(Z)
{
    fR = sqrt(X*X + Y*Y);
    fPhi = atan2(Y,X);
}