#include "Point.h"
#include "TObject.h"

ClassImp(Point)

Point::Point():TObject(),
    fZ(0.0),
    fR(0.0),
    fPhi(0.0)
{
    // Default constructor
}

Point::Point(const double& X,const double& Y,const double& Z):TObject(),
    fZ(Z)
{
    fR = sqrt(X*X + Y*Y);
    fPhi = atan2(Y,X);
}