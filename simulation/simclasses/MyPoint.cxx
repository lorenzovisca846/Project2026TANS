#include "MyPoint.h"
#include "TObject.h"

ClassImp(Point)

MyPoint::MyPoint(double X, double Y, double Z):
    TObject(),
    fX(X),fY(Y),fZ(Z)
{
    fR = sqrt(X*X + Y*Y);
    fPhi = atan2(Y,X);
}