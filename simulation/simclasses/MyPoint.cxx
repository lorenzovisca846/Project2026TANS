#include "MyPoint.h"
#include "TObject.h"

ClassImp(Point)

MyPoint::MyPoint(double X, double Y, double Z):
    TObject(),
    fX(X),fY(Y),fZ(Z),
    fTrackID(-1)
{
    fR = sqrt(X*X + Y*Y);
    fPhi = atan2(Y,X);
}

MyPoint::MyPoint(double X, double Y, double Z, int trackID):
    TObject(),
    fX(X),fY(Y),fZ(Z),
    fTrackID(trackID)
{
    fR = sqrt(X*X + Y*Y);
    fPhi = atan2(Y,X);
}