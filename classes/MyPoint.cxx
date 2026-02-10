#include "MyPoint.h"
#include "TObject.h"

ClassImp(MyPoint)

MyPoint::MyPoint(double R, double Phi, double Z):
    TObject(),
    fR(R),fPhi(Phi),fZ(Z),
    fTrackID(-1)
    {}

MyPoint::MyPoint(double R, double Phi, double Z, int trackID):
    TObject(),
    fR(R),fPhi(Phi),fZ(Z),
    fTrackID(trackID)
    {}