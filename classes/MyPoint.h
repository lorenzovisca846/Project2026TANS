#ifndef MYPOINT_H
#define MYPOINT_H

#include "TObject.h"
#include <cmath>

class MyPoint : public TObject
{
    public:
        MyPoint():TObject(),fZ(0.0),fR(0.0),fPhi(0.0),fTrackID(-1){}
        MyPoint(double R, double Phi, double Z);
        MyPoint(double R, double Phi, double Z, int trackID);

        double GetX() const {return fR*cos(fPhi);}
        double GetY() const {return fR*sin(fPhi);}
        double GetZ() const {return fZ;}

        double GetR() const {return fR;}
        double GetPhi() const {return fPhi;}

        double GetZSmeared() const {return fZ_smeared;}
        double GetPhiSmeared() const {return fPhi_smeared;}

        double GetXSmeared() const {return fR*cos(fPhi_smeared);}
        double GetYSmeared() const {return fR*sin(fPhi_smeared);}

        int GetTrackID() const {return fTrackID;}

        void SetSmeared(double zs, double phis) {fZ_smeared = zs; fPhi_smeared = phis;}

    private:
        double fZ;              // true z coordinate
        double fR;              // radial coordinate
        double fPhi;            // true azimut angle

        double fZ_smeared;      // smeared z coordinate
        double fPhi_smeared;    // smeared azimut angle

        int fTrackID;           // track ID

    ClassDef(MyPoint,1)
};



#endif //MYPOINT_H