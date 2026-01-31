#ifndef MYPOINT_H
#define MYPOINT_H

#include "TObject.h"

class MyPoint : public TObject
{
    public:
        MyPoint():TObject(),fX(0.0),fY(0.0),fZ(0.0),fR(0.0),fPhi(0.0){}
        MyPoint(double X, double Y, double Z);

        double GetX() const {return fX;}
        double GetY() const {return fY;}
        double GetZ() const {return fZ;}

        double GetR() const {return fR;}
        double GetPhi() const {return fPhi;}

    private:
        double fX;              // x coordinate
        double fY;              // y coordinate
        double fZ;              // z coordinate

        double fR;              // radial coordinate
        double fPhi;            // azimut angle

    ClassDef(MyPoint,1)
};



#endif //MYPOINT_H