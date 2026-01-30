#ifndef MYPOINT_H
#define MYPOINT_H

#include "TObject.h"

class MyPoint : public TObject
{
    public:
    MyPoint();
    MyPoint(const double& X,const double& Y,const double& Z);

    double GetX() const {return fX;}
    double GetY() const {return fY;}
    double GetZ() const {return fZ;}

    double GetR() const {return fR;}
    double GetPhi() const {return fPhi;}

    private:
        double fX;
        double fY;
        double fZ;

        double fR;
        double fPhi;

    ClassDef(MyPoint,1)
};



#endif //MYPOINT_H