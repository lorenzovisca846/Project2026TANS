#ifndef POINT_H
#define POINT_H

#include "TObject.h"

class Point : public TObject
{
    public:
    Point();
    Point(const double& X,const double& Y,const double& Z);

    double GetZ() const {return fZ;}
    double GetR() const {return fR;}
    double GetPhi() const {return fPhi;}

    private:
        double fZ;
        double fR;
        double fPhi;

    ClassDef(Point,1)
};



#endif //POINT_H