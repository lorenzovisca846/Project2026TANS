#ifndef CYLINDER_H
#define CYLINDER_H

#include "TObject.h"
#include "string.h"
using namespace std;

class Cylinder : public TObject
{
public:
    Cylinder();
    Cylinder(const double& R, const double& L, const double& W, string Material);

    double GetR() const { return fR;}
    double GetL() const { return fL;}
    double GetW() const { return fW;}
    double GetX0() const { return fX0;}

private:
    double fR;                  // radius (cm)
    double fL;                  // length (cm)
    double fW;                  // width (cm) 
    double fX0;                 // radiation length (cm)

ClassDef(Cylinder,1)
};

#endif // CYLINDER_H