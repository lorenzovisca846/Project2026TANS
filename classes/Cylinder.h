#ifndef CYLINDER_H
#define CYLINDER_H

#include "TObject.h"
#include "string.h"
using namespace std;

/*
NOTA: soprattutto nei layer dei rivelatori la traiettoria non è radiale, nel calcolo del MS
          approssimiamo comunque la lunghezza del tragitto come lo spessore del layer?
*/

class Cylinder : public TObject
{
public:
    Cylinder();
    Cylinder(double R, double L, double W, string Material, bool msEnabled);

    //Se decidiamo di ereditare da Cylinder con beampipe e layer va messo virtual
    ~Cylinder();

    double GetR() const { return fR;}
    double GetL() const { return fL;}
    double GetW() const { return fW;}
    double GetMS() const { return fMsCoefficient;}

private:
    double fR;                  // radius (cm)
    double fL;                  // length (cm)
    double fW;                  // width (cm) 
    double fMsCoefficient;      // multiple scattering coefficient (GeV/c)
    bool fMsEnabled;

ClassDef(Cylinder,1)
};

#endif // CYLINDER_H