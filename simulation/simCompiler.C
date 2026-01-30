#include "TSystem.h"
#include "TString.h"


void simCompiler(TString myopt="fast")
{
    TString opt;
    if(myopt.Contains("force"))
        opt = "kfg";
    else
        opt = "kg";

    gSystem->CompileMacro("simclasses/Cylinder.cxx",opt.Data());
    gSystem->CompileMacro("simclasses/MyPoint.cxx",opt.Data());
    gSystem->CompileMacro("simclasses/SimRandom.cxx",opt.Data());
    gSystem->CompileMacro("simclasses/Particle.cxx",opt.Data());
    gSystem->CompileMacro("simulation.C",opt.Data());
}