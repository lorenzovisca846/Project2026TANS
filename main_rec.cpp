#include "classes/Config.h"
#include <iostream>

int main(int argc, char** argv)
{

    //================================= Config parameters =================================
    string cFile = "inputConfig.txt";
    if (argc > 1) cFile = argv[1];
    string configFile = "../config/" + cFile;

    TEnv *configEnv = new TEnv(configFile.c_str());

    SimRandom* simrand = new SimRandom(configEnv->GetValue("Seed",0));
    delete gRandom;
    gRandom = simrand;

    Config config(simrand, configEnv);
    config.Print();

    delete configEnv;

    return 0;
}