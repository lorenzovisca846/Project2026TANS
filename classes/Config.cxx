#include "Config.h"
ClassImp(Config)

Config::Config(SimRandom* srnd, TEnv* configEnv):
    TObject(),
    fSimrand(srnd)
    {

    beamPipeRadius      = configEnv->GetValue("BeamPipeRadius", 3.0);
    beamPipeThickness   = configEnv->GetValue("BeamPipeWidth", 0.08);
    layer1Radius        = configEnv->GetValue("Layer1Radius", 4.0);
    layer1Thickness     = configEnv->GetValue("LayerWidth", 0.02);
    layer2Radius        = configEnv->GetValue("Layer2Radius", 7.0);
    layer2Thickness     = configEnv->GetValue("LayerWidth", 0.02);
    detectorLength      = configEnv->GetValue("Length", 27.0);
    vertexZSigma        = configEnv->GetValue("SigmaZ", 5.3);
    vertexXYSigma       = configEnv->GetValue("SigmaXY", 0.01);
    nEvents             = configEnv->GetValue("Events", 10000);
    multiplicityMin     = configEnv->GetValue("Minimum", 1);
    multiplicityMax     = configEnv->GetValue("Maximum", 69);
    noiseRateLayer      = configEnv->GetValue("NoiseRate", 5.0);
    noiseMaxLayer       = configEnv->GetValue("MaxNoise", 20);

    smearZ              = configEnv->GetValue("SmearZ", 0.0120);
    smearRPhi           = configEnv->GetValue("SmearRPhi", 0.0030);
    
    deltaPhiCut         = configEnv->GetValue("DeltaPhiCut", 0.0046875);
    runningWindowSize   = configEnv->GetValue("RunningWindowSize", 0.5);

    inputFileName       = configEnv->GetValue("outputName", "simulation_output.root");

    }

void Config::Print()
{
    cout << "\n=== CONFIGURAZIONE DEL DETECTOR ===" << endl;
    cout << "Geometria:" << endl;
    cout << "  - BeamPipe: R=" << beamPipeRadius << " cm, spessore=" << beamPipeThickness << " cm" << endl;
    cout << "  - Layer1: R=" << layer1Radius << " cm, spessore=" << layer1Thickness << " cm" << endl;
    cout << "  - Layer2: R=" << layer2Radius << " cm, spessore=" << layer2Thickness << " cm" << endl;
    cout << "  - Lunghezza totale: " << detectorLength << " cm" << endl;
    cout << "\nParametri vertice:" << endl;
    cout << "  - sigma_z = " << vertexZSigma << " cm" << endl;
    cout << "  - sigma_xy = " << vertexXYSigma << " cm" << endl;
    cout << "\nRisoluzione detector:" << endl;
    cout << "  - sigma_z = " << smearZ*10000 << " um" << endl;
    cout << "  - sigma_rphi = " << smearRPhi*10000 << " um" << endl;
    cout << "\nSimulazione:" << endl;
    cout << "  - Eventi: " << nEvents << endl;
    cout << "  - Molteplicita': [" << multiplicityMin << ", " << multiplicityMax << "]" << endl;
    cout << "  - Taglio DeltaPhi: " << deltaPhiCut << " rad" << endl;
    cout << "====================================\n" << endl;
}