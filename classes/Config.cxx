#include "Config.h"
ClassImp(Config)

Config::Config(SimRandom* srnd, TEnv* configEnv):
    TObject(),
    fSimrand(srnd)
    {

    beamPipeRadius      = configEnv->GetValue("BeamPipeRadius", 3.0);
    beamPipeThickness   = configEnv->GetValue("BeamPipeWidth", 0.08);
    beamPipeMaterial    = configEnv->GetValue("BeamPipeMaterial", "Be");
    layer1Radius        = configEnv->GetValue("Layer1Radius", 4.0);
    layer1Thickness     = configEnv->GetValue("LayerWidth", 0.02);
    layer2Radius        = configEnv->GetValue("Layer2Radius", 7.0);
    layer2Thickness     = configEnv->GetValue("LayerWidth", 0.02);
    layerMaterial       = configEnv->GetValue("LayerMaterial", "Si");
    detectorLength      = configEnv->GetValue("Length", 27.0);
    vertexZSigma        = configEnv->GetValue("SigmaZ", 5.3);
    vertexXYSigma       = configEnv->GetValue("SigmaXY", 0.01);
    vertexZedges        = configEnv->GetValue("Zedges_uniform", 5.3);

    nEvents             = configEnv->GetValue("Events", 10000);
    multiplicityMin     = configEnv->GetValue("Minimum", 1);
    multiplicityMax     = configEnv->GetValue("Maximum", 69);
    noiseRateLayer      = configEnv->GetValue("NoiseRate", 5.0);
    noiseMaxLayer       = configEnv->GetValue("MaxNoise", 20);
    gentypes            = configEnv->GetValue("Generation", "ghp");
    msEnabled           = configEnv->GetValue("MultScattering", false);
    noiseEnabled        = configEnv->GetValue("NoiseEnabled", false);

    smearZ              = configEnv->GetValue("SmearZ", 0.0120);
    smearRPhi           = configEnv->GetValue("SmearRPhi", 0.0030);
    
    deltaPhiCut         = configEnv->GetValue("DeltaPhiCut", 0.0046875);
    runningWindowSize   = configEnv->GetValue("RunningWindowSize", 0.5);

    inputFileName       = configEnv->GetValue("outputName", "simulation_output.root");

    multminZoom         = configEnv->GetValue("MultZoom_Min", 5);
    multmaxZoom         = configEnv->GetValue("MultZoom_Max", 10);

    displayerrfull      = configEnv->GetValue("Residuals_all_mult", false);
    displayerrselect    = configEnv->GetValue("Residuals_select_mult", false);
    errZlimit           = configEnv->GetValue("Residuals_zlim", 0.1);   

    displayefffull      = configEnv->GetValue("Efficiency_mult_allZ", false);
    displayeff1sigma    = configEnv->GetValue("Efficiency_mult_1sigma", false);
    displayeff3sigma    = configEnv->GetValue("Efficiency_mult_3sigma", false);

    displayresfull      = configEnv->GetValue("Resolution_mult_allZ", false);
    displayres1sigma    = configEnv->GetValue("Resolution_mult_1sigma", false);
    displayres3sigma    = configEnv->GetValue("Resolution_mult_3sigma", false);

    displayeffZvrt      = configEnv->GetValue("Efficiency_Zvert", false);
    displayresZvrt      = configEnv->GetValue("Resolution_Zvert", false);

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