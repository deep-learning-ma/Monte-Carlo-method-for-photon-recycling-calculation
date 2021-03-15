#ifndef SIMULATION_HEADER
#define SIMULATION_HEADER
// Root Headers
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TRandom3.h>

// c++ Headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "Material.h"

#define M_PI 3.14159265358979323846

struct EMObject
{
    double lambda;
    double d;
    double theta;
    double phi;
    bool isPhoton;
    bool isEIpair;
    bool isEliminated;
    bool isOut;
    int RecycleTime;
    int ReflectionTime;
    double TransverseX;
    double TransverseY;
};

double Reflectivity(double n1, double n2, double theta);

class Simulation
{
public:
    Material Perovskite;
    int NSim;
    double Thickness;
    double LightCentral;
    double LightResolution;
    double LightPosition;
    double LightPositionResolution;
    double DiffusionLength;
    double DiffusionResolution;
    double PLQEint;
    double nout;
    bool DoPhotonRecycling;
    int seed;

    // process variable
    double AverageRecycleTime;
    double AverageReflectionTime;
    double AverageTransverseD;
    int Nout;
    int NLightOut;
    int NEIOut;

    TH1D *IncidentLambda;
    TH1D *OutputLambda;
    TH1D *OutputPosition;
    TH1D *RecycleTime;
    TH1D *TransverseD;
    TH1D *ReflectionTime;

    double RecycleTime_Distribution[1000];
    double TransverseD_Distribution[1000];

    bool isOutput;

    TRandom3 My_Rnd;

    Simulation();
    virtual ~Simulation();
    void SetMaterialProperty(const char *nfilename, const char *kapafilename, const char *alphafilename, double PLQE, double thickness);
    void SetInitialLight(double peak, double fwhm);
    void SetInitialLight(double peak, double fwhm, double position, double reso);
    void SetSeed(int Seed);
    void ifPhotonRecycling(bool option);
    void RunSimulation(int times);
    void Output(const char *outputfilename);
    double SampleUniformOnSphere(double &theta, double &phi, double thetamin = 0, double thetamax = M_PI);
};

#endif