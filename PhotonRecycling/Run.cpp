// Root Headers
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

// c++ Headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "Simulation.h"
#include "Material.h"

using namespace std;

double Thickness_PLQY(double thick)
{
    double ThicknessTest[6] = {0.02, 0.09, 0.335, 0.895, 1.463, 2.45};
    double PLQYTest[6] = {0.34677, 0.28182, 0.68065, 0.83561, 0.85404, 0.85302};

    for (int i = 0; i < 5; i++)
    {
        if (thick > ThicknessTest[i] && thick < ThicknessTest[i + 1])
        {
            double k = (PLQYTest[i] - PLQYTest[i + 1]) / (ThicknessTest[i] - ThicknessTest[i + 1]);
            double b = (PLQYTest[i] * ThicknessTest[i + 1] - PLQYTest[i + 1] * ThicknessTest[i]) / (ThicknessTest[i + 1] - ThicknessTest[i]);
            return k * thick + b;
        }
    }

    if (thick > 2.44)
    {
        return 0.8535;
    }

    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage:  " << argv[0] << " ncurve.dat kapacurve.dat alphacurve.dat" << endl;
        return 0;
    }

    bool doTest = false;

    if (!doTest)
    {
        for (int i = 0; i < 50; i++)
        {
            Simulation SimulationWithPhotonRecycling;
            double tempthick = 0.06 * (i + 1);
            SimulationWithPhotonRecycling.SetMaterialProperty(argv[1], argv[2], argv[3], Thickness_PLQY(tempthick), tempthick);
            SimulationWithPhotonRecycling.SetInitialLight(770, 40, 0, 0);
            SimulationWithPhotonRecycling.SetSeed(i + 1);
            SimulationWithPhotonRecycling.ifPhotonRecycling(true);
            SimulationWithPhotonRecycling.RunSimulation(1e6);
            SimulationWithPhotonRecycling.Output("output.root");
        }
    }
    else
    {
        Simulation SimulationWithPhotonRecycling;
        SimulationWithPhotonRecycling.SetMaterialProperty(argv[1], argv[2], argv[3], 0, 0.06);
        SimulationWithPhotonRecycling.SetInitialLight(770, 40, 0.03, 0);
        SimulationWithPhotonRecycling.SetSeed(11111);
        SimulationWithPhotonRecycling.ifPhotonRecycling(true);
        SimulationWithPhotonRecycling.RunSimulation(1e6);
        SimulationWithPhotonRecycling.Output("output.root");
    }

    return 1;
}
