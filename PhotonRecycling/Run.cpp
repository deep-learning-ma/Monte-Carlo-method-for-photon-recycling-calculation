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

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage:  " << argv[0] << " ncurve.dat kapacurve.dat alphacurve.dat" << endl;
        return 0;
    }

    Simulation SimulationWithPhotonRecycling;
    SimulationWithPhotonRecycling.SetMaterialProperty(argv[1], argv[2], argv[3], 0.9, 2.45); // n, kapa, alpha, PLQY and thickness of material
    SimulationWithPhotonRecycling.SetInitialLight(770, 40, 0, 0); // initial light wave length, resolution, initial position, position resolution
    SimulationWithPhotonRecycling.SetSeed(11111); // Random seed
    SimulationWithPhotonRecycling.ifPhotonRecycling(true); // simulate with photon recycling or not
    SimulationWithPhotonRecycling.RunSimulation(1e6); // total events
    SimulationWithPhotonRecycling.Output("output.root");

    return 1;
}
