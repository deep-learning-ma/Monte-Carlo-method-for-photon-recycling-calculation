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

    bool doTest = true;

    if (!doTest)
    {
        // TH2D *RecycleTime_Thickness = new TH2D("RecycleTime_Thickness", "RecycleTime_Thickness", 50, 0, 1, 20, 0, 20);
        // TH2D *TransverseD_Thickness = new TH2D("TransverseD_Thickness", "TransverseD_Thickness", 50, 0, 1, 100, 0, 5);
        for (int i = 0; i < 50; i++)
        {
            Simulation SimulationWithPhotonRecycling;
            double tempthick = 0.06 * (i + 1);
            SimulationWithPhotonRecycling.SetMaterialProperty(argv[1], argv[2], argv[3], 0, tempthick);
            SimulationWithPhotonRecycling.SetInitialLight(770, 40, 0, 0);
            SimulationWithPhotonRecycling.SetSeed(i + 1);
            SimulationWithPhotonRecycling.ifPhotonRecycling(true);
            SimulationWithPhotonRecycling.RunSimulation(1e6);
            SimulationWithPhotonRecycling.Output("output.root");
            /*
            for (int j = 0; j < RecycleTime_Thickness->GetNbinsY(); j++)
            {
                RecycleTime_Thickness->SetBinContent(i + 1, j + 1, SimulationWithPhotonRecycling[i].RecycleTime_Distribution[j]);
            }
            for (int j = 0; j < TransverseD_Thickness->GetNbinsY(); j++)
            {
                TransverseD_Thickness->SetBinContent(i + 1, j + 1, SimulationWithPhotonRecycling[i].TransverseD_Distribution[j]);
            }
            */

            // cout << "Simulation " << i << " end" << endl;
        }

        /*
        TFile *outfile = new TFile("output.root", "update");
        outfile->cd();
        RecycleTime_Thickness->Write();
        TransverseD_Thickness->Write();
        outfile->Close();
        */
    }
    else
    {
        Simulation SimulationWithPhotonRecycling;
        SimulationWithPhotonRecycling.SetMaterialProperty(argv[1], argv[2], argv[3], 0, 2.45);
        SimulationWithPhotonRecycling.SetInitialLight(770, 40, 0, 0);
        SimulationWithPhotonRecycling.SetSeed(11111);
        SimulationWithPhotonRecycling.ifPhotonRecycling(true);
        SimulationWithPhotonRecycling.RunSimulation(1e6);
        SimulationWithPhotonRecycling.Output("output.root");
    }

    return 1;
}
