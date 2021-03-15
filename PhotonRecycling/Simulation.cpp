#include "Simulation.h"
using namespace std;

double Reflectivity(double n1, double n2, double theta)
{
    if (fabs(n1 / n2 * sin(theta)) >= 1)
    {
        return 1;
    }

    double Rs, Rp;
    Rs = pow((n1 * cos(theta) - n2 * sqrt(1 - pow(n1 / n2 * sin(theta), 2))) / (n1 * cos(theta) + n2 * sqrt(1 - pow(n1 / n2 * sin(theta), 2))), 2);
    Rp = pow((n2 * cos(theta) - n1 * sqrt(1 - pow(n1 / n2 * sin(theta), 2))) / (n2 * cos(theta) + n1 * sqrt(1 - pow(n1 / n2 * sin(theta), 2))), 2);

    if (0.5 * (Rs + Rp) > 1)
        cout << Rs << " " << Rp << " " << n1 << " " << theta / M_PI * 180 << " " << 0.5 * (Rs + Rp) << endl;

    return 0.5 * (Rs + Rp);
}

double Simulation::SampleUniformOnSphere(double &theta, double &phi, double thetamin, double thetamax)
{
    theta = -1;
    while (theta <= thetamin || theta >= thetamax)
    {
        double x = My_Rnd.Gaus(0, 1);
        double y = My_Rnd.Gaus(0, 1);
        double z = My_Rnd.Gaus(0, 1);
        double r = sqrt(x * x + y * y + z * z);
        theta = acos(z / r);
    }

    phi = My_Rnd.Uniform(0, 2 * M_PI);
}

Simulation::Simulation()
{
    NSim = 1e6;    // Total number of simulated incident photon
    Thickness = 3; // the thickness of the material, unit is um
    // Assume the initial light strength is Gaussian distributed
    LightCentral = 770;           // nm
    LightResolution = 40 / 2.355; // nm
    LightPosition = Thickness / 2;
    LightPositionResolution = Thickness / 6;
    // The diffusion length and its resolution
    DiffusionLength = 0.5;     // um
    DiffusionResolution = 0.2; // um

    PLQEint = 0.95; // PLQE
    nout = 1;       // outside material's refraction index

    seed = 11111;

    DoPhotonRecycling = false;
    isOutput = false;

    AverageRecycleTime = 0;
    AverageReflectionTime = 0;
    AverageTransverseD = 0;
    Nout = 0;
    NLightOut = 0;
    NEIOut = 0;
}

Simulation::~Simulation()
{
    if (!isOutput)
    {
        delete IncidentLambda;
        delete OutputLambda;
        delete OutputPosition;
        delete RecycleTime;
        delete TransverseD;
        delete ReflectionTime;
    }
}

void Simulation::SetMaterialProperty(const char *nfilename, const char *kapafilename, const char *alphafilename, double PLQE, double thickness)
{
    // Initial The Material
    ifstream infile_n, infile_kapa, infile_alpha;
    infile_n.open(nfilename, ios::in);
    infile_kapa.open(kapafilename, ios::in);
    infile_alpha.open(alphafilename, ios::in);
    double nx[1000], ny[1000], kapax[1000], kapay[1000], alphax[1000], alphay[1000];
    int nN = 0, kapaN = 0, alphaN = 0;

    while (infile_n.peek() != EOF)
    {
        infile_n >> nx[nN] >> ny[nN];
        nN++;
    }

    while (infile_kapa.peek() != EOF)
    {
        infile_kapa >> kapax[kapaN] >> kapay[kapaN];
        kapaN++;
    }

    while (infile_alpha.peek() != EOF)
    {
        infile_alpha >> alphax[alphaN] >> alphay[alphaN];
        alphaN++;
    }

    infile_n.close();
    infile_kapa.close();
    infile_alpha.close();

    Perovskite.Setn(nx, ny, nN);
    Perovskite.Setkapa(kapax, kapay, kapaN);
    Perovskite.Setalpha(alphax, alphay, alphaN);

    PLQEint = PLQE;
    Thickness = thickness;
}

void Simulation::SetInitialLight(double peak, double fwhm)
{
    LightCentral = peak;
    LightResolution = fwhm / 2.355;
    DoPhotonRecycling = false;
}

void Simulation::SetInitialLight(double peak, double fwhm, double position, double reso)
{
    LightCentral = peak;
    LightResolution = fwhm / 2.355;

    // when simulated photon recycling
    LightPosition = position;
    LightPositionResolution = reso;
    DoPhotonRecycling = true;
}

void Simulation::ifPhotonRecycling(bool option)
{
    DoPhotonRecycling = option;
}

void Simulation::SetSeed(int Seed)
{
    My_Rnd.SetSeed(seed);
}

void Simulation::RunSimulation(int times)
{
    NSim = times;

    //Initial variables
    vector<EMObject> OutPhoton; // storage of output photon
                                // Begin the Loop
    IncidentLambda = new TH1D("IncidentLambda", "IncidentLambda", 200, LightCentral - 5 * LightResolution, LightCentral + 5 * LightResolution);
    OutputLambda = new TH1D("OutputLambda", "OutputLambda", 200, LightCentral - 5 * LightResolution, LightCentral + 5 * LightResolution);
    OutputPosition = new TH1D("OutputPosition", "OutputPosition", 200, 0, Thickness);
    RecycleTime = new TH1D("RecycleTime", "RecycleTime", 20, 0, 20);
    TransverseD = new TH1D("TransverseD", "TransverseD", 100, 0, 5);
    ReflectionTime = new TH1D("ReflectionTime", "ReflectionTime", 100, 0, 100);
    IncidentLambda->Sumw2();
    OutputLambda->Sumw2();
    OutputPosition->Sumw2();
    TransverseD->Sumw2();
    ReflectionTime->Sumw2();
    RecycleTime->Sumw2();

    bool doDebug = false;

    for (int i = 0; i < NSim; i++)
    {
        // Firstly, create an incident photon
        EMObject IncidentObject;
        IncidentObject.lambda = My_Rnd.Gaus(LightCentral, LightResolution);
        if (DoPhotonRecycling)
        {
            double tempd;
            if (fabs(LightPosition) < 1e-8)
            {
                IncidentObject.d = Thickness;
                // only lower half of the sphere
                // SampleUniformOnSphere(IncidentObject.theta, IncidentObject.phi, 0.50001 * M_PI, M_PI); // here theta is the angle between its direction and positive z axis
                // simulate photon injection
                SampleUniformOnSphere(IncidentObject.theta, IncidentObject.phi, M_PI - 0.750245244, M_PI); // here theta is the angle between its direction and positive z axis
            }
            else
            {
                tempd = My_Rnd.Gaus(LightPosition, LightPositionResolution);
                while (tempd >= Thickness || tempd <= 0)
                {
                    tempd = My_Rnd.Gaus(LightPosition, LightPositionResolution);
                }
                IncidentObject.d = tempd;
                // On the whole sphere
                SampleUniformOnSphere(IncidentObject.theta, IncidentObject.phi); // here theta is the angle between its direction and positive z axis
            }
        }
        else
        {
            PLQEint = 0;
            IncidentObject.d = Thickness;
            // only lower half of the sphere
            SampleUniformOnSphere(IncidentObject.theta, IncidentObject.phi, 0.50001 * M_PI, M_PI); // here theta is the angle between its direction and positive z axis
        }

        IncidentObject.isPhoton = true;
        IncidentObject.isEIpair = false;
        IncidentObject.isEliminated = false;
        IncidentObject.isOut = false;
        IncidentObject.RecycleTime = 0;
        IncidentObject.ReflectionTime = 0;
        IncidentObject.TransverseX = 0;
        IncidentObject.TransverseY = 0;

        IncidentLambda->Fill(IncidentObject.lambda);

        if (doDebug)
        {
            cout << "Generate Photon : lambda " << IncidentObject.lambda << " Position : " << IncidentObject.d << " Theta : " << IncidentObject.theta << " Phi :" << IncidentObject.phi << endl;
            cout << "n material : " << Perovskite.n(IncidentObject.lambda) << endl;
        }

        while (!IncidentObject.isOut && !IncidentObject.isEliminated)
        {
            if (IncidentObject.isPhoton)
            {
                // generate the penetration length
                // double length = My_Rnd.Exp(IncidentObject.lambda / (4 * M_PI * Perovskite.kapa(IncidentObject.lambda))) * 0.001;
                double length = My_Rnd.Exp(1 / Perovskite.alpha(IncidentObject.lambda));

                // cout << IncidentObject.lambda/(4 * M_PI * Perovskite.kapa(IncidentObject.lambda)) *0.001 << endl;
                // cout << Perovskite.alpha(IncidentObject.lambda) << endl;
                if (doDebug)
                    cout << " Photon d: " << IncidentObject.d << " Theta : " << IncidentObject.theta << " Phi :" << IncidentObject.phi << endl;

                double tempd = IncidentObject.d;
                IncidentObject.d += length * cos(IncidentObject.theta);
                if (doDebug)
                    cout << "Photon Penetration Length : " << length << " Reflectivity : " << Reflectivity(Perovskite.n(IncidentObject.lambda), nout, IncidentObject.theta - M_PI) << endl;
                // here considered multiple reflection
                int tempreflectiontime = 0;
                int LastReflectionPosition = 0;
                while (IncidentObject.d <= 0 || IncidentObject.d >= Thickness)
                {
                    if (IncidentObject.d <= 0)
                    {
                        // if the photon reach the lower surface, decide wheter it reflects or not
                        double possibility = My_Rnd.Rndm() - Reflectivity(Perovskite.n(IncidentObject.lambda), nout, IncidentObject.theta - M_PI);
                        if (possibility > 0)
                        {
                            // the photon goes out
                            if (doDebug)
                                cout << " Photon goes out " << endl;

                            OutputPosition->Fill(tempd);
                            NLightOut++;
                            IncidentObject.isOut = true;
                            IncidentObject.d = 0;
                            break;
                        }
                        else
                        {
                            IncidentObject.d = -IncidentObject.d;
                            IncidentObject.theta = M_PI - IncidentObject.theta;
                            IncidentObject.ReflectionTime++;
                            tempreflectiontime++;
                            LastReflectionPosition = -1;
                            if (tempreflectiontime == 1)
                            {
                                IncidentObject.TransverseX += tempd * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                                IncidentObject.TransverseY += tempd * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                            }
                            else
                            {
                                IncidentObject.TransverseX += Thickness * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                                IncidentObject.TransverseY += Thickness * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                            }

                            if (doDebug)
                            {
                                cout << " Photon reflected on the lower surface" << endl;
                                cout << " Photon Position : " << IncidentObject.TransverseX << " " << IncidentObject.TransverseY << " " << IncidentObject.d << endl;
                            }
                        }
                    }
                    if (IncidentObject.d > Thickness)
                    {
                        // if the photon reach the upper surface, just reflect back and change into EI pair
                        IncidentObject.d = 2 * Thickness - IncidentObject.d;
                        IncidentObject.theta = M_PI - IncidentObject.theta;
                        IncidentObject.ReflectionTime++;
                        tempreflectiontime++;
                        LastReflectionPosition = 1;
                        if (tempreflectiontime == 1)
                        {
                            IncidentObject.TransverseX += (Thickness - tempd) * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                            IncidentObject.TransverseY += (Thickness - tempd) * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                        }
                        else
                        {
                            IncidentObject.TransverseX += Thickness * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                            IncidentObject.TransverseY += Thickness * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                        }

                        if (doDebug)
                        {
                            cout << " Photon reflected on the upper surface " << endl;
                            cout << " Photon Position : " << IncidentObject.TransverseX << " " << IncidentObject.TransverseY << " " << IncidentObject.d << endl;
                        }
                    }
                }

                if (tempreflectiontime >= 1)
                {

                    if (LastReflectionPosition == -1)
                    {
                        IncidentObject.TransverseX += IncidentObject.d * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                        IncidentObject.TransverseY += IncidentObject.d * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                    }
                    else if (LastReflectionPosition == 1)
                    {
                        IncidentObject.TransverseX += (Thickness - IncidentObject.d) * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                        IncidentObject.TransverseY += (Thickness - IncidentObject.d) * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                    }
                }
                else
                {
                    IncidentObject.TransverseX += fabs((tempd - IncidentObject.d)) * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                    IncidentObject.TransverseY += fabs((tempd - IncidentObject.d)) * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                }

                if (doDebug)
                {
                    cout << " Photon Reflection Times : " << tempreflectiontime << endl;
                    cout << " Photon Position : " << IncidentObject.TransverseX << " " << IncidentObject.TransverseY << " " << IncidentObject.d << endl;
                }

                if (IncidentObject.isOut)
                {
                    continue;
                }
                IncidentObject.isPhoton = false;
                IncidentObject.isEIpair = true;
                SampleUniformOnSphere(IncidentObject.theta, IncidentObject.phi);
            }

            if (IncidentObject.isEIpair)
            {
                double RateofElimination = My_Rnd.Rndm() - PLQEint;
                if (RateofElimination > 0)
                {
                    IncidentObject.isEliminated = true;
                    if (doDebug)
                        cout << " Eliminated " << endl;
                    continue;
                }

                if (doDebug)
                    cout << " EIpair d : " << IncidentObject.d << " theta : " << IncidentObject.theta << " phi : " << IncidentObject.phi << endl;

                // generate the penetration length
                double length = My_Rnd.Gaus(DiffusionLength, DiffusionResolution);
                while (length < 0 || length > DiffusionLength + 3 * DiffusionResolution)
                {
                    length = My_Rnd.Gaus(DiffusionLength, DiffusionResolution);
                }

                double tempd = IncidentObject.d;
                IncidentObject.d += length * cos(IncidentObject.theta);
                // change of model, on the surface, the pair just reflects
                // It may be reflected multiple times
                int tempreflectiontime = 0;
                int LastReflectionPosition = 0;

                if (doDebug)
                    cout << " EI pair penetration length " << length << endl;

                while (IncidentObject.d <= 0 || IncidentObject.d > Thickness)
                {
                    if (IncidentObject.d <= 0)
                    {
                        IncidentObject.d = -IncidentObject.d;
                        IncidentObject.ReflectionTime++;
                        tempreflectiontime++;
                        LastReflectionPosition = -1;
                        if (tempreflectiontime == 1)
                        {
                            IncidentObject.TransverseX += tempd * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                            IncidentObject.TransverseY += tempd * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                        }
                        else
                        {
                            IncidentObject.TransverseX += Thickness * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                            IncidentObject.TransverseY += Thickness * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                        }

                        if (doDebug)
                        {
                            // cout << " EIpair reflection on the lower surface " << endl;
                            // cout << " EIpair Position : " << IncidentObject.TransverseX << " " << IncidentObject.TransverseY << " " << IncidentObject.d << endl;
                        }
                    }
                    if (IncidentObject.d > Thickness)
                    {
                        IncidentObject.d = 2 * Thickness - IncidentObject.d;
                        IncidentObject.ReflectionTime++;
                        tempreflectiontime++;
                        LastReflectionPosition = 1;
                        if (tempreflectiontime == 1)
                        {
                            IncidentObject.TransverseX += (Thickness - tempd) * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                            IncidentObject.TransverseY += (Thickness - tempd) * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                        }
                        else
                        {
                            IncidentObject.TransverseX += Thickness * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                            IncidentObject.TransverseY += Thickness * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                        }

                        if (doDebug)
                        {
                            // cout << " EIpair reflection on the upper surface " << endl;
                            // cout << " EIpair Position : " << IncidentObject.TransverseX << " " << IncidentObject.TransverseY << " " << IncidentObject.d << endl;
                        }
                    }
                }

                if (tempreflectiontime >= 1)
                {

                    if (LastReflectionPosition == -1)
                    {
                        IncidentObject.TransverseX += IncidentObject.d * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                        IncidentObject.TransverseY += IncidentObject.d * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                    }
                    else if (LastReflectionPosition == 1)
                    {
                        IncidentObject.TransverseX += (Thickness - IncidentObject.d) * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                        IncidentObject.TransverseY += (Thickness - IncidentObject.d) * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                    }
                }
                else
                {
                    IncidentObject.TransverseX += fabs(tempd - IncidentObject.d) * sin(IncidentObject.theta) * cos(IncidentObject.phi);
                    IncidentObject.TransverseY += fabs(tempd - IncidentObject.d) * sin(IncidentObject.theta) * sin(IncidentObject.phi);
                }

                if (doDebug)
                {
                    cout << " EIpair Reflection Times : " << tempreflectiontime << endl;
                    cout << " EIpair Position : " << IncidentObject.TransverseX << " " << IncidentObject.TransverseY << " " << IncidentObject.d << endl;
                }
                IncidentObject.isPhoton = true;
                IncidentObject.isEIpair = false;
                IncidentObject.lambda = My_Rnd.Gaus(LightCentral, LightResolution);
                SampleUniformOnSphere(IncidentObject.theta, IncidentObject.phi);
                IncidentObject.RecycleTime++;
            }
        }

        if (IncidentObject.isOut)
        {
            // OutPhoton.push_back(IncidentObject);
            OutputLambda->Fill(IncidentObject.lambda);
            RecycleTime->Fill(IncidentObject.RecycleTime);
            TransverseD->Fill(sqrt(IncidentObject.TransverseX * IncidentObject.TransverseX + IncidentObject.TransverseY * IncidentObject.TransverseY));
            ReflectionTime->Fill(IncidentObject.ReflectionTime);
            AverageRecycleTime += IncidentObject.RecycleTime;
            AverageReflectionTime += IncidentObject.ReflectionTime;
            AverageTransverseD += sqrt(IncidentObject.TransverseX * IncidentObject.TransverseX + IncidentObject.TransverseY * IncidentObject.TransverseY);
            Nout++;
            // cout << sqrt(IncidentObject.TransverseX * IncidentObject.TransverseX + IncidentObject.TransverseY * IncidentObject.TransverseY) << endl;
        }

        if (doDebug)
            cout << endl;
    }

    AverageRecycleTime /= (double)Nout;
    AverageReflectionTime /= (double)Nout;
    AverageTransverseD /= (double)Nout;

    RecycleTime->Scale(1 / RecycleTime->Integral());
    TransverseD->Scale(1 / TransverseD->Integral());
    for (int i = 0; i < RecycleTime->GetNbinsX(); i++)
    {
        RecycleTime_Distribution[i] = RecycleTime->GetBinContent(i + 1);
    }

    for (int i = 0; i < TransverseD->GetNbinsX(); i++)
    {
        TransverseD_Distribution[i] = TransverseD->GetBinContent(i + 1);
    }
}

void Simulation::Output(const char *outputfilename)
{
    /*
    cout << " The average recycling time is " << AverageRecycleTime << endl;
    cout << " The average reflection time is " << AverageReflectionTime << endl;
    cout << " N out : " << Nout << endl;
    cout << " N photon out : " << NLightOut << endl;
    cout << " N EI out : " << NEIOut << endl;
    */

    cout << Thickness << " " << Nout << " " << AverageRecycleTime << " " << AverageTransverseD << " " << PLQEint << " " << OutputLambda->GetMean() << endl;

    TFile *outfile = new TFile(outputfilename, "RECREATE");
    outfile->cd();
    IncidentLambda->Write();
    OutputLambda->Write();
    OutputPosition->Write();
    RecycleTime->Write();
    TransverseD->Write();
    ReflectionTime->Write();
    outfile->Close();
    isOutput = true;

    delete IncidentLambda;
    delete OutputLambda;
    delete OutputPosition;
    delete RecycleTime;
    delete TransverseD;
    delete ReflectionTime;
}
