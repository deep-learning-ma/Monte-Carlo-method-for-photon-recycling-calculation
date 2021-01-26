void Plot()
{
    TFile *infile = new TFile("output.root", "READ");
    TH1D *IncidentLambda, *OutputLambda, *OutputPosition;
    TH2D *RecycleTime_Thickness;
    TH2D *TransverseD_Thickness;

    IncidentLambda = (TH1D *)infile->Get("IncidentLambda");
    OutputLambda = (TH1D *)infile->Get("OutputLambda");
    OutputPosition = (TH1D *)infile->Get("OutputPosition");
    RecycleTime_Thickness = (TH2D *)infile->Get("RecycleTime_Thickness");
    TransverseD_Thickness = (TH2D *)infile->Get("TransverseD_Thickness");

    /*
    TCanvas *c = new TCanvas("c", "canvas", 1000, 800);
    c->cd();
    TPad *pad = new TPad("pad", "pad", 0, 0, 1, 1);
    pad->Draw();
    pad->cd();
    IncidentLambda->SetTitle("Intensity");
    IncidentLambda->SetStats(0);
    IncidentLambda->Draw("E0");
    OutputLambda->SetStats(0);
    OutputLambda->SetLineColor(2);
    OutputLambda->Draw("E0 same");

    TLegend *legend = new TLegend(0.70, 0.70, 0.9, 0.9);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetLineWidth(0);
    legend->AddEntry(IncidentLambda, "Incident", "lpfe");
    legend->AddEntry(OutputLambda, "Output", "lpfe");
    legend->Draw("same");
    */

    ofstream outfile;
    outfile.open("WithPR.dat", ios::out);

    for (int i = 0; i < TransverseD_Thickness->GetNbinsX(); i++)
    {
        for (int j = 0; j < TransverseD_Thickness->GetNbinsY(); j++)
        {
            outfile << TransverseD_Thickness->GetBinContent(i + 1, j + 1) << " ";
        }
        outfile << endl;
    }
}