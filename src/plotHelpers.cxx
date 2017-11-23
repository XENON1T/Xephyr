#include "TGraph.h"
#include "TH1F.h"
#include "TTree.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TEventList.h"
#include "TCanvas.h"
#include <iostream>
#include "TDirectory.h"
#include "TFile.h"

using namespace std;

#ifndef PLOT_HELPERS
#define PLOT_HELPERS

namespace plotHelpers
{

//! \brief produces a TGraph of quantiles from branch entry of a tree
TH1F giveQuantiles(TTree *tree, double percent[], double quantiles[], int nquanta, TString var, TString cut = "")
{

    // filling the list of cuts
    tree->Draw(">>list", cut);
    TEventList *list = (TEventList *)gDirectory->Get("list");
    tree->SetEventList(list);
    double min = tree->GetMinimum(var.Data());
    double max = tree->GetMaximum(var.Data());

    // making histo with 10'000 bins, enough to make the computation binning independent
    TH1F distro_histo("distro_histo", var + " :: " + cut, 10000, min, max);
    tree->Draw(var + ">>distro_histo");
    distro_histo.GetQuantiles(nquanta, quantiles, percent);

    // some printing
    cout << "quantiles:";
    for (int i = 0; i < nquanta; i++)
        cout << TString::Format(" %1.2f \t", percent[i]);
    cout << endl
         << "          ";
    for (int i = 0; i < nquanta; i++)
        cout << TString::Format(" %1.3f \t", quantiles[i]);
    cout << endl;

    tree->SetEventList(0); // removing the list
    delete list;

    return distro_histo;
}



TGraphAsymmErrors giveTSquantiles(TTree *tree, double *mu_list, int mu_size, TString OutDir, double wimpMass)
{
    TGraphAsymmErrors *q_mu_90 = new TGraphAsymmErrors(mu_size);

    TString mass = TString::Itoa(wimpMass, 10);

    TFile *out = new TFile(OutDir + "ts_distros_quantiles_m" + mass + ".root", "RECREATE");

    const int nq = 4;
    double xq[nq] = {0.5, 0.88, 0.90, 0.92}; // 2% unc
    double yq[nq] = {0.};                    // quantiles to be filled

    TCanvas *c1 = new TCanvas();
    c1->Print("quantiles_m"+ mass+".pdf[");

    for (int i = 0; i < mu_size; i++)
    {

        TString cut = "q_mu>= -0.01 && mass ==" + mass + " && mu_fit ==" + TString::Format("%1.2f", mu_list[i]);
        TH1F temp = giveQuantiles(tree, xq, yq, nq, "q_mu", cut);

        // storing the distro in file pdf
        temp.Draw();
        c1->Print("quantiles_m"+ mass+".pdf");

        q_mu_90->SetPoint(i, mu_list[i], yq[2]);
        q_mu_90->SetPointError(i, 0., 0., yq[2] - yq[1], yq[3] - yq[2]);
    }

    c1->Print("quantiles_m"+ mass+".pdf]");
    delete c1;

    out->cd();
    q_mu_90->Write("quantiles_m" + mass);
    out->Close();

    return *q_mu_90;
}


}
#endif