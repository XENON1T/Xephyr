#include "plotHelpers.h"

namespace plotHelpers {
//! \brief produces a TGraph of quantiles from branch entry of a tree
TH1F giveQuantiles(TTree *tree, double percent[], double quantiles[], int nquanta, TString var, TString cut)
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
    c1->Print(OutDir + "quantiles_m"+ mass+".pdf[");

    for (int i = 0; i < mu_size; i++)
    {
	std::cout << "------->  mu " << mu_list[i] << std::endl;

        TString cut = "q_mu>= -0.01 && mass ==" + mass + " && mu_fit ==" + TString::Format("%1.2f", mu_list[i]);
        TH1F temp = giveQuantiles(tree, xq, yq, nq, "q_mu", cut);

        // storing the distro in file pdf
        temp.Rebin(100);
        temp.Draw();
        c1->Print(OutDir + "quantiles_m"+ mass+".pdf");

        q_mu_90->SetPoint(i, mu_list[i], yq[2]);
        q_mu_90->SetPointError(i, 0., 0., yq[2] - yq[1], yq[3] - yq[2]);
    }

    c1->Print(OutDir + "quantiles_m"+ mass+".pdf]");
    delete c1;

    out->cd();
    q_mu_90->Write("quantiles_m" + mass);
    out->Close();

    return *q_mu_90;
}


TGraphAsymmErrors sensitivity(TTree *tree, TString OutDir, double wimpMass[], int N_mass){
                         // -2       -1        0    +1        +2        sigmas
    double percents[5] = {0.022750, 0.158655, 0.5, 0.841345, 0.977250 };
    double quantiles[5] = { 0. };

    TFile *fOut = new TFile(OutDir + "sensitivity.root", "RECREATE");
    
    TGraphAsymmErrors median_and_one_sigma(N_mass);
    TGraphAsymmErrors median_and_two_sigma(N_mass);
    
    TGraphAsymmErrors XS_median_and_one_sigma(N_mass);
    TGraphAsymmErrors XS_median_and_two_sigma(N_mass);

    if(!OutDir.Contains("/"))
        OutDir += "/";

 
    TCanvas *c1 = new TCanvas();
    c1->Print(OutDir + "limitDistros.pdf[");

    for (int massItr=0; massItr < N_mass; massItr++){
        TString mass = TString::Format("%1.0f", wimpMass[massItr]);
    
        std::cout << "------->  mass " << wimpMass[massItr] << std::endl;

        // limit distro in NS
        TH1F temp = giveQuantiles(tree,percents, quantiles, 5, "mu_limit", "mass ==" + mass);

        // some silly drawing of distros...
        temp.Rebin(100);
        temp.Draw();
        c1->Print(OutDir + "limitDistros.pdf");
        

        // Filling Limit TGraph
        median_and_one_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        median_and_one_sigma.SetPointEYhigh(massItr, quantiles[3]); // + 1 sigma
        median_and_one_sigma.SetPointEYlow(massItr, quantiles[1]);  // - 1 sigma

        median_and_two_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        median_and_two_sigma.SetPointEYhigh(massItr, quantiles[4]); // + 2 sigma
        median_and_two_sigma.SetPointEYlow(massItr, quantiles[0]);  // - 2 sigma

        // Limit distro in X section
        temp = giveQuantiles(tree,percents, quantiles, 5, "limit", "mass ==" + mass);
        
        // some silly drawing of distros...
        temp.Rebin(100);
        temp.Draw();
        c1->Print(OutDir + "limitDistros.pdf");
        
        // Filling Limit TGraph
        XS_median_and_one_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        XS_median_and_one_sigma.SetPointEYhigh(massItr, quantiles[3]); // + 1 sigma
        XS_median_and_one_sigma.SetPointEYlow(massItr, quantiles[1]);  // - 1 sigma

        XS_median_and_two_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        XS_median_and_two_sigma.SetPointEYhigh(massItr, quantiles[4]); // + 2 sigma
        XS_median_and_two_sigma.SetPointEYlow(massItr, quantiles[0]);  // - 2 sigma
    }

    fOut->cd();
    median_and_one_sigma.Write("sensitivity_ns_1s");
    median_and_two_sigma.Write("sensitivity_ns_2s");

    XS_median_and_one_sigma.Write("sensitivity_xsec_1s");
    XS_median_and_two_sigma.Write("sensitivity_xsec_2s");
    fOut->Close();
    
    c1->Print(OutDir + "limitDistros.pdf]");
    delete c1;

    return XS_median_and_one_sigma;
}


void addHisto(TH2F *histo, TH2F *h_toBeAdded, double scalefactor ){
      int Nx = histo->GetNbinsX();
      int Ny = histo->GetNbinsY();
      for (int x=1; x <= Nx ; x++){
        for (int y=1; y <= Ny; y++){
          double new_content = histo->GetBinContent(x,y)  + h_toBeAdded->GetBinContent(x,y) * scalefactor;
         histo->SetBinContent(x,y, new_content ) ;
        }
      }
}


}
