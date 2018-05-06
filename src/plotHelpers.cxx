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
        if(quantiles[i] > 1.E-3) cout << TString::Format(" %1.3f \t", quantiles[i]);
        else cout << TString::Format(" %.2e \t", quantiles[i]);
    cout << endl;

    tree->SetEventList(0); // removing the list
    delete list;

    return distro_histo;
}

double givePval(TTree *tree, double value, TString var, TString cut){

    // filling the list of cuts
    tree->Draw(">>list", cut);
    TEventList *list = (TEventList *)gDirectory->Get("list");
    tree->SetEventList(list);
    double min = tree->GetMinimum(var.Data());
    double max = tree->GetMaximum(var.Data());

    TCanvas c1;
    // making histo with 10'000 bins, enough to make the computation binning independent
    TH1F distro_histo("distro_histo", var + " :: " + cut, 1000000, min, max);
    tree->Draw(var + ">>distro_histo");
    int bin = distro_histo.FindBin(value);
    double p_val = distro_histo.Integral(bin, -1) / distro_histo.Integral();
    c1.Print("pval_distro_example.png");

    return p_val;

}

TGraphAsymmErrors giveTSquantiles(TTree *tree, double *mu_list, int mu_size, TString OutDir, double wimpMass)
{
    TGraphAsymmErrors *q_mu_90 = new TGraphAsymmErrors(mu_size + 1);

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

        if(i == 0 )  q_mu_90->SetPoint(0, 0, yq[2] ); // just to set mu =0 to the closest computed value
        q_mu_90->SetPoint(i +1, mu_list[i], yq[2]);
        q_mu_90->SetPointError(i +1, 0., 0., yq[2] - yq[1], yq[3] - yq[2]);
    }

    c1->Print(OutDir + "quantiles_m"+ mass+".pdf]");
    delete c1;

    out->cd();
    q_mu_90->GetXaxis()->SetTitle("#mu_{test} [similar to #events]");
    q_mu_90->GetYaxis()->SetTitle("Alternative Hypo LLR 90\% quantile");
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
        median_and_one_sigma.SetPointEYhigh(massItr, quantiles[3] - quantiles[2]); // + 1 sigma
        median_and_one_sigma.SetPointEYlow(massItr, quantiles[2] - quantiles[1]);  // - 1 sigma

        median_and_two_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        median_and_two_sigma.SetPointEYhigh(massItr, quantiles[4] - quantiles[2] ); // + 2 sigma
        median_and_two_sigma.SetPointEYlow(massItr, quantiles[2] - quantiles[0]);  // - 2 sigma

        // Limit distro in X section
        temp = giveQuantiles(tree,percents, quantiles, 5, "limit", "mass ==" + mass);
        
        // some silly drawing of distros...
        temp.Rebin(100);
        temp.Draw();
        c1->Print(OutDir + "limitDistros.pdf");
        
        // Filling Limit TGraph
        XS_median_and_one_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        XS_median_and_one_sigma.SetPointEYhigh(massItr, quantiles[3] - quantiles[2]); // + 1 sigma
        XS_median_and_one_sigma.SetPointEYlow(massItr, quantiles[2] - quantiles[1]);  // - 1 sigma

        XS_median_and_two_sigma.SetPoint(massItr, wimpMass[massItr], quantiles[2]);
        XS_median_and_two_sigma.SetPointEYhigh(massItr, quantiles[4] - quantiles[2]); // + 2 sigma
        XS_median_and_two_sigma.SetPointEYlow(massItr, quantiles[2] - quantiles[0]);  // - 2 sigma
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

TGraphAsymmErrors pulls(TTree *tree, TString OutDir, TString name, TString cut="",  bool isUnCond =true ){
        

        double percents[3] = { 0.158655, 0.5, 0.841345 };
        double quantiles[3] = { 0. };
        int numberOfParams = 0;
        vector <string> *name_params = 0;
        tree->SetBranchAddress("name_params",&name_params);
        
        tree->GetEntry(1);
        cout << "numberOfParams " << numberOfParams << "  " << name_params->at(0) << endl;
        unsigned int Nparam = name_params->size();

        TGraphAsymmErrors pulls(Nparam);
        TGraphAsymmErrors sigma1(Nparam);

        TCanvas *c1 = new TCanvas();
        c1->Print(OutDir + name +"_pulls.pdf[");
    
        for(unsigned int i =0; i < Nparam; i++) {
            TString var = "cond_params[" + TString::Itoa(i, 10) + "]";
            if(isUnCond)  var = "un" + var;

            TH1F temp = giveQuantiles(tree,percents, quantiles, 3, var, cut);
            if(!isUnCond) temp.SetTitle(TString(name_params->at(i)) + " Conditional Fit and " + cut) ;
            else temp.SetTitle(TString(name_params->at(i)) + " Unconditional Fit and " + cut) ;

            temp.Rebin(100);
            temp.Draw();
            c1->Print(OutDir + name +"_pulls.pdf");
            c1->Print(OutDir + name + TString(name_params->at(i) + ".png"));
            pulls.SetPoint(i, i + 2, quantiles[1]);
            pulls.SetPointEYhigh(i, quantiles[2] - quantiles[1]);
            pulls.SetPointEYlow(i,  quantiles[1] - quantiles[0] );
            sigma1.SetPoint(i, i+2, 0.);
            sigma1.SetPointEYhigh(i, 1);
            sigma1.SetPointEYlow(i, 1);

        }

            TAxis *x = pulls.GetXaxis();
            TAxis *x2 = sigma1.GetXaxis();
            for(unsigned int i=0; i  < Nparam; i++) {
                x->SetBinLabel(x->FindBin(i +2), name_params->at(i).c_str() );
                x2->SetBinLabel(x->FindBin(i +2), name_params->at(i).c_str() );
            }

            pulls.GetYaxis()->SetRangeUser(-2,2);
            pulls.SetTitle("");
            c1->Print(OutDir + name +"_pulls.pdf]");
        TFile output(OutDir + name +"_pulls.root", "RECREATE");
        pulls.Write(name+"_pulls");
        sigma1.Write("one_sigma");
        output.Close();
        return pulls;
        // put file ----> 
        // TH1F giveQuantiles(TTree *tree, double percent[], double quantiles[], int nquanta, TString var, TString cut)
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


double giveBackgroundLikelihood(pdfLikelihood *likeHood, double s1, double s2) {
	  	  
    //just sum up components
    TH2F bkgPdf = likeHood->bkg_components[0]->getInterpolatedHisto();

    for(unsigned int k=1; k < likeHood->bkg_components.size(); k++){

	    TH2F temp_bkgPdf (likeHood->bkg_components[k]->getInterpolatedHisto());
	  	bkgPdf.Add(&temp_bkgPdf);
		//cout << TString::Format("\t%s  = %f events", temp_bkgPdf.GetName(), temp_bkgPdf.Integral()) << endl;
	}

    return bkgPdf.GetBinContent(bkgPdf.FindBin(s1,s2));

}

double giveSignalLikelihood(pdfLikelihood *likeHood, double s1, double s2) {

    TH2F signalPdf(likeHood->signal_component->getInterpolatedHisto());
    double sigma = likeHood->getPOI()->getCurrentValue();

    return signalPdf.GetBinContent(signalPdf.FindBin(s1,s2)) * sigma * likeHood->getSignalMultiplier();
}


}
