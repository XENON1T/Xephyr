#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TEventList.h"
#include "TCanvas.h"
#include <iostream>
#include "TDirectory.h"
#include "TFile.h"
#include "XeLikelihoods.h"

using namespace std;

#ifndef PLOT_HELPERS
#define PLOT_HELPERS

namespace plotHelpers
{

//! \brief produces a TGraph of quantiles from branch entry of a tree
TH1F giveQuantiles(TTree *tree, double percent[], double quantiles[], int nquanta, TString var, TString cut = "");


TGraphAsymmErrors giveTSquantiles(TTree *tree, double *mu_list, int mu_size, TString OutDir, double wimpMass);


TGraphAsymmErrors sensitivity(TTree *tree, TString OutDir, double wimpMass[], int N_mass);


void addHisto(TH2F *histo, TH2F *h_toBeAdded, double scalefactor );

}
#endif
