#ifndef SignalDef_c
#define SignalDef_c

#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

// This file just describes a set of functions that returns info about
// the path to signal files and the histogram names based on some inputs,
// like for example the mass.

// You can RUN the full SR1 analysis on your signal model jsu by modifying theese three 
// simple functions. 



TString getSignalHistoName(double primaryParameter, int likelihoodType, double other){
    // Returns the histogram name of for the signal as a function of a primary parameter (in case of SR1 is mass)
    // as a function of the volume type [0=egg_volume, 1=U_volume, 2=Wall_volume]
    // as a function of an additional parameter (not used in SR1) 

    // In SR1 For the signal we use the same histogram (anyway the signal shape in S1S2 doesn't change). 
    // The number of event rescaling according to the volume is done aftewards in the code. 
    // You just need to change this line with your histo name.
    // In SR1 the histo name does not depend on any of the parameter, these are just for your convenience.

    // /!\ IMPORTANT:  the histogram must be normalized to: number_of_events / year / tonne
    // (year = 365 days)  
    return  TString("signal_rescaled") ;

}

TString getSignalFile(double primaryParameter, int likelihoodType, double other){
    // Read comments of function above ;)

    // This function just returns the file name (including file path wrt $XEPHYR_DIR)
    // as a function of the mass of the particle
    // other parameter are not used, but are there for your convenience (in case you want to modify)
    
    TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
    TString path  =  "xephyr_examples/SR1Like/data/";
    TString Signal_file = "SignalModel_SR1_GeV_FromRunSourceCombinedFit180420.root";
    TString m = TString::Format("%1.0f", primaryParameter);     // mass
    Signal_file.Insert( Signal_file.Index("GeV"), m );

    Signal_file = xeDir + path + Signal_file;  

    return Signal_file;

}
    

double getSignalUncertainty(double primaryParameter, int likelihoodType, double other ){

    // This is used in the likelihood definition to assign a rate uncertainty to the signal
    // this is a relative rate uncertainty defined at 1sigma, so for example if the returned value
    // would be 0.2  then you would get that the error on the signal rate would be 20% of its rate.

    // In SR1 we read this from a file, you can do it in whatever way you prefer.
    double mass = primaryParameter;
    TString xeDir(gSystem->Getenv("XEPHYR_DIR"));  // get handle on XEPHYR_DIR

    TFile f_uncertainty(xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/lax_1.5.1_egg/models_extended/Acceptance_RunSourceCombined180420_SR1_WIMP.root");
    TGraph *signal_unc = (TGraph*)f_uncertainty.Get("Graph");
    
    return signal_unc->Eval(mass);
}

#endif