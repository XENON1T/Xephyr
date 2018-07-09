#include "XeLikelihoods.h"
#include "XePdfObjects.h"
#include "dataHandler.h"
#include "TSystem.h"
#include <vector>

///---------------     PARAMETERS ------------         ///



pdfLikelihood* getTheLikelihood(){

 TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
 TString inputDir = xeDir + "/xephyr_examples/likelihood1D/data/";

 
 pdfComponent *bkg = new pdfComponent("bkg", inputDir + "modelsAndData.root");

 bkg->autoLoad(); // this loads automatically all shape uncertainty find in the file with prefix "bkg"

 scaleSys *s1 = new scaleSys("Scaling",0.1); // relative scale uncertainty of 10%
   // s1->setMinimum(-1.);
   // s1->setMaximum(1.);   // if you want to limit the parameter to 1 sigma variation (in this case 10% in bkg rate)
   // s1->setType(FREE_PARAMETER); // with this the parameter will not be constrained by a gaussian.

 bkg->addScaleSys(s1);   


 pdfComponent *Signal_model = new pdfComponent( "Signal" , inputDir + "modelsAndData.root" );
 // Signal_model->setScaleFactor(  ); // in case you want to scale to a specific exposure

 scaleSys *s2 = new scaleSys("SignalScale", 0.2); // 20% uncertainty on signal rate, coming for example from theory or acceptances.
 Signal_model->addScaleSys(s2);
 
 // data loading                    // name  // file name                     // tree name 
 dataHandler *data = new dataHandler("data", inputDir + "modelsAndData.root" , "exampleTree"); 
  
  // pdfLikelihood::Constructor(string NameOfLikelihoood, double "val of signal hypothesis parameter" for exampl mass)  
  pdfLikelihood *pl = new pdfLikelihood("LikelihoodExaxmple", 1.);
    
    pl->setExperiment(1);                    // only used for combination
    pl->addBkgPdfComponent(bkg, false);      // adding the background component, the bool tells if the background should be considered in the safeguard
    					     // here no safeguard
    pl->setSignalPdf(Signal_model);       
    pl->setSignalDefaultNorm(1.E-45);	     // setting the cross section of the histogram
    pl->setWithSafeGuard(false);             // no safeguard.

    pl->setDataHandler(data) ;               // adding the data
      
    return pl;



}

