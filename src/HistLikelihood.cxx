#include "HistLikelihood.h"


HistLikelihood::~HistLikelihood(){}

HistLikelihood::HistLikelihood(XeRun *r, SignalModel1D *SignalPtr, TString bkgStr, TString DataStr): ProfileLikelihood() , RunComponent()

{ 

if(printLevel > 1) cout << "HistLikelihood::HistLikelihood() - DEBUG : ENTERING" << endl;

   setRun(r);

   setName("Hist Profile likelihood on "+run->getName());

   sig        = SignalPtr;


   bkg        =  (TH1F*) run->getDataFile()->Get(bkgStr);

   obsData    =  (TH1F*) run->getDataFile()->Get(DataStr);


   if(bkg!=NULL)     asimovData =  (TH1F*) bkg->Clone("asimovData") ; //default

   if(obsData!=NULL) data       =  (TH1F*) obsData->Clone("TempData");  // default

   checkHistos();
  
   initialize(); 

   run->setPValue(this);

   run->markAsChanged(); 
  
   run->checkIt();

   checkRunCompatibility();


   if(printLevel > 1) cout << "HistLikelihood::HistLikelihood() - DEBUG : LEAVING" << endl;
   
}

void HistLikelihood::checkHistos(){
	
	if(bkg==NULL || data==NULL || sig==NULL || obsData ==NULL || asimovData ==NULL){
		cout<< "HistLikelihood::checkHistos() - ERROR  at least one histo is bad " << endl;
		exit(100);
	}

	int nBins = bkg->GetXaxis()->GetNbins();
	
	if(sig->getDefaultSignal()->GetNbinsX() != nBins || data->GetNbinsX() != nBins || 
		obsData->GetNbinsX() != nBins || asimovData->GetNbinsX() != nBins){

		cout<< "HistLikelihood::checkHistos() - ERROR one of the histo has different binning" << endl;
		exit(100);
	}
}

void HistLikelihood::initialize(){

	// default is no nuissance parameter
	isLEffTValueFit		=  false;
	isQyTValueFit		=  false;
	isSystBkgTValueFit	=  false;
	isSigAccFit		=  false;
	isStatSignal		=  false;
	isStatBkgTValueFit	=  false;

	bkgScaleFactor          = 1.;
	bkgScaleUnc             = 1.;

 	SigmaAccSignal          = 0.;

       //initializing nuissance parameter
       ls = new SigmaParameter();
       addParameter(ls);
       lkSystBkg = new TSystBkgParameter(run->getNumber());
       addParameter(lkSystBkg);
       le = new TLEffParameter();
       le->setMaximum(1.);
       le->setMinimum(-1.);
//       le->setStep(1.);
       addParameter(le);
       lqy = new TQyParameter();
       lqy->setMaximum(run->getQy()->getMaximumTval());
       lqy->setMinimum(run->getQy()->getMinimumTval()); 
       //lqy->setStep(run->getQy()->getStepSize()); 
       addParameter(lqy);

       // signal acceptance nuissance
       lkacc = new LKParameter(PAR_SIGN_ACCEPTANCE_TVALUE, NUISANCE_PARAMETER, "signalAccRan",0.,0.01,-2.,2.);
       addParameter(lkacc);  

        activateParameter(lkSystBkg,isSystBkgTValueFit);
        activateParameter(le,isLEffTValueFit);
        activateParameter(lqy,isQyTValueFit);
        activateParameter(lkacc,isSigAccFit);

  	for(int b = 0; b < bkg->GetNbinsX() + 2 ; b++){
		
        	lkStatBkgs.push_back(new TGaussParameter(b, run->getNumber()));
        	addParameter(lkStatBkgs[b]);
    		lkStatBkgs[b]->setInitialValue( 0. );
		double Nevent  = bkg->GetBinContent(b+1);
		double rel_unc = sqrt(Nevent) / Nevent + 0.2 ; // stat plus sys 20%
		if(b < bkg->GetNbinsX()) 
			((TGaussParameter*)lkStatBkgs[b])->setUncertainty( rel_unc );
		else 
			((TGaussParameter*)lkStatBkgs[b])->setUncertainty( 0.02 );  // stat uncertainty on signal
//    		lkStatBkgs[b]->setStep( 0.1 );
    		activateParameter(lkStatBkgs[b],isStatBkgTValueFit);
	}



}

double HistLikelihood::getWimpMass() {return sig->getCurrentMass();}



void HistLikelihood::setData(int dataType) { 

   switch (dataType) {

	case DM_CUT_DATA:
		data = obsData;
		break;

	case ASIMOV_DATA:
		data = asimovData;
		break;

	default :
		data = obsData; 
	 	cout<< "HistLikelihood::setData - ERROR - Wrong data Type, set to obsData " << endl;

   }

}


void HistLikelihood::generateAsimov(double mu_prime) {

	asimovData->Reset();

 	double scaleFactorSignal = mu_prime * getSignalMultiplier()  ;
	
  	for(int b=1; b <= bkg->GetNbinsX()  ; b++) {
		
		double nb = bkg->GetBinContent(b) * bkgScaleFactor ;
		double ns = sig->getDefaultSignal()->GetBinContent(b) * scaleFactorSignal;
		double content = nb + ns ;
		asimovData->SetBinContent(b, content);
	}
}



void HistLikelihood::generateToyDataset(double seed, double mu_prime) { cout << "ERROR - HistLikelihood::generateToyDataset is not supported" << endl;}

double HistLikelihood::getSignalDefaultNorm() { return sig->getDefaultNorm();}

double HistLikelihood::getSignalMultiplier() {return sig->getXsecMultiplier() ; }


double HistLikelihood::nSignalPerCm2() { 
 
   return  ( sig->getDefaultSignal()->Integral() / sig->getDefaultNorm() );


}



double HistLikelihood::computeTheLogLikelihood() {

//Retriving Parameter of Interest
  double sigma=getParameterValue(PAR_SIGMA);

//Retriving Nuissace parametercurrent fitting  t_value
  double t_NormBkg 	  = 	0.;
  double t_Leff		  = 	0.;
  double t_Qy 	 	  = 	0.;
  double t_epsilon_Bkg_b  = 	0.;
  double t_sig_acc  	  = 	0.;
  double t_stat_signal[2] = {0.,0.};

  if(isLEffTValueFit)      t_Leff 	    =   getParameterValue(PAR_LEFF_TVALUE);
  if(isQyTValueFit)	   t_Qy		    =   getParameterValue(PAR_QY_TVALUE);
  if(isSystBkgTValueFit)   t_NormBkg   	    =   getParameterValue(PAR_SYST_BKG_TVALUE);
  if(isSigAccFit)          t_sig_acc   	    =   getParameterValue(PAR_SIGN_ACCEPTANCE_TVALUE);
  if(isStatSignal)         t_stat_signal[0] =   getParameterValue(PAR_GAUSS_TVALUE + bkg->GetNbinsX() );
  if(isStatSignal)         t_stat_signal[1] =   getParameterValue(PAR_GAUSS_TVALUE + bkg->GetNbinsX() + 1 );

//Check limit status
  if(std::isnan(t_Leff)){
    cout<<"Fit is unstable for NP Leff "<<endl;
    return VERY_SMALL;
  }
  if(std::isnan(t_Qy)){
    cout<<"Fit is unstable for NP Qy "<<endl;
    return VERY_SMALL;
  }
  if(std::isnan(t_NormBkg)){
    cout<<"Fit is unstable for NP NormBkg"<<endl;
    return VERY_SMALL;
  }

  if(std::isnan(t_sig_acc)){
    cout<<"Fit is unstable for NP signAcc"<<endl;
    return VERY_SMALL;
  }
  

  double LL=0;

//  if(isSystBkgTValueFit)    LL 	= LL - t_NormBkg*t_NormBkg / 2.; 
  if(isLEffTValueFit)       LL 	= LL - t_Leff*t_Leff / 2.; 
  if(isSigAccFit)           LL 	= LL - t_sig_acc*t_sig_acc / 2.; 
  if(isStatSignal)          LL  = LL - t_stat_signal[0]*t_stat_signal[0] / 2. - t_stat_signal[1]*t_stat_signal[1] / 2.;


  // loop over bins 
  for(int b=1; b <= bkg->GetNbinsX()  ; b++) {

  	double N_s           = sigma *  sig->getXsecMultiplier() * sig->getDefaultSignal()->GetBinContent(b);  // normalized to 10^-45 cm^2

	if(isLEffTValueFit)  N_s = N_s * sig->getLeffmultiplier( t_Leff, b );  // apply Leff
 
        if(isSigAccFit)      N_s = N_s * ( 1. + t_sig_acc * SigmaAccSignal) ;  // THIS SHOULD NOT BE DONE LIKE THIS FIXME but getting sigma directly for sys object 

	if(isStatSignal ){
	     if(b < 7) N_s = N_s * ( 1. + t_stat_signal[0] * 0.02) ;  // THIS SHOULD NOT BE DONE LIKE THIS FIXME but getting sigma directly for sys object 
	     else      N_s = N_s * ( 1. + t_stat_signal[1] * 0.02) ;  // THIS SHOULD NOT BE DONE LIKE THIS FIXME but getting sigma directly for sys object 
	}

  	double N_co_obs      = bkg->GetBinContent(b);
 
	// stat term
	double t_stat_Bkg_b  = 0.;
        if(isStatBkgTValueFit) t_stat_Bkg_b  = getParameterValue(PAR_GAUSS_TVALUE+b-1);  // stat parameter start from 0
  	double N_co   = N_co_obs + t_stat_Bkg_b * ( sqrt(N_co_obs) + 0.2 * N_co_obs ) ;  // stat + 20% sys
  
  	double N_b    = ( t_NormBkg * bkgScaleUnc + bkgScaleFactor )* N_co;

  	double N_obs  = data->GetBinContent(b); 



       //Hagar style 
       if(N_s + N_b <= 0.) return VERY_SMALL;
       LL += N_obs*log(N_s + N_b) -1.*( N_s + N_b )  ;
  
       if(isStatBkgTValueFit && N_co <=0.)  return VERY_SMALL;
  
       if(isStatBkgTValueFit)  LL +=  -1.*t_stat_Bkg_b * t_stat_Bkg_b / 2.;
       //if(isStatBkgTValueFit)  LL += N_co_obs*log(N_co) - N_co;
  }

  return LL;
}




bool HistLikelihood::update() { 

  cout << "HistLikelihood:: Update() was called " << endl; 

  if(!ls->inCombinedMode())  ls->setMinuitUnit(1.);

  activateParameter(lkSystBkg,isSystBkgTValueFit);

  activateParameter(le,isLEffTValueFit);

  activateParameter(lqy,isQyTValueFit);

  activateParameter(lkacc,isSigAccFit);

  for(int b=0; b < bkg->GetNbinsX() + 2; b++) {
    //epsilons nuissance parameter (stat. uncertainties)
    lkStatBkgs[b]->setInitialValue( 0. );
    ((TStatBkgParameter*)lkStatBkgs[b])->setStatError( 1. );
    lkStatBkgs[b]->setStep( 0.1 );
    activateParameter(lkStatBkgs[b],isStatBkgTValueFit);

  }

  return true;

}
