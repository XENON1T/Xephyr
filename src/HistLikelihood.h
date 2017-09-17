#ifndef HISTLIKELIHOOD_H
#define HISTLIKELIHOOD_H

#include "XeCore.h"
#include "XePhys.h"
#include "XeStat.h"
#include "XeAnalysis.h"
#include "TF1.h"
#include <cmath>
#include <csignal>
#include <iostream>
#include "TGraph2D.h"

class HistLikelihood :virtual public ProfileLikelihood, virtual public RunComponent {

  public :
	~HistLikelihood();

	HistLikelihood(XeRun *run, SignalModel1D *SignalPtr, TString bkgStr, TString DataStr);

  	double getWimpMass();

  	void setData(int dataType);

  	void generateAsimov(double mu_prime) ;

  	void generateToyDataset(double seed, double mu_prime);

  	double getSignalDefaultNorm();

  	double getSignalMultiplier(); 

	void   setSignalMultiplier(double val) { sig->setXsecMultiplier(val); };

	bool update();
	
	double nSignalPerCm2();

  	double computeTheLogLikelihood();

        void initialize();
	
 /**
     * set the flag to fit LEff t-value
     * @param  fit  Flag specifying whether LEff should be considered as a nuisance parameter
 */
       void fitLEffTValue(bool fit) {isLEffTValueFit = fit;} ;
 
 /**
     * set the flag to fit Qy t-value
     * @param  fit  Flag specifying whether Qy should be considered as a nuisance parameter
 */
       void fitQyTValue(bool fit) {isQyTValueFit = fit;};

 /**
     * set the flag to fit the systematic background t-value
     * @param  fit  Flag specifying whether systematic Bkg should be considered as a nuisance parameter
 */
    void fitSystBkgTValue(bool fit){isSystBkgTValueFit = fit;};

 /**
     * set the flag to fit the systematic background t-value
     * @param  fit  Flag specifying whether systematic Bkg should be considered as a nuisance parameter
 */
    void fitStatBkgTValue(bool fit){isStatBkgTValueFit = fit;};

    void fitStatSignal(bool fit){ isStatSignal = fit;};

    void   fitSignalAcc(bool doFit) {isSigAccFit = doFit;};

    void   setBkgScaleFactor(double bkgSF) { bkgScaleFactor = bkgSF;};
    void   setBkgScaleFactorUnc(double unc) { bkgScaleUnc   = unc;};

    void   setSigmaAccSignal(double set) { SigmaAccSignal = set; };
    double getSigmaAccSignal() {return SigmaAccSignal;};

        // inputs 
	SignalModel1D *sig;

	TH1F          *bkg;
	
	TH1F	      *data;

	TH1F          *obsData;

	TH1F 	      *asimovData;

   
    protected:

	//flags	  
	 bool    isLEffTValueFit;               /*!<Do we want to fit LEff t-value?*/
	 bool    isQyTValueFit;                 /*!<Do we want to fit Qy   t-value?*/
	 bool    isSystBkgTValueFit; /*!<Do we want to fit background syst t-value?*/
	 bool    isSigAccFit; /*!<Do we want to fit background syst t-value?*/
	 bool    isStatBkgTValueFit; /*!<Do we want to fit background stat t-value?*/
	 bool    isStatSignal; /*!<Do we want to fit background stat t-value?*/
	 double  SigmaAccSignal;

	
	//Handy for nuissance parameter
    	 LKParameter     *le;             /*!< pointer to likelihood LEff parameter*/
	 LKParameter     *lqy;	     /*!<pointer to likelihood Qy parameter*/
    	 LKParameter     *ls;             /*!<pointer to likelihood sigma parameter*/
    	 LKParameter     *lkSystBkg;      /*!< pointer to likelihood background parameter*/
    	 LKParameter     *lkacc;         /*!< pointer to likelihood signal acceptance parameter*/
    	 vector <LKParameter*>     lkStatBkgs;
              /*!< pointers to likelihood background in band t-value parameter*/
	

	double 		bkgScaleFactor;
	double 		bkgScaleUnc;

	void checkHistos();

	ClassDef(HistLikelihood,1)

};








#endif
