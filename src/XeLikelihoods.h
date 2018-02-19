#ifndef XELIKELIHOODS_H
#define XELIKELIHOODS_H

#include "XeStat.h"
#include <TString.h>
#include "TGraph2D.h"
#include "TNtuple.h"
#include <cmath>
#include <csignal>
#include <iostream>
#include <cstdio>
#include <vector>
#include <utility>      // std::pair, std::make_pair
#include "dataHandler.h"

#include "XePdfObjects.h"
#include "plotHelpers.h"

using namespace std;

/**
 * \class pdfLikelihood
 * \brief implements 2D pdf profiled likelihood fitting.
 * 
 * It expects as inputs pdfObjects. It also implements safeguard.
 */
class pdfLikelihood : public ProfileLikelihood {

  public :

	//FIXME  TODO put some comments

	~pdfLikelihood();
	pdfLikelihood(TString name, double wimpMass);


    void initialize();

	void setData(int dataType);

	double computeTheLogLikelihood();

  	void generateAsimov(double mu_prime) ;

	void generateAndUseAsimov(double mu_prime) {generateAsimov(mu_prime); useAsimovData() ;}

  	void generateToyDataset(double seed, double mu_prime);

  	double getSignalDefaultNorm() { return siganlDefaultNorm; } ;

  	double getSignalMultiplier()  { return XsecMultiplier;    } ; 

	void   setSignalMultiplier(double val) { XsecMultiplier = val; } ; 

	void   setTreeIndex(int index);

	void   setSignalDefaultNorm(double val) { siganlDefaultNorm = val; } ; 

        double getCurrentNs();

	void   checkInputs();

	// add a bkg pdf, if its shape should be considered in safeguard method
	// Safeguarded should be true, default is false.
	void addBkgPdfComponent(pdfComponent *addMe , bool Safeguarded = false); 
	
	void setSignalPdf(pdfComponent *signal)  {signal_component = signal; };

	void setDataHandler(dataHandler *d )     {dmData = d; useDMData(); } ;

	void setAsimovData(dataHandler *d) {asimovData = d;  useAsimovData(); } ;

	int useAsimovData();

	int useDMData();
	
	bool isNegativeAnywhere(TH2F histo); 
	
	void setCalibrationData(dataHandler *d )     {calibrationData = d; } ;

	void setWithSafeGuard(bool b)  {withSafeGuard = b;} ;

	void drawAllOnProjection(bool isS1Projection);

    /** \brief prints a summary of all bkg and signal events with current parameter choice 
	 * 
	*/
	void printEventSummary();

	vector<string> getTrueParamsNames() { return data->getTrueParamsNames(); };
	
	vector<double> getTrueParams()      { return data->getTrueParams(); };

	pdfComponent* getBkgComponent(TString name);

	// Get a handler for plotting comparisons between models
        // used by drawAllOnProjection()
	histoCompare getModelCompare();

	/**
	 * returns (1-epsilon)Fb(x,y) + epsilon*Fs(x,y) with all 
	 * uncertainties interpolated and added.
	 * this term is a full bkg pdf including also the non safeguarded bkgs
	 * TO BE USED IN THE "physics" likelihood
	*/
	TH2F getSafeguardedBkgPdf();

	/**
	 *  returns (1-epsilon)Fb(x,y) + epsilon*Fs(x,y) for all bkg components that
	 * are considered for dafeguard, with all uncertainties 
	 * interpolated and added. TO BE USED in the fit to calibration
	 */
	TH2F getSafeguardedBkgPdfOnly();

	vector <pdfComponent*> bkg_components;
        
	//! \brief says if a bkg component needs to be safeguarded.
	vector <bool>    safeguarded_bkg_components; 

	// is an additional component to be added in the fit to calibration for
	// safeguard. Used in this case for the AC that is different between data and Rn
	// this component has to be scaled to Bkg_data ---> Nbkg/Ncal
	TH2F *safeguardAdditionalComponent;

	void setAdditionalSafeGuardComponent(TH2F *h){ safeguardAdditionalComponent = h;};

	double LLsafeGuard();	

	//! \brief returns the number of bkg components that ask for safeguard
	int   numberOfSafeguarded();	
	
	pdfComponent           *signal_component;

	dataHandler	       		*data;
	
	dataHandler            *calibrationData;	

	dataHandler            *asimovData;

	dataHandler            *dmData;
	
	double                 siganlDefaultNorm;   //corresponding cross section of histos

	double                 XsecMultiplier;   // handy multiplier to keep sigma [0,50]
	bool                   withSafeGuard;

	double                 wimp_mass;


	LKParameter            *safeGuardParam;

	double                  safeguard_scaling;

	//This is needed for compatibility, ancestral xephyr roots.
	//FIXME: move getWimpMass to Asymptotics
  	double getWimpMass() {return wimp_mass;};
     //------------------------------//



};





#endif
