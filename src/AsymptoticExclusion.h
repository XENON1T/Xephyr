#ifndef AsymptoticExclusion_h
#define AsymptoticExclusion_h

#include "XeStat.h"
#include "XePdfObjects.h"
#include "TGraphAsymmErrors.h"
#include "XeUtils.h"

using namespace ROOT;
using namespace Math;

/**
  * Computes limits with the asymptotic formulae, default is CLs.
  * Output: root files with histograms for a single mass point, several files can be hadd togheter. This is to promote parallelization.
  * The mass point has to be defined previously using the appropriated XeRun method
*/
class AsymptoticExclusion : public XeObject{
    public:

	~AsymptoticExclusion();
	AsymptoticExclusion(ProfileLikelihood *pl, double confidenceLevel);

	ProfileLikelihood* getProfileLikelihood() {return plike;};
	TGraphAsymmErrors  getExpectedMedian()    {return expectedMedian;};
	TGraphAsymmErrors  getExpected2Sigma()    {return expected2Sigma;};
        TGraph             getObservedLimit()	  {return observedLimit;};
        double 		   getMass()		  {return plike->getWimpMass();};
	TString		   getName() 		  {return name;};

    /**
      * This will trigger the XeRun associated to produce an Asimov dataset with the 
        specified signal strenght mu_prime. The generated asimov is set as current data for 
	the profile likelihood. 
    */
        void generateAndSetAsimov(double mu_prime);

    /**
      * Set DM data as current data for the profile likelihood.
    */
	void setRealData();

    /**
      * This will trigger the XeRun associated to produce an a toy dataset of bkg + signal
	with the specified signal strenght mu_prime. The generated dataset is set as current data for 
	the profile likelihood. 
    */
        void setToyDataset(double seed, double mu_prime);
 
    /**
      * Compute the value of the test statistic q_tilde for a given signal strenght mu.
      * @param useStoredFit impose to use the conditional fit outcome already computed,
	usefull when computing limits with data.
    */
	double computeQTestStat(double mu, bool useStoredFit = false);

    /**
      * Compute the expected CLs limits with bands. It stores it in the graphs getExpectedMedian ecc..
    */
	void computeSensitivity();

     /**better way of producing sensitivity **/
	double computeSensitivityHagar();


    /**
      * Limit wrapper formula, this is based on the "asymptotic paper" formula 89, and includes CLs.
      * returns the limit for the specified band (zero in case of median) in terms of the signal 
      strenght.
      * @param sigma, is the asimov sigma related to mu =0
      * @param CL, confidence level
      * @param N, sigma band, can take positive and negative value
    */
	double computeExpectedLimit(double sigma_0, double N, double CL);

    /**
      * returns the log(likelihood(mu,theta)) for the conditional fit with signal strenght mu
    */
	double conditionalFit(double mu);

    /**
      * returns the -2log(likelihood(mu,theta)) for the unconditional fit and saves mu_hat best fit
	to an internal variable retrivable with "getMu_hat()", it aslo store the outcome  of the fit
	into "logD" member.
    */
	double unconditionalFit();

    /** 
      * returns signal strenght best fit, if the fit has not been runned returns -9999
    */ 	
        double getMu_hat();

    /**
      * Compute the variance of sigma_hut around sigma_true with Asimov dataset
      * @param mu_prime, true assumed value of signal strenght
    */
	double computeSigmaAsimov(double mu_prime);


    /**
      * Set the number of scan points for the limit compyutation, default is 100.
    */
       void setNscanPoints(int nPoints){nScanPoints = nPoints;};



    /**
      * Write to file the histograms
    */
        void writeToFile(TString prefix ="limit_");

        TGraph getSigmaScan(){return sigmaScan;};

	TGraph getqTestScan() {return qTestScan;};
	

     /**
       * Compute Observed limit with CLs and asymptotics
     */

       void computeLimits();

      /**
	* Compute pvalues
      */
	double compute_pval_s_plus_b(double q_obs);
	double compute_pval_b(double q_obs, double mu_val, double sigma_val);

        void setScanMin(double mi) {scanMin = mi;};
        void setScanMax(double ma) {scanMax = ma;};

	void LikelihoodScan();

	void setQTilde(bool doSet) {useQtilde = doSet;};  

       /**
    	* Will plot the limits Graphs as a function of another variable. A 2D Histo will be also 
	filled  with the limit reported for mass VS additionalAxis that can be plotted as a COLZ
	**/
	void setAlternativeXAxisVal(double x) { AlternativeX = x;};

	void setAlternative2DHistoRange(int Nx, double xmin, double xmax, int Ny, double ymin, double ymax);

	TGraph  getqTestScanData() {return qTestScanData;};

	double getObsLimitCLS(){return Obslimit;}
	double getObsLimitnoCLS(){return ObslimitnoCLS;}
	

    private:
    /**
      * compute set of asimov sigma for limit computation of Pb
    */
        void computeSetAsimovSigma();

	ProfileLikelihood *plike;		// handler on profile likelihood

	double Obslimit, ObslimitnoCLS;
		
	TGraphAsymmErrors    expectedMedian;	//Graph holding the informations
	TGraphAsymmErrors    expected2Sigma;
	TGraph		     observedLimit;
	TGraph		     observedLimitNoCLS;
	TGraph		     sigmaScan;
	TGraph		     qTestScan;
	TGraph		     qTestScanData;
	TH1F		     pulls_uncond;
	TH1F		     pulls_cond;
	TH2F		     Histo2D_mass_vs_x;   //2D histo of the observed limits as a 
						 //function of mass and the additional variable X, 
						//specified with setAlternativeXAxisVal()

	int		     nScanPoints;	// Number of points for which mu is computed in limit calculation, default 100.

	TString 	     name;

	vector<double>       asimovSigmaSet;
	vector<double>       muStepsSet;

	double 		     cl;			// confidence level DEFAULT = 0.1
	double 		     stored3Sigma_mu;
	double 		     storedMedian_mu;

	double 		     scanMin;
	double 		     scanMax;	

	bool 		     useQtilde;

	double               AlternativeX;

};


#endif
