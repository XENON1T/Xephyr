#ifndef XeTest_h
#define XeTest_h
#include "XeAnalysis.h"

  enum EbLKTestParams{FB,X0};

/**
 * Example of ML fit of the parameters of a distribution made of
 *    an exponential distribution plus a flat background
 */
class ebLKTest : public LikelihoodFromDataSet {

  protected:
    double bMax;

  public:

   ~ebLKTest();
    ebLKTest(DataSet* da);
    double computeIndividualLog(double *val);

    ClassDef(ebLKTest,1)

};

/**
 * Test and Demonstration package
 */
class XeTest : virtual public S1S2Display, virtual public XeStat
            , virtual public XeCore {

  public :

 ~XeTest();
  XeTest();


 /**
     * Plot rate for Xenon target, spin independent
     * @param mass Wimp mass
     * @param sigma wimp-nucleon cross section
     * @param A    NATURAL, DEPLETED, or nuclear number
 */
  void plotSigmaSI(double mass=50.*GeV, double sigma=2.E-45, int A=NATURAL);

 /**
     * Plots the expected number of events for Xenon100(r) Run 10
     * and the theoretical exclusion limit in case of no candidates
     * @param runNumber    Run number to be mimicked
 */
  void plotLimitsSI(int runNumber);


 /**
     * Plots S1 and Er distributions for  SI, various LEff t-values
     * @param mass Wimp mass
     * @param runNumber    Run number to be mimicked
 */
  void plotS1SI(double mass=50.*GeV,int runNumber=10);

 /**
     * Plots S1 and Er distributions for  various models of SD
     * @param mass Wimp mass
     * @param runNumber    Run number to be mimicked
 */
  void plotS1SD(double mass=50.*GeV,int runNumber=10);


 /**
     * Plot Energy recoil for various SD parametrizations
     * @param mass Wimp mass
     * @param sigma wimp-nucleon cross section
 */
 void plotSigmaSD(double mass=50.*GeV,double sig=1.E-45*cm2);

 /**
     * Plot theoretical cross section limits for various SD parametrizations
     * @param runNumber    Run number to be mimicked
 */
  void plotLimitsSD(int runNumber=10) ;

  /**
    * Thorough comparison of FF and Structure functions for SD parametrizations
  */
   void plotFFSD();

/**
 * ML fit of the parameters of a distribution made of
 *    an exponential distribution plus a flat background
 */

   void testLK();

/**
 * shows how the Exclusion class runs on simple Signal+background PL
 * @param sMin Lower bound of the interval
 * @param sMax Upper bound of the interval
*/
   void testExclusionPL(int sMin=0,int sMax=10);
 
/**
 * Compute the coverage of Poisson confidence levels
 */
 void studyPoissonCoverage();

/**
 *  Print a table of combined PL
 */
 void tableCombinedPL();

/**
 *  Print a table of combined CLs thru Fisher method
 *  @param withCls Apply CLs correction?
 */
 void tableCombinedCLs(bool withCLs=true);

/**
 *  Print a table of simple CLs exclusoin
 *  @param withCls Apply CLs correction?
 */
 void tableCLs(bool withCLs=true);

/**
 * Compute the CLs limits making poisson a maximum likekihood
 *  @param withCls Apply CLs correction?
 * @param sigToE  Conversion from cm2 to number of events
 */
  void PLexclusionComputedBackground(bool withCLs=true, double sigToE=1.);

/**
 * Demonstrates Yellin's method
 * @param  s signal
 * @param  b background
 * @param n  number of points
 */
 void testYellin(int s=5, int b=0, int n=100);



/**
 * Plot cut acceptances, for both smeared and unsmeared cuts
*/
 void plotAcceptance();

/**
 * Technical method for interactive operation: setup an analysis
 * @param mode either PL_ANALYSIS or CUTS_ANALYSIS
 * @param doFitSystBkTvalue do we also fit the syst unc. on background?
 * @param doFitLeffTValue Do we want to fit LEFF t-value?
 * @param withS1Shape consider S1 shape?
 */

 void setupTheAnalysis(int mode=PL_ANALYSIS
                      , bool doFitSystBkTvalue=true
                      , bool doFitLeffTValue=true
                      , bool withS1Shape=true
                      );

/**
 * technical method for interactive operation: setup a run
 * @param publishedBackground use PRL background instead of computed one
 * @param mass Wimp Mass for initial print-out
 * @param runNumber Where to take data from
 */
 void setupTheRun(bool publishedBackground=true
                 , double mass=DEFAULT_WIMP_MASS,int runNumber=10);
 
/**
 * Print the result of the setup
 * @return  ok or not
 */
 bool printTheSetup();

/**
 * Show data, S1S2 bands for signal and background
 * @param mass Wimp Mass for initial print-out
 * @param runNumber Where to take data from
 * @param analysisMode PL_ANALYSIS or CUTS_ANALYSIS
 * @param displayMode S2_VS_S1, _OVER_S1_VS_S1, or FLATTENED_S2_VS_S1
 */
  void showRunData( bool publishedEr=true 
                  , double mass=DEFAULT_WIMP_MASS, int runNumber=10
                  , int analysisMode=PL_ANALYSIS
                  , int displayMode=FLATTENED_S2_VS_S1);

/** 
 * Compute limit of sigma SI for runs 8 and 10 Profile likelihood
 * @param s1 with s1 profile likelihood
 * @param publishedBkg use PRL background instead of computed one
 * @param doFitSystBkTvalue do we also fit the syst unc. on background?
 * @param specialMass if >0., print detailed computation for this mass
 * @param runNumber Where to take data from
 * @param limitSI compute the spin independent limit
 * @param limitSD compute the spin dependent limit
 */
  void runAnalysis(bool s1=true, bool publishedBkg=true
                  ,bool doFitSystBkTvalue=false
                  ,double specialMass=0.,int runNumber=10
                  ,bool limitSI=true,bool limitSD=false);

/**
 * Compute SI limit for combinations of runs
 * @param  nr number of runs (2->8,10  3->8,10,12)
 */
  void combineRuns(int nr=2);



/** 
 * Plot log likelihood as a function of number of events 
 * @param s1 with s1 profile likelihood
 * @param publishedBkg use PRL background instead of computed one
 * @param doFitSystBkTvalue do we also fit the syst unc. on background?
 * @param specialMass if >0., do it for one mass only
 * @aram nMasses number of masses to be drawns
 * @param runNumber Where to take data from
 */
  void plotLogLL(bool s1=true, bool publishedBkg=true
                ,bool doFitSystBkTvalue=false
                ,double specialMass=0.,int nMasses=0
                ,int runNumber=10);


/**
 * Compute Sensitivity bands of sigma 
 * @param s1 with s1 profile likelihood
 * @param publishedBkg use PRL background instead of computed one
 * @param nSimulations number of runs per mass points
 * @param runNumber Where to take data from
 */
 void computeSensitivityBands( bool s1=true, bool publishedBkg=true
                             , int nSimulations=100, int runNumber=10);
 
/**
 * Show how far we are from Wilk's theorem
 * @param number of masses to study (default: all of them)
 * @param nSimulations number of runs per mass points
 * @param s1 with s1 profile likelihood
 * @param publishedBkg use PRL background instead of computed one
 * @param runNumber Where to take data from
 */
void studyQValuesForBackground( int nm=0, int nSimulations=100000
                              , bool s1=true, bool publishedBkg=true
                              , int runNumber=10) ;

  Wimp             *wimp          ; /*!< a WIMP*/
  Target           *target        ; /*!< a mixture of Xenon isotopes */
  SIInteraction    *interSI       ; /*!< a Spin Independent interaction*/
  SDInteraction    *interSD       ; /*!< a Spin Dependent interaction*/
  XeRun            *run           ; /*!< a xenon run */
  OneDimSimulatedDataSet *simData ; /*!< simulated data set*/
  Exclusion        *exclusion     ; /*!< A machinery for exclusion finding*/
  XeDist           *dist          ; /*!<A distribution*/
  PValue           *pv            ; /*!< pvalue used in many examples*/
  S1S2PL           *pl            ; /*!< P.L. used in many examples*/
  S1S2Bands        *bands         ; /*!< a collection of S1/S2 bands*/
  XeSetOfCuts      *cuts          ; /*!< Set of cuts applied to the data  */
  XeGraph          *gr            ; /*!< An XeGraph used in many opprtunities*/
  XeMultiGraph     *mg            ; /*!< An XeMultigraph used a lot*/

  ClassDef(XeTest,1)
};
#endif

