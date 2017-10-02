#ifndef XeStat_h
#define XeStat_h
#include <fstream>
#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "XeMath.h"
#include "Math/DistFunc.h"
#include <TRandom3.h>
#include "XeCore.h"

enum content         { SPECTRUM
                     , VALUES
                     } ;

enum forceType       { DO_NOT_FORCE
                     , FORCE_TRUE
                     , FORCE_FALSE
                     } ;

enum analysisType    { NO_ANALYSIS
                     , PL_ANALYSIS
                     , CUTS_ANALYSIS
                     } ;

enum parameterType   { PARAMETER_OF_INTEREST
                     , NUISANCE_PARAMETER
                     , FIXED_PARAMETER
                     , FROZEN_PARAMETER
		     , FREE_PARAMETER
                     , N_PARAMETER_TYPES
                     } ;

enum sentivityModes  { MINUS_TWO_SIGMAS
                     , MINUS_ONE_SIGMA
                     , MEDIAN
                     , PLUS_ONE_SIGMA
                     , PLUS_TWO_SIGMAS
                     , N_SENSITIVITY_MODES
                     } ;

enum systematicModes { ONE_SIGMA_BELOW      =  0   // has to be zero
                     , CENTRAL
                     , ONE_SIGMA_ABOVE
                     , N_SIGMA_MODES        // should be 3
                     } ;

enum rejectionFactor { REJECT995
                     , REJECT9975
                     , REJECT999
                     , REJECT_ALL
                     , N_REJECTION_MODES
                     } ;

enum  basicParameter { PAR_SIGMA              = -1
		     , PAR_SYST_BKG_TVALUE    = -2
                     , PAR_LEFF_TVALUE        = -3 
                     , PAR_QY_TVALUE          = -4 
                     , PAR_EFFICIENCY_TVALUE  = -5
		     , PAR_SIGN_ACCEPTANCE_TVALUE = -6 
		     , PAR_NOT_ASSIGNED = -900 
		     , PAR_STAT_BKG_TVALUE    = -120 /*  run from -120 to -105*/
		     , PAR_GAUSS_TVALUE       = -220 /*  run from -220 to -120*/
                     } ;

enum   ciMode        { CI_UP
                     , CI_LOW	
                     , CI_TWO_SIDED
                     , CLS_UP
                     , CLS_LOW
                     , CLS_TWO_SIDED
                     , N_CI_MODES
                     } ;

enum sigmaModes   { ESTIMATED, UPPER_LIMIT };
enum sigmaUnits   { SIGMA_UNIT     , EVENT_UNIT };

static const  double DEFAULT_CL                    =            0.9 ;
static const  double DEFAULT_SIMULATED_X0          =            0.1 ;  
static const  double DEFAULT_SIMULATED_SIGMA       =            0.1 ;  
static const  double DEFAULT_SIMULATED_MU_GAUSS    =            0.5 ;  
static const  double DEFAULT_SIMULATED_XMAX_B      =            1.0 ;
static const  double LOGSQR2PI                     = 0.918938533205 ;
static const  double DEFAULT_EVENT_UPPER_LIMIT     =           20.0 ;



/**
   * Collection of static methods defining analysis methods.
   * Relevant flags: PL vs Cls, print Level
*/

class XeStat : virtual public XeMath , virtual public XeObject {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    
 /**
  * Set the analysis mode
  * @param mode must be either CUTS_ANALYSIS or PL_ANALYSIS
 */
    static void   setAnalysisMode(int mode); 

/**
  * Get the analysis mode
  * @return NO_ANALYSIS, CUTS_ANALYSIS or PL_ANALYSIS
 */
    static int    getAnalysisMode();

/**
 * Set print trace level
 * @param level print level
 */  
    static void   setPrintLevel(int level);

 /**
     * Check the anaylsis mode and set it if neccesary
     * @param name object name to be printed in case of error 
     * @param requested either NO_ANALYSIS, CUTS_ANALYSIS or PL_ANALYSIS
 */
    static bool   checkAnalysisMode(string name, int requested);
 

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
/**
 * Check wether the AnalysisMethod has been defined
 * @param verbose print a messae if not
 */
    static bool   isAnalysisDefined(bool verbose=false);
 
 /**
 * Is the current analysis cut based?
 */  
    static bool   isCutsBased();
 
 /**
 * Is the current analysis profile likelihood
 */  
    static bool   isPL();

/**
 * Get current print level
 */
    static int    getPrintLevel();

/**
 * Return defaut Lin/Log mode for sigma unit
 * @param unit SIGMA_UNIT or EVENT_UNIT
 */
    static int    getSigmaLinLog(int unit);

/**
 * Create default range if needed. Does not work in SIGMA_UNIT unit mode.
 * @param Original range, if null will create a new range;
 * @param unit SIGMA_UNIT or EVENT_UNIT
 */
    static XeRange*  newSigmaRange(XeRange*range=NULL, int unit=EVENT_UNIT);

/**
 * Return name of current analysis mode
 */
    static string getAnalysisModeName();

/**
 * Return name of  a given analysis mode
 * @param mode NO_ANALYSIS , PL_ANALYSIS or CUTS_ANALYSIS
 */
    static string getAnalysisModeName(int mode);

/** 
 * Return systematic error name
 * @param syst ONE_SIGMA_BELOW, CENTRAL, or ONE_SIGMA_ABOVE
 */
    static string getSystematicModeName(int syst);

/**
 * Return unit name for a given sigma unit
 * @param unit SIGMA_UNIT or EVENT_UNIT
 */
    static string getSigmaUnitName(int unit);

/**
 * Return sigma label for a given sigma unit
 * @param unit SIGMA_UNIT or EVENT_UNIT
 */
    static string getSigmaLabel(int unit);

/**
 * Return name a given mode (estimated or upper limit)
 * @param mode  ESTIMATED or UPPER_LIMIT 
 */
    static string getSigmaModeName(int mode);


/**
 * Return sigma limit label for a given sigma unit
 * @param unit SIGMA_UNIT or EVENT_UNIT
 */
    static string getUpperSigmaLabel(int unit);

    virtual ~XeStat();
    XeStat(string name);
    XeStat();

  protected:

    static int    printLevel; 
    static int    analysisMode;

} ;








 
/**
   * Description of a likelihood parameter. It is evaluated by the Likelihood class
*/
class LKParameter :  public XeStat {

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  public  :


    static string   getTypeName(int type);

    virtual        ~LKParameter();
   
/**
 * Empry constructor for ROOT
 */    
   LKParameter();

 /**
     * Regular Constructor 
     * @param id  Its identifier
     * @param type  Can be PARAMETER_OF_INTEREST 
     * , NUISANCE_PARAMETER , FIXED_PARAMETER , FROZEN_PARAMETER.  
     * FROZEN means temporary fixed
     * @param  nam  name
     * @param initial initial value
     * @param step step for minimization
     * @param  min  minimum value
     * @param  max  maximum value
 */
    LKParameter(int id, int type,string nam,double initial,double step
               ,double min, double max);

    void    printInitial();
    virtual void    printCurrent(bool withErrors);
    void    setId(int id);
    void    setType(int typ);
    void    setInitialValue(double i);
    void    setCurrentValue(double c);
    void    setCurrentValueInMinuitUnits(double c);


 /**
     * set the unit so that Minuit deals with parameters of order of mag 1
     * @param s  typical unit

 */
    void    setT0value(double val) {t0 = val;};
    double  getT0value() {return t0;};
    void    setMinuitUnit(double s=1.);
    void    setStep(double st);
    void    setSigma(double sig);
    void    setSigmaInMinuitUnits(double sig);
    void    setMinimum(double mi);
    void    setMaximum(double max);
    int     getType();
    int     getId();
    double  getCurrentValue();
    double  getCurrentValueInMinuitUnits();
    double  getInitialValue();
    double  getInitialValueInMinuitUnits();
    double  getMaximum();
    double  getMaximumInMinuitUnits();
    double  getMinimum();
    double  getMinimumInMinuitUnits();
    double  getStep();
    double  getStepInMinuitUnits();
    double  getSigma();
    string  getTypeName();

    double  getInitialSigma() { return initialSigma;};

    //return the gaussian constraint on the t-value
    double getLLGausConstraint();
    
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    

 /**
     * Set the fact that this parameter is or isn't common to various likelihoods calculators
     * @param c Is it common
 */
    void    setCommon(bool c=true);
    bool    isActive();
    bool    isOfInterest();
    bool    isCommon();

 /**
     *  freeze the parameter in profile likelihood double tag computation
     * @param doFreeze   Freeze is or not?
 */
    void    freeze(bool doFreeze);
                    
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    bool    compares(LKParameter* par, bool print);
    void    printHeader();
    bool    inCombinedMode();
    void    initialize();
    void    setCombinedMode(bool mode=true);

  protected :

    double  t0;
    int     type;
    int     id;
    bool    common;
    bool    combinedMode;     
    double  initialValue;
    double  step;
    double  minimum;
    double  maximum;
    double  currentValue;
    double  sigma;
    double  initialSigma;
    double  MinuitUnit;

     
} ;

/**
   * Predifined parameter: sigma
*/
class SigmaParameter : public LKParameter {
  public :
    SigmaParameter();
   ~SigmaParameter();

     
} ;

/**
   * Predefined parameter: t-value for systematic background uncertainty
*/
class TSystBkgParameter : public LKParameter {
  public :
    TSystBkgParameter(int run);
   ~TSystBkgParameter();

  static string getTheName(int run);

} ;

class TGaussParameter : public LKParameter {
  public :
    TGaussParameter(int id, int run);
   ~TGaussParameter();

    void   setUncertainty(double unc) { uncertainty = unc;};
    double getUncertainty()           { return uncertainty; };

    static string getTheName(int b, int run);

    double uncertainty;

} ;

/**
   * Predefined parameter: t-value for statistical background uncertainty
*/
class TStatBkgParameter : public LKParameter {

  public :
    TStatBkgParameter();
    TStatBkgParameter(int band, int run);
   ~TStatBkgParameter();


    void setStatError(double err);
    void printCurrent(bool withError);
  protected :

    static string getTheName(int b, int run);
    int band;
    double stat_error;
} ;


/**
   * Predefined parameter: t-LEff
*/
class TLEffParameter : public LKParameter {
  public :
    TLEffParameter();
   ~TLEffParameter();


} ;

/**
   * Predefined parameter: t-Qy
*/
class TQyParameter : public LKParameter {
  public :
    TQyParameter();
   ~TQyParameter();


} ;


/**
   * Predefined parameter: t-efficiency
*/
class TEfficiencyParameter : public LKParameter {
  public :
    TEfficiencyParameter();
   ~TEfficiencyParameter();


} ;


 /**
     * A likelihood object, consisting of parameters.
     * This is a virtual class
 */
class Likelihood :  virtual public XeStat {


   public:
    
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @return  Pure virtual method returning the log likelihood
     * @param   none  All parameters values set thru LKParameter
 */

     virtual  double computeTheLogLikelihood()=0;
     virtual ~Likelihood();

/**
 * Constructor
 * @param name name of the object
 */
     Likelihood(string name);

/**
 * Empty constructor for root
 */
     Likelihood();

 
/**
     * Method returning current estimated sigma
     * @return  current value of estimated sigma
 */
     double  getSigmaHat();   

 /**
 * Return the log likelilood for the estimated sigma
 */
     double  getLogD();

/**
 * Check wheher a parameter exists
 */
     bool     parameterExists(int p);

/**
 * Add a parameter
 * @param id either its hard coded id or the keyword AUTO
 * @param type NUISANCE_PARAMETER,PARAMETER_OF_INTEREST or FIXED_PARAMETER
 * @param name parameter name
 * @param initialVal  initial value
 * @param step        step for minimisation
 * @param mi          minimum value
 * @param ma          maximum value
 */
     void     addParameter(int id, int type,string nam,double initialVal
                         ,double step,double mi, double ma);

/**
 * add an existing parameter
 * @param  param pointer to existing LKParameter
 * @param  id    requested id, or SAME (the one coded in the parameter), or AUTO
 */
     void     addParameter(LKParameter* param,int id=SAME);

/**
 * Replace a parameter. If exists, delete the existing one
 * @param  param pointer to existing LKParameter
 */
   void     replaceParameter(LKParameter* param);

/**
 * Add a parameter with its own id, wheter already added or not
 * @param  param pointer to existing LKParameter
 */
    void     addParameterTolerant(LKParameter* param);

 /**
 * Activate (add or remove) a parameter with its own id
 * @param  param pointer to existing LKParameter
 * @param active  activate(true) or disable (true)
 */
     void     activateParameter(LKParameter* param, bool active=true);

 /**
 * Remove a parameter with its own id
 * @param  param pointer to existing LKParameter
 * @param  tolerant  don't issue a warning if it does not exist
 */
     void     removeParameter(LKParameter* param, bool tolerant=false);

 /**
 * Remove a parameter giving its id
 * @param  id  identifier
 * @param  tolerant  don't issue a warning if it does not exist
 */
     void     removeParameter(int id, bool tolerant=false);

     void     listParameters(); 
     void     printInitialParameters(); 
     void     printResultParameters(); 
     void     printCurrentParameters(); 
     void     ignoreParameter(int id);
     void     setParameterType(int id,int type);
     void     setParameterValue(int id,double v);
     void     setInitialValue(int id,double v);
     double   getParameterValue(int id);
     double   maximize(bool freezeParametersOfInterest);
     double   maximizeNumerically(int numberOfToys , bool freezeParametersOfInterest);
     void     setSeed(double Inputseed) {seed = Inputseed;};
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
     double   seed;               
     int      getNTotalParameters();
     int      getNActiveParameters();
     int      getNParametersOfInterest();
     int      getNMinuitParameters();
     int      getNNuisanceParameters();
     int      getNParameters(int type);
     void     resetParameters();
     LKParameter* getParameter(int id);
     map<int,LKParameter*>* getParameters();

     TH1F     getPullsHisto();  // return a TH1F with current values of all NP


 /**
     * return the shape LL that a set of values comes from a mixture of tabulated 
     * parent distribution. 
     * Makes sure that no bin has a negative expetancy!
     * @return Log Likelihood
     * @param bins list of bin number of input values
     * @param nDists number of distributions
     * @param  dists tabulated distributions 
     * @param  norm  normalized weights (i.e. normalized to 1.)
 */
     static double shapeLikelihood( vector<int>* bins, int nDist, double** dists
                                    , double *norm );

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     void     setCurrentValues(const double* v,const double* ers=NULL);
     void     setCurrentValuesInMinuitUnits(const double*v,const double*e=NULL);
     void     setCombinedMode(bool mode=true);
     void     printInitialHeader();
     void     printCurrentHeader();
     bool     checkConsistency();
     bool     inCombinedMode();
     int      mapMinuitParameters(bool freezeParametersOfInterest);
     int      getNParametersForChi2();
     void     forceNParametersOfInterest(int nF);
     void     clearTheParameters();

   protected:

/** 
 * Technical routine for constructors
 */
     void     setup();

     int                   currentId; 
     int                   nParametersOfInterest;
     int                   nActiveParameters;
     int                   nNuisanceParameters;
     int                   forcedNParametersOfInterest;
     bool                  combinedMode;     
     vector<LKParameter*>  MinuitParameters;
     map<int,LKParameter*> parameters;
    
     double       sigmaHat; /*!< Saved value of estimated sigma */

     double       LogD;   /*!< Saved value of Log Likelihood at estimated sigma*/
     
     void                  clear();
     bool                  checkParameter(int p, bool shouldExist);


} ;

typedef map<int,LKParameter*>::iterator ParameterIterator;

#define TRAVERSE_PARAMETERS(it) \
  for(ParameterIterator it=parameters.begin(); it!=parameters.end(); it++) 




/**
   *  Profile likelihood calculator
*/

class ProfileLikelihood : virtual public Likelihood {
    
  /* virtual class */ 

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
/**
 * Empty constructor for root
 */
      ProfileLikelihood();

/**
 * Constructor
 * @param name name of the object
 */
      ProfileLikelihood(string name);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
  
 
    virtual ~ProfileLikelihood();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/


 /**
     * print All flags and parameters
 */
    void    printFlagsAndParameters();


 /**
     *  estime cross section and save corresponding likelihood 
 */
     void    estimateCrossSection();
     bool  initialize() {return true;}; //FIXME ALE - place holder.
 

 /**
 * @return  an XeGraph of Likelihood for a sigma range
 * @param er Event range
 */
    XeGraph* newGraphOfLogLikelihood(EventRange* er=NULL);

/**
  * @return a TGraph of the likelihood for n_points sigma equal steps between min and max
  * @param n_points is the number of scan points for which the likelihood is maximized
*/
    TGraph* getGraphOfLogLikelihood(int n_points);
    //produce a graph of scan qValues scanning on mu, if forceMuhatComp=true then 
    //will recompute the minimum even if it was already computed and stored.
    TGraph* getGraphAtQval(double qval, bool forceMuhatComp=true);
    double getSigmaAtQval(double qval, bool forceMuhatComp=true);
    double returnLimitHagar( double cl, bool forceMuhatComp=true);

    //double getMaximum()
/**
  * @return a TGraph of parameter post fit values for different values of sigma
  * @param n_points is the number of sigma scan points 
  * @param param_index is the index value of the parameter
*/

    TGraph* getGraphOfParameter(int n_points,  int param_index);

    TGraph* getLikelihoodScanOfParameter( int n_points, LKParameter * par);

  virtual double getWimpMass();

  virtual void setData(int dataType)=0;

  virtual void generateAsimov(double mu_prime)=0 ;

  virtual void generateToyDataset(double seed, double mu_prime)=0;

  virtual double getSignalDefaultNorm()=0;
  virtual double getSignalMultiplier()=0; 
  virtual void   setSignalMultiplier(double val)=0; 


  protected : 

/** 
 * Technical routine for constructors
 */
    void    setup();
                 
    LKParameter *sigPar;  /*!< Pointer to main parameter of interest */


} ;

/**
   *  Combined Profile likelihood calculator for multi-experiments
*/

class CombinedProfileLikelihood : public ProfileLikelihood {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
           CombinedProfileLikelihood(string name);
          ~CombinedProfileLikelihood();
    void   combine(ProfileLikelihood* pl);
    double computeTheLogLikelihood();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    ProfileLikelihood* getProfile(int experiment);

 /**
     * print All flags and parameters
 */
    void               printFlagsAndParameters();
//    void               setWimpMass(double);
    double 	       getWimpMass();

    void setData(int dataType);

    void generateAsimov(double mu_prime) ;

    void generateToyDataset(double seed, double mu_prime);

    double getSignalDefaultNorm();
    double getSignalMultiplier();
    void   setSignalMultiplier(double val); 


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
             
   // bool    update();
    

  protected:

    bool    checkPValue();   
    map<int,ProfileLikelihood* > exps;
    int                          nCommon;
    double                       sigToEvents;


};

typedef map<int,ProfileLikelihood*>::iterator expIterator;

#define TRAVERSE_EXPERIMENTS(it) \
 for(expIterator it=exps.begin(); it!=exps.end(); it++) 


#endif



