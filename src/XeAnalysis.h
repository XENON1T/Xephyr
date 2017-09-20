#ifndef XeAnalysis_h
#define XeAnalysis_h

#include "TF1.h"
#include "XeCore.h"
#include "XeStat.h"
#include "XePhys.h"
#include <cmath>
#include <csignal>
#include <iostream>
#include "TGraph2D.h"


enum RunNumber { RUN_08 =  8
               , RUN_10 = 10
               , RUN_12 = 12
               } ;


static const int N_MAX_BANDS            = 20      // Maximum number of Bands
               , DEFAULT_N_BANDS_RUN_10 = 12
               , DEFAULT_N_BANDS_RUN_8  = 12
	       , DEFAULT_N_SLICES_RUN_8 = 27
               ;

const double  
        DEFAULT_MIN_S1_RUN_08           = 3.0*p_e // lower cut
      , DEFAULT_MIN_S1_RUN_10           = 3.0*p_e 
      , DEFAULT_MIN_S1_RUN_12           = 3.0*p_e 
      , DEFAULT_MAX_S1_RUN_08           = 30.*p_e // upper cut
      , DEFAULT_MAX_S1_RUN_10           = 30.*p_e
      , DEFAULT_MAX_S1_RUN_12           = 30.*p_e
      , DEFAULT_SIGMA_PMT               = 0.5*p_e // Fluctuations of response
      , DEFAULT_LIGHT_YIELD_RUN_08      = 2.2     // Pe/KeV at 122 KeV
      , DEFAULT_LIGHT_YIELD_RUN_10      = 2.28    // Pe/KeV at 122 KeV
      , DEFAULT_SE                      = 0.58    // Field quenching for e
      , DEFAULT_SR                      = 0.95    // Field quenching for nr
      , DEFAULT_US1_MIN                 = 1.*p_e  // Efficiency thresh. for uS1
      , DEFAULT_S2_OVER_S1_MIN          = 50./30.*unity
      //, DEFAULT_S2_OVER_S1_MIN          = 20.*unity
      , DEFAULT_S2_OVER_S1_MAX          = 4500./3.*unity
      //, DEFAULT_S2_OVER_S1_MAX          = 500.*unity
      , DEFAULT_FLATTENED_MAX           =  1.
      , DEFAULT_FLATTENED_MIN           = -DEFAULT_FLATTENED_MAX
      , DEFAULT_FLATTENED_FIT_MAX       =  0.45
      , DEFAULT_FLATTENED_FIT_MIN       = -DEFAULT_FLATTENED_FIT_MAX
      , DEFAULT_FLATTENED_ANOMALOUS_MAX =  0.5
      , DEFAULT_FLATTENED_ANOMALOUS_MIN = -DEFAULT_FLATTENED_ANOMALOUS_MAX
      , DEFAULT_LEFF_ER_MIN             =  3.0    // value under which LEff=0.
      ;

//  LEff Tabulation

const int N_LEFF_RUNS_8_AND_10 = 43;

enum   LEffTabulation { SIDE_LEFF_TABULATED = 25
                      , N_LEFF_TABULATED   =  2*SIDE_LEFF_TABULATED +1
                      } ;

static const  double LEFF_TVALUE_STEP=.1
                  , LEFF_TVALUE_MAX=SIDE_LEFF_TABULATED*LEFF_TVALUE_STEP
                  ;

static const  bool DEFAULT_EQUALLY_FILLED_S1_SLICES  = false
		, DEFAULT_DOMCSIGNALBAND             = false
                , DEFAULT_SPREAD_LEAKAGE_OVER_BANDS = true
                , DEFAULT_ER_GAUSSIAN_FLAT_IN_S1    = false
                ;

enum DATA_TYPE { 
     NO_DATA
   , REAL_DATA
   , SIMULATED_DATA  
   } ;


enum S1S2dataType { 
     UNSPECIFIED_DATA = - 1
   , AM_BE_DATA               /*!< Neutron calibration data*/
   , E_GAMMA_DATA             /*!< Gamma calibration data*/
   , DM_DATA                  /*!< Dark matter data*/

   , AM_BE_CUT_DATA           /*!< Neutron calibration data after cuts*/
   , E_GAMMA_CUT_DATA         /*!< Gamma calibration data after cuts*/
   , DM_CUT_DATA              /*!< Dark matter data after cuts*/

   , ER_BACKGROUND            /*!< Gaussian ER background*/
   , ER_LEAKAGE               /*!< Anomalous leakage*/
   , NR_BACKGROUND            /*!< Neutron background*/
   , ALL_BACKGROUNDS          /*!< All backgrounds*/
   , DM_SIMULATED_DATA        /*!< Dark matter simulated data, generated from BKG only, 
				you have to fill this calling XeRun::generateData()*/

   , DM_CUT_DATA_INTEGRATED
   , ASIMOV_DATA	     /*!, Asimov dataset generated with the function XeRun::generateAsimovData(mu)*/	
   , ALL_BACKGROUNDS_INTEGRATED
   , DEFAULT_SIGNAL   
   , DEFAULT_SIGNAL_INTEGRATED
   , FIRST_BACKGROUND = ER_BACKGROUND
   , LAST_BACKGROUND  = ALL_BACKGROUNDS
   , N_BACKGROUNDS    = 1 + ALL_BACKGROUNDS
   , N_DATA_TYPES     = 1 + DM_DATA
   , N_CUT_DATA_TYPES = 1 + DM_CUT_DATA
   , N_ALL_DATA_TYPES = 1 + DEFAULT_SIGNAL_INTEGRATED
   } ;

enum S1S2shape { SHAPE_SIGNAL
               , SHAPE_ER_BACKGROUND
               , SHAPE_ER_LEAKAGE
               , SHAPE_NR_BACKGROUND
               , N_SHAPES
               , SHAPE_DATA = N_SHAPES
               } ;


enum BACKGROUND_MODEL_TYPE {
     ELECTRON_BACKGROUND_MODEL
   , NR_BACKGROUND_MODEL
   };
   

enum S1S2_COLUMNS {
     S1_COLUMN
   , S2_COLUMN
   , MAX_NUMBER_OF_EVENTS = 5000000 //Ale it was 50000  and Ambe was not applying cuts!!!!!!
   } ;

// For cuts

static const  double        S2_MIN_RUN_8  = 300. ;
static const  double        S2_MIN_RUN_10 = 150. ;

enum  cutMode            { CUT_UNKNOWN
                         , CUT_BELOW
                         , CUT_ABOVE
                         } ;

enum  selectionCutMode   { SELECTION_CUT_UNKNOWN
                         , SELECTION_CUT_ON_UNSMEARED_S1
                         , SELECTION_CUT_ON_SMEARED_S1
                         } ;

enum  selectionCutNature { USER_GENERAL_CUT
                         , USER_SELECTION_CUT
                         , USER_SELECTION_CUT_ON_S1
                         , USER_SELECTION_CUT_ON_UNSMEARED_S1
                         , S1_STOT_SELECTION_CUT
                         , S2_PEAK_S0_SELECTION_CUT
                         , S1_COIN2_SELECTION_CUT
                         , DARK_MATTER_CUT
                         } ;

class XeRun;
class SignalModel;
class MCSignalModel;
class ModelImporter;

/**
 * Run component
 * Describes objects associated to an XeRun
 */
class RunComponent : virtual public XeObject {

  public:

/**
 * Void constructor for ROOT
 */
   RunComponent();
   
/**
 * Regular constructor
 * @param name name of the object
 */
   RunComponent(string name);

/** 
 * Virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  virtual bool isRunCompatible(int runNumber);

/**
 * Do check run compatibility
 * @param warn Do print a warning if not compatible
 */
  bool checkRunCompatibility(bool warn=true);


/**
 * Destructor
 */
   virtual ~RunComponent();

/**
 * Return run number
 */
     int      getRunNumber();

/** 
 * Return the run it applies to
 * @return a pointer to the XeRun object
 */
     XeRun*   getRun();


/**
 * Set the Run
 * @param run Pointer to XeRun
 */
   void setRun(XeRun *run);

 /**
     * print flags for tracing purpose.
     * actually print the flags of the parent XeRun
 */
    void traceTheFlags();


   protected: 

    XeRun* run;
    int    runNumber;

   ClassDef(RunComponent,1)

};


/**
 * Detector component (for documentation purpose only)
 */
class DetectorComponent : virtual public XeMath, virtual public XeObject {

  public:

/**
 * Void constructor for ROOT
 */
   DetectorComponent();
   
/**
 * Regular constructor
 * @param name name of the object
 */
   DetectorComponent(string name);

/**
 * Destructor
 */
   virtual ~DetectorComponent();

   ClassDef(DetectorComponent,1)
};

class Qy :  public DetectorComponent {

    public:	
	Qy();

/**
   * return minimum t-value for LEff
*/
   double getMinimumTval(){return -maxTval;};	

/**
   * return maximum t-value for LEff
*/
   double getMaximumTval(){return maxTval;};	

/**
   * set minimum t-value for LEff
*/
   void  setMaximumTval(double max) {maxTval=max;};	

/**
   * set number of steps from [0., max] from [min, 0.] have the same configuration,
   this assumes that min = -max. (you cannot actually do differently).
*/
   void  setNstepsHalfway(double Nsteps) {TvalStepSize= maxTval /Nsteps; };	

   void  setStepSize(double steps) {TvalStepSize= steps; };	

/**
   * get step size 
*/
   double  getStepSize() {return TvalStepSize;};	


  protected :
 
    double     tValue;  
    double     maxTval;
    double     TvalStepSize;

	

};


/**
   * Description of LEffective.
   * This is a virtual class
*/
class LEff : public DetectorComponent {

  public:

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
/**
     * base class creator
*/
             LEff();

 /**
 * Compute graph of LEFF
 * @return a graph (actually XeGraph) of LEffective 
 * @param  er Energy Range (actually ErRange)
 * @param  t  Nuisance parameter modified LEff, in units of "sigma equivalent"
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph *newGraphOfLEff(ErRange *er=NULL, double t=0., int plot=NONE);

  /**
 * Compute Multi graph of LEFF for canonical values of LEff t-values
 * @return a Multi-graph (actually XeMultiGraph) of LEffective 
 * @param  er Energy Range (actually ErRange)
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
   XeMultiGraph *newMultiGraphOfLEff(ErRange *er=NULL,int step=1,int plot=NONE);  

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
 /**
     * compute LEff t-value for a given tabulation index
     * @return The default LEff class. Returns NULL if unknwoen run number
     * @param runNumber Run number
 */
    static LEff* newDefault(int runNumber);

 /**
     * compute LEff t-value for a given tabulation index
     * @param T_INDEX index
 */
    static double tabulated(int T_INDEX);

 /**
     * return name of LEff t-value for a given tabulation index
     * @param T_INDEX index
 */
    static string tabulatedName(int T_INDEX);

/**
   * return minimum t-value for LEff
*/
   double getMinimumTval(){return -maxTval;};	

/**
   * return maximum t-value for LEff
*/
   double getMaximumTval(){return maxTval;};	

/**
   * set minimum t-value for LEff
*/
   void  setMaximumTval(double max) {maxTval=max;};	

/**
   * set number of steps from [0., max] from [min, 0.] have the same configuration,
   this assumes that min = -max. (you cannot actually do differently).
*/
   void  setNstepsHalfway(double Nsteps) {TvalStepSize= maxTval /Nsteps; };	

   void  setStepSize(double steps) {TvalStepSize= steps; };	

/**
   * get step size 
*/
   double  getStepSize() {return TvalStepSize;};	



 /**
     * Set the interpolate mode wrt t-value
     * @param  mode  can be either LINEAR or EXPONENTIAL
 */
    static void              setTValueMode(int mode);
    static LinearRange*      getTRange();
    static pair<int,double>  getBinAndFraction(double t);


    void          setTValue(double t);
    void          setErMin(double e);
    double        getTValue();
    double        getErMin();
    double        getLEff(double Er);
    double        getLEff(double Er, double t);
    XeGraph      *newGraphOfCurrentLEff(ErRange *er=NULL,int plot=NONE);  
    virtual      ~LEff();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    virtual  double getLEffFromTable(double Er,int mode)=0;

  protected :
 
    static LinearRange *tRange;
    static int          tValueMode;

    double     tValue;  
    double     ErMin;  
    double     maxTval;
    double     TvalStepSize;

    ClassDef(LEff,1)

};

/**
   * Implementation of LEffective with latest 2011 measurements.
   * This is a real class
*/
class LEff2011 : public LEff {

  public:


    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    LEff2011();
   ~LEff2011();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    double getLEffFromTable(double Er,int mode=CENTRAL);
    static double  _leff[N_ER_POINTS][N_SIGMA_MODES];

    ClassDef(LEff2011,1)

};

/**
   * Implementation of LEffective as coded in runs 8 and 10
   * This is a real class
*/
class LEffRuns8And10 : public LEff {

  public:


    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    LEffRuns8And10();
   ~LEffRuns8And10();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    double getLEffFromTable(double Er,int mode=CENTRAL);
    static double  _leff[N_SIGMA_MODES][N_LEFF_RUNS_8_AND_10];
    static double  _er[N_LEFF_RUNS_8_AND_10];

    ClassDef(LEffRuns8And10,1)

};


class RunFlattener : virtual public Flattener, virtual public RunComponent {

  public:

    ~RunFlattener();

/**
 * Zero constructor for TObject
 */
     RunFlattener();

/**
 * Constructor for a given run
 * @param run Pointer to the XeRun
 */
     RunFlattener(XeRun* run);

/**
 * Compute the flattened log(s2/s1)
 * @return flattened log s2/s1
 * @param   S1 original S1
 * @param   S2 original S2
 *                 */
    double   flatten(double S1, double S2);

 /**
 * Return S1 given S1 and flattened s2/s1
 * @return S2 from flattened log s2/s1
 * @param   S1 original S1
 * @param   flat flattened quantity
 */
    double   unflatten(double S1, double flat);

/**
 * Internal method to establsh the name
 */
    static string getTheName(XeRun* run);

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);

  protected: 

    double eBandMin;
    double eBandMax;
    TF1*   eBandFlat;


    ClassDef(RunFlattener,1)
};

/**
   * Description of a cut, either a preselection or a selection one.
   * Cuts are applied to real data and return "pass or fail".
   * They also return acceptance on S1 only.
   * This is a virtual class
*/

class XeCut :  virtual public S1S2Object, virtual public RunComponent {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
 * Regular constructor
 * @param  run the run (actually XeRun) on which the cut applies
 * @param  name name of the cut
 */
    XeCut(XeRun* run, string name);

/**
 * Empty constructor for root
 */
     XeCut();

 /**
     * Do S1 and S2 pass the cut? (virtual)
     * @return a boolean, true if passes
     * @param S1  measured S1
     * @param S2  measured S2
 */
     virtual  bool     passIt(double S1,double S2);

/**
 * produce a graph and optionally plots it
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     virtual  XeGraph* newGraph(int plot=NONE);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     virtual ~XeCut();
 /**
     * return  cut nature
     * @return nature: 
     * USER_GENERAL_CUT , USER_SELECTION_CUT , USER_SELECTION_CUT_ON_S1
     *  USER_SELECTION_CUT_ON_UNSMEARED_S1 , S1_STOT_SELECTION_CUT
     *  S2_PEAK_S0_SELECTION_CUT , S1_COIN2_SELECTION_CUT , DARK_MATTER_CUT
 */
     int      getNature();

/**
 * Draw the cut, being and S1S2Object object
 */
     void     draw();

/**
 * Enable the cut
 * @param enabled it is to be enabled?
 */
     void     enable(bool enabled);

/**
 * Is he cut enabled?
 * @return true or false
 */
     bool     isEnabled();

 /**
     * Do S1 and S2 pass the cut? If not enabled, return true
     * @return a boolean, true if passes
     * @param S1  measured S1
     * @param S2  measured S2
 */
     bool     passes(double S1,double S2);

/**
 * Return the funtion used to define the cut, if any
 * @return pointer to the cut definition function, or NULL if not such a thing
 */
     TF1*     getCutFunction();

 /**
     * return  the Overall Acceptance.
     * @return acceptance for Dark Matter cuts, 1. for preselection cuts;
*/
     double   getOverallAcceptance();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
        
/**
 * Return maximum 'X' value, as requested by being an S1S2Object object
 */
    double maxX();

/**
 * Return maximum 'Y' value, as requested by being an S1S2Object object
 */
    double maxY();

/**
 * Return minimum 'X' value, as requested by being an S1S2Object object
 */
    double minX();

/**
 * Return minimum 'Y' value, as requested by being an S1S2Object object
 */
    double minY();


              
    protected :  
       
     bool     enabled;
     int      cutMode;
     int      nature;
     TF1     *cutFunction;
     double   overallAcceptance;
    
     static constexpr   int    N_PE   = 40 ;
     static constexpr   double PE_MIN =  1.;
     static constexpr   double PE_MAX = 40.;

      ClassDef(XeCut,1)
};

/**
   * Set of cuts
*/

class XeSetOfCuts : public S1S2Object, public vector<XeCut*>
                  , public RunComponent  {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     XeSetOfCuts();
 /**
 * create a mutiple graph (actually XeMultiGraph) of all cuts in S1 or uS1
 * @return pointer to the newly created XeMultiGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     XeMultiGraph*   newMultiGraph(int plot=NONE);

/**
     * create the Overall Acceptance (1. for preselection cuts);
     * @return pointer to the newly created XeMultiGraph
*/
     double   getOverallAcceptance();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

 /** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);

                   
                    ~XeSetOfCuts();
/**  
 * Draw the cuts in a standard S1,S2 plane
 */
     void            draw();
  
/**
 * Return maximum 'X' value, as requested by being an S1S2Object object
 */
    double        maxX();

/**
 * Return maximum 'Y' value, as requested by being an S1S2Object object
 */
    double        maxY();

/**
 * Return minimum 'X' value, as requested by being an S1S2Object object
 */
    double        minX();

/**
 * Return minimum 'Y' value, as requested by being an S1S2Object object
 */
    double        minY();



/**  
 *  remove all the cuts
 */
     void            reset();
  
 
/**  
 *  adds one cut to the set 
 *  @param cut the XeCut to be added
 */
     void            addCut(XeCut* cut);
  
/**  
 *  adds one set of cutc to the current set 
 *  @param cuts the XeSetOfCut to be added
 */
     void            addCuts(XeSetOfCuts* cuts);

/** 
 * enable one cut
 * @param i  cut index, counting from 0
 * @param enabled enable or disable
 */
     void            enable(int i, bool enabled=true);

/** 
 * disable one cut
 * @param i  cut index, counting from 0
 */
     void            disable(int i);

/** 
 * enable one cut
 * @param enabled enable or disable
 */
     void            enableAll(bool enabled=true);

/** 
 * disable all the cuts
 */
     void            disableAll();

/**
 * Do s1 and s2 pass the cut?
 *  @param s1 s1 value
 *  @param s2 s2 value
 */
     bool            passes(double s1,double s2);

/**
 * enable first n cuts, disable all the rest
 * @param n how many need to be enabled
 */
     void            enableFirst(int n);

/**
 * is a cut enabled?
 * @param i  cut index, counting from 0
 */
     bool            isEnabled(int i);
/**
 * return number of cuts
 */
     int             getNCuts();

/**
 * Return name of a cut
 * @param i  cut index, counting from 0
 */
     string          getCutNameBySequence(int i);

/**
 * Get a cut by its sequence number
 * @return pointer to the requested cut
 * @param i  cut index, counting from 0
 */
     XeCut*          getCutBySequence(int i);

/**
 * get first cut of a given nature
 * @return pointer to first cut of the given nature
 * @param nat cut nature :
 * USER_GENERAL_CUT , USER_SELECTION_CUT , USER_SELECTION_CUT_ON_S1
 *  USER_SELECTION_CUT_ON_UNSMEARED_S1 , S1_STOT_SELECTION_CUT
 *  S2_PEAK_S0_SELECTION_CUT , S1_COIN2_SELECTION_CUT , DARK_MATTER_CUT
 */
     XeCut*          getCutByNature(int nat);

 /**
     * create a default set of dark Matter cuts
     * @return a pointer to the newly created set of dark Matter cuts
     * @param  run  The given run
 */
     static XeSetOfCuts* newDefaultXeSetOfDarkMatterCuts(XeRun* run);


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
/**  
 *  print the set of cuts
 *  @param level print level
 */
     bool            printIt(int level=1);
                    
   protected:
 
     virtual void    setTheName();
     void            init();
     bool            checkIndex(int i);

      ClassDef(XeSetOfCuts,1)
};

/**
   * Selection cut.
   * Also a virtual class
*/

class SelectionCut : public XeCut{


  public:

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    virtual ~SelectionCut();
             SelectionCut();
 /**
     * @param  run XeRun on which this applies
     * @param  mode either SELECTION_CUT_ON_UNSMEARED_S1
     *              or SELECTION_CUT_ON_SMEARED_S1
     * @param  variable whether the acceptance can be modified
     *          by a nuisance parameter          
 */

     SelectionCut(XeRun* run,int mode, bool variable);
 /**
 * create  a graph of the acceptance
 * @return  a pointer to the newly created  XeGraph
 * @param  pe Range in photo electron (actually S1Range, NULL=default one)
 * @param  t  Nuisance parameter modifying the acceptance, in sigma equivalent
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph*  newGraphOfAcceptance(S1Range* pe=NULL, double t=0.,int plot=NONE);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                   
/**  
 *  set "t-value" (how many sigmas are we above or below the average?)
 *  @param level t t-value
 */
    void             setTValue(double t);

/**
 * Are we in smeared or unsmeared mode?
 * @return SELECTION_CUT_UNKNOWN , SELECTION_CUT_ON_UNSMEARED_S1 , SELECTION_CUT_ON_SMEARED_S1
 */
    int              getSmearMode();
/** 
 * Is this cut driven by a "t-value" ?
*/
    bool             isVariable();

/**
 * Return the acceptance as a function of (un)smeared S1 for the central value
 * @param S1  (un)smeared s1
 */
    double           getAcceptance(double S1);

/**
 * Return the acceptance as a function of (un)smeared S1
 * @param S1  (un)smeared s1
 * @param t "t-value"
 */
    double    getAcceptance(double S1,double t);

/** 
 * return the currently set "t-value"
 */
    double    getTValue();
/**
 * get the name of the smear mode ("smeared", "unsemared")
 */ 
    string    getSmearModeName();


 /**
 * create  a graph of the acceptance for the current t-value
 * @return  a pointer to the newly created  XeGraph
 * @param  pe Range in photo electron (actually S1Range, NULL=default one)
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph*  newGraphOfCurrentAcceptance(S1Range* pe=NULL,int plot=NONE);   

/**
 * Set the effect of the "t-value"
 * @param mode either LINEAR or LOG
 */
    static void      setTValueMode(int mode);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    void             tabulate();
                    
    
/**  
 *  print the selection cut
 *  @param level print level
 */
    bool             printIt(int level=1);
   protected: 
 
    static int       tValueMode;
    int              smearMode;
    double           tValue;
    bool             variable;
    bool             tabulated;

    double          _acc[N_ER_POINTS][N_SIGMA_MODES];
  
    double           getAcceptance(double S1,int mode);
    virtual double   computeAcceptance(double S1,int mode)=0;

      ClassDef(SelectionCut,1)
 
};

/**
   * Selection cut expressed as a function of S1.
   * A virtual class
*/

class SelectionCutS1 : public SelectionCut {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public:
    virtual ~SelectionCutS1();
             SelectionCutS1();
/**
 * Constructor of this virtual class
 * @param run current XeRun
 * @param variable The acceptance can be described by a "t-value"
 */
             SelectionCutS1(XeRun* run, bool variable);

 
      ClassDef(SelectionCutS1,1)
};

/**
   * Selection cut expressed as a function of unsmeared S1.
   * A virtual class
*/

class SelectionCutUnsmearedS1 : public SelectionCut {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public:
    virtual ~SelectionCutUnsmearedS1();
             SelectionCutUnsmearedS1();
 /**
 * Constructor of this virtual class
 * @param run current XeRun
 * @param variable The acceptance can be described by a "t-value"
 */
             SelectionCutUnsmearedS1(XeRun* run, bool variable);

      ClassDef(SelectionCutUnsmearedS1,1)

};

/**
   * the 's1sTotCut' preselection cut
*/

class s1sTotCut : public SelectionCutS1 {

  public:

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
/**
 * Regular constructor
 * @param run the XeRun it applies to
 */
    s1sTotCut(XeRun* run);
    s1sTotCut();

/**
 * Create a graph to plot it un the current S1S2 representation
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraph(int plot=NONE);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
/**
 * get default S1Min
 * @param run run number
 */
    static   double defaultS1Min(int run);
                    
/**
 * get default S1Max
 * @param run run number
 */
    static   double defaultS1Max(int run);

            ~s1sTotCut();
                    
/**
 * get S1Min
 * @param minS1  Minimum S1
 */
    void     setMinS1(double minS1);
                    
/**
 * get S1Max
 * @param maxS1  Maximum S1
 */
    void     setMaxS1(double maxS1);

/**
 * return current minimum S1
 */
    double   getMinS1();

/**
 * return current maximum S1
 */
    double   getMaxS1();

 /** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);


/**
 * Do s1 and s2 pass this cut?
 */
    bool     passIt(double s1,double s2);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

    double   computeAcceptance(double S1,int mode);
          
    void     setTheName();
  protected:

    double   MinS1;
    double   MaxS1;
    LEff    *leff;


      ClassDef(s1sTotCut,1)

};

/**
   * the 's2peaks0Cut' preselection cut
*/

class s2peaks0Cut : public SelectionCutUnsmearedS1 {
  
  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
             s2peaks0Cut();
                    
/**
 * Regular constructor
 * @param run the XeRun it applies to
 */
             s2peaks0Cut(XeRun* run);

/**
 * Create a graph to plot it un the current S1S2 representation
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraph(int plot=NONE);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
            ~s2peaks0Cut();
 
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

 /** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);

          
  protected :

    double             minS2;
    static   double   _acc_run_08[N_PE_POINTS];

/**
 * Get default s2Min
 * @param run number
 */
    static   double s2Min(int run);

/*
 * Compute the acceptance
 * @param uS1 unsmeared S1
 * @param mode must be CENTRAL
 */
    double   computeAcceptance(double uS1,int mode);
    void     setTheName();

      ClassDef(s2peaks0Cut,1)

};

/**
   * the 'S1coin2Cut' preselection cut
*/

class S1coin2Cut : public SelectionCutS1 {
  
  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
            S1coin2Cut();
                    
/**
 * Regular constructor
 * @param run the XeRun it applies to
 */
            S1coin2Cut(XeRun* run);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
           ~S1coin2Cut();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  protected :

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);


/*
 * Compute the acceptance
 * @param S1 smeared S1
 * @param mode must be CENTRAL
 */
    double  computeAcceptance(double S1, int mode);

      ClassDef(S1coin2Cut,1)

};

/**
   * This is for run 10, and excludes s1coin2
*/

class OtherSelectionCutsS1 : public SelectionCutS1 {

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  public :

            OtherSelectionCutsS1();
/**
 * Regular constructor
 * @param run the XeRun it applies to
 */
            OtherSelectionCutsS1(XeRun* run);
           ~OtherSelectionCutsS1();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  
  protected :


/*
 * Compute the acceptance
 * @param S1 smeared S1
 * @param mode ONE_SIGMA_BELOW, CENTRAL, or ONE_SIGMA_ABOVE
 */
    double  computeAcceptance(double S1,int mode);

      ClassDef(OtherSelectionCutsS1,1)

} ;

/**
   * This is for run 8, and includes s1coin2
*/

class AllSelectionCutsS1 : public SelectionCutS1 {
  
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  public :

/**
 * Regular constructor
 * @param run the XeRun it applies to
 */
            AllSelectionCutsS1(XeRun* run);
            AllSelectionCutsS1();
           ~AllSelectionCutsS1();
          
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    static double   _acc_run_08[N_PE_POINTS];

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);

  protected :

/*
 * Compute the acceptance
 * @param S1 smeared S1
 * @param mode must be CENTRAL
 */
    double  computeAcceptance(double S1,int mode);

      ClassDef(AllSelectionCutsS1,1)

} ;


/**
   * A set of selection cuts
*/

class XeSetOfSelectionCuts : public XeSetOfCuts  {

   public :
        
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     XeSetOfSelectionCuts();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     virtual      ~XeSetOfSelectionCuts();

/**
 * Add a selection cut
 * @param cut pointer to the selection cut to be added
 */
     void          addSelectionCut(SelectionCut* cut);

/**
 * Cumulate an already existing set of selection cuts
 * @param sc pointer to the set of selection cuts to be cumulated
 */
     void          addSelectionCuts(XeSetOfSelectionCuts* sc);
                    
                   
/**  
 *  print the remaining count of cuts onr by one, in a set of selection cuts
 *  @param level print level
 *  @param remaing list of remaining counts (0=before, i=after cut "i-1"
 */
     bool          printIt(int level,vector<double>* remaining);

/**
 * set the "t-value" for a given cut
 * @param i  cut sequence number, starting at 0
 * @param t "t-value"
 */
     void          setTValue(int i, double t);

/**
 * get the overall acceptance for all cuts on either S1 or unsmeared S1
 * @param S1 either unsmeared or smeared S1, depensing on parameter 
 * @param smeared  are we dealing with smeared or unsmeared S1?
 */
     double        getAcceptance(double S1,bool smeared);

/**
 * get overall max cut-value on S1
 */
     double        getMaxS1();

/**
 * get overall min cut-value on S1
 */
     double        getMinS1();

/**
 * get selection cut by sequence
 * @return pointer to the desired selection cut
 * @param i  cut sequence number, starting at 0
 */
     SelectionCut* getSelectionCutBySequence(int i);
/**
 * get selection cut by nature
 * @return pointer to the desired selection cut
 * @param nature nature of the cnut. Can be:
 * USER_GENERAL_CUT , USER_SELECTION_CUT , USER_SELECTION_CUT_ON_S1
 *  USER_SELECTION_CUT_ON_UNSMEARED_S1 , S1_STOT_SELECTION_CUT
 *  S2_PEAK_S0_SELECTION_CUT , S1_COIN2_SELECTION_CUT , DARK_MATTER_CUT
*/
     SelectionCut* getSelectionCutByNature(int nature);

 /**
     * return a default set of selection cuts
     * @param  run  The given run
 */
     static XeSetOfSelectionCuts* newDefault(XeRun* run);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     void    setTheName();


  /**  
 *  print the Set of selection cuts
 *  @param level print level
 */
     bool          printIt(int level=1);
   protected:
 
     int            mode;

      ClassDef(XeSetOfSelectionCuts,1)
} ;


/**
   *  Cut for cut-based analysis.
   *  A virtual class
*/

class DarkMatterCut : public XeCut {

  /* virtual class */

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @param  run run (actually XeRun) on which the cut apply
     * @param  name name of the cut
     * @param  efficiency background rejection
 */

    DarkMatterCut(XeRun* run, string name,double efficiency);
 /**
     * @param  run run (actually XeRun) on which the cut apply
     * @param  name name of the cut
     * @param  mode either REJECT995, REJECT9975, REJECT999, REJECT_ALL
 */
    DarkMatterCut(XeRun* run, string name,int mode);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

            DarkMatterCut();
    virtual ~DarkMatterCut();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  protected :

    double   efficiency;
    static   string getTheName(string name,int mode);
    static   string getTheName(string name,double efficiency);

      ClassDef(DarkMatterCut,1)

};

/**
   * Dark Matter cut based on S1S2.
   * A virtual class
*/

class S1S2SingleCut : public DarkMatterCut {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/

  /**
     * @param  run run (actually XeRun) on which the cut apply
     * @param  name name of the cut
     * @param  above  is the cut above or below
     * @param  efficiency background rejection
 */

                 
             S1S2SingleCut(XeRun* run,string name,bool above,double efficiency);
   /**
     * @param  run run (actually XeRun) on which the cut apply
     * @param  name name of the cut
     * @param  above  is the cut above or below
     * @param  mode either REJECT995, REJECT9975, REJECT999, REJECT_ALL
 */
    S1S2SingleCut(XeRun* run, string name,bool above, int mode);

    S1S2SingleCut();

/**
 * Create a graph to plot it un the current S1S2 representation
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraph(int plot=NONE);

 /**
     * return Do S1 and S2 pass the cut? (virtual)
     * @return a boolean, true if passes
     * @param S1  measured S1
     * @param S2  measured S2
 */
    bool     passIt(double S1,double S2);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    virtual ~S1S2SingleCut();

/**
 * Describing the S2 cut value
 * @return the value according to which S2 is cuy
 * @return S1 S1-value
 */
    virtual  double s2Cut(double S1)=0;

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

                  
   protected :

    bool       above;
/**
 * Define as a cut above or below somthing
 * @param above true is accept stuff above
 */
    void       setMode(bool above);    
    int        spaceMode;
    Flattener* flattener;

      ClassDef(S1S2SingleCut,1)

};


/**
   * Real Implementation of S1/S2 Dark Matter cut
*/

class S1OverS2Cut : public S1S2SingleCut {

  public  :
  
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/

                   
    S1OverS2Cut(); 
    /**
     * @param  run run (actually XeRun) on which the cut apply
     * @param  mode either REJECT995, REJECT9975, REJECT999, REJECT_ALL
    */
    S1OverS2Cut(XeRun* run, int mode); 

    /**
     * @param  run run (actually XeRun) on which the cut apply
     * @param  efficiency  background rejection 
    */
                   
    S1OverS2Cut(XeRun* run, double efficiency);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
   ~S1OverS2Cut();

/**
 * Return the s2 cutting value
 * @param s1 S1 value
 */

    double   s2Cut(double s1);
 /**
     * return Do S1 and S2 pass the cut? (virtual)
     * @return a boolean, true if passes
     * @param S1  measured S1
     * @param S2  measured S2
 */
    bool     passIt(double s1, double s2);

/** 
 * define the ROOT funtion defining the cut
 * @param mode  Either REJECT995, REJECT9975, or REJECT999
 */
    void     setCutFunction(int mode);    

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);



    double getTheOverallAcceptance();
    void   setRunAndMode(XeRun* r, int mode);
   
  protected :

    int      rejectionMode;
    static   string getOfficialName(int run); 

      ClassDef(S1OverS2Cut,1)

};


class S1S2Bands;
 
/**
   * DataSet in S1 S2
*/

class S1S2Data : virtual public S1S2Object, virtual public DataSet
               , virtual public RunComponent{

  public:

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    

    S1S2Data();
 /**
     * Constructor ab nihilo 
     * @param   name  name of data set
     * @param   run   run (actually XeRun) where data come from
     * @param   gr    Graph where S1S2 are taken from 
 */
    S1S2Data(string name,int type,XeRun* run,TGraph* gr=NULL);

 /**
     * Constructor from an existing data set, after further cuts
     * @param original  original S1S2Data before cuts
     * @param cuts      set of preselection cuts
 */

    S1S2Data(S1S2Data* original, XeSetOfCuts* cuts);

 /**
     * Constructor from an existing data set, after dark matter cut
     * @param original  original S1S2Data before cuts
     * @param cut       dark matter cut
 */
    S1S2Data(S1S2Data* original, DarkMatterCut* cut);

  /**
     * fill bands and slices and mark events with their band and slice number
     * @param   bands S1S2Bands to be filled from which the events are read
     * @param   mark     mark the band and slice numbers for each event
 */
    void fillAndMarkBandsAndSlices(S1S2Bands* bands,bool mark);
 
    virtual ~S1S2Data();


 /**
     * Multiply everything by a constant 
     * @param  ratio multiplication factor
 */
    void normalize(double ratio);

 /**
     * Multiply everything by a constant factor such as the total content will
     * be as requested
     * @param   events requested total content
 */
    void normalizeToEvents(double events);

/**
 * Add an entry corresponding to s1,s2 measurements
 * @param s1  S1 value
 * @param s2  S2 value
*/
    void add(double s1, double s2);

/**
 * Add an entry corresponding to s1,s2 measurements, and record band/slice indices
 * @param s1  S1 value
 * @param slice  slice index (couting from 0)
 * @param s2  S2 value
 * @param band band index (counting from 0)
*/
    void add(double s1,int slice, double s2,int band);

/**
 * Add content of a Tgraph
 * @param gr pointer to the graph whose values are S1 and log(S2/S1)
 */
    void add(TGraph* gr);

/**
 * Add content of another data set after applying a set of cuts
 * @param original pointer to the original data set
 * @param cuts pointer to the set of cuts
 */
    void add(S1S2Data* original, XeSetOfCuts* cuts);

/**
 * Add content of another data set after applying a dark matter cut
 * @param original pointer to the original data set
 * @param cuts pointer to the Dark matter cuts
 */
    void add(S1S2Data* original, DarkMatterCut* cut);

/**
 * Create from content of a Tgraph
 * @param gr pointer to the graph whose values are S1 and log(S2/S1)
 */
    void make(TGraph* gr);

/**
 * Create from another data set after applying a set of cuts
 * @param original pointer to the original data set
 * @param cuts pointer to the set of cuts
 */
    void     make(S1S2Data* original, XeSetOfCuts* cuts);

/**
 * crate from another data set after applying a dark matter cut
 * @param original pointer to the original data set
 * @param cuts pointer to the Dark matter cuts
 */
    void     make(S1S2Data* original, DarkMatterCut* cut);

/**
 * Define the corresponding data set after cuts
 * @parm cut pointer to the S1S2Data data set after cuts
 */ 
    void     setCutData(S1S2Data* cut);

/** 
 * Draw the data set, as requested by this being S1S2Object
 */
    void     draw();

/**
 * Reset (empty all contents)
 */
    void     reset();
                    
/**
 *  Print the data set with the 1st bands
 *  @param level print level
 * @param  bandmax index of highest band to be printed (default: no limit)
 */
    void     printComputedBands(int level=1, int bandMax=ALL);

/**
 * Print a summary of the content , befrore and after cut, and normalisation
 * @param header Entry title in the table
 */
    void     printSummary(string header);

/** 
 * Return the data type
 * @return One of : AM_BE_DATA, E_GAMMA_DATA , DM_DATA , AM_BE_CUT_DATA  
 * , E_GAMMA_CUT_DATA, DM_CUT_DATA, ER_BACKGROUND, ER_LEAKAGE    
 * , NR_BACKGROUND, ALL_BACKGROUNDS, DM_CUT_DATA_INTEGRATED, 
 * ALL_BACKGROUNDS_INTEGRATED , DEFAULT_SIGNAL, DEFAULT_SIGNAL_INTEGRATED
 */
    int      getDataType();

/**
 * Return the default color while drawing
 */
    int      defaultColor();

/**
 * Return number of events in the data set
 */
    int      getNEvents();

/**
 * Return band index of a given event (if defined)
 * @return index of band, counted from zero
 * @param ev event number ,counted from zero
 */
    int      getBand(int ev);

/**
 * Return slice index of a given event (if defined)
 * @return index of slice, counted from zero
 * @param ev event number ,counted from zero
 */
    int      getSlice(int ev);

/**
 * Get normalization factor
 */
    double   getNormalization();

/**
 * Get normalized number of events
 */
    double   getNEventsNormalized();

/**
 * Return S1 value of a given event
 * @param ev event number ,counted from zero
 */
    double   getS1(int ev);

/**
 * Return S2 value of a given event
 * @param ev event number ,counted from zero
 */
    double   getS2(int ev);

/**
 * Return S2/S1  of a given event
 * @param ev event number ,counted from zero
 */
    double   getS2overS1(int ev);

/**
 * get "X" of a given event, when drawing
 * @param ev event number ,counted from zero
 */
    double   getX(int ev);

/**
 * get "Y" of a given event, when drawing
 * @param ev event number ,counted from zero
 */
    double   getY(int ev);

/**
 * get max "X" value
 */
    double   maxX();

/**
 * get max "Y" value
 */
    double   maxY();

/**
 * get max S1
 */
    double   maxS1();

/**
 * get max S2
 */
    double   maxS2();

/**
 * get max S2/S1
 */
    double   maxS2OverS1();

/**
 * get min "X" value
 */
    double   minX();

/**
 * get min "Y" value
 */
    double   minY();

/*
 * get min S1
 */
    double   minS1();

/**
 * get min S2
 */
    double   minS2();

/**
 * get min S2/S1
 */
    double   minS2OverS1();

/**
 * get all S1 values
 * @return  pointer to array
 */
    double*  getS1();

/** 
 * get all S2 values
 * @return  pointer to array
 */
    double*  getS2();

/**
 * get name of data type
 */
    string   getTypeName();

/**
 * get corresponding data set after cut
 * @return pointer to the S1S2Data XeObject
 */
    S1S2Data *getCutData();

/**
 * Create a graph
 * @return pointer to the newly created XeGraph XeObject
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph  *newGraph(int plot=NONE);

/**
 * Create a TTree with important quantities in.
 * They are S1,S2, logS2S, flattened.
 * If reference bands are given, slice and band are also given
 * @return a pointer to a newly created TTree
 * @param referenceBands  Referance bands used to compute slice and band numver
 */
    TTree    *newTree(S1S2Bands* referenceBands=NULL);

/**
 * Return S1 bins with equal content, with automatic S1 edges
 * @return pointer to a newly created EquiContentBins XeBins XeObject
 * @param nS1 number of S1 slices
 */
    EquiContentBins  *newEquiContentS1Bins(int ns1);

/**
 * Return S1 bins with equal content, giving S1 edges
 * @return pointer to a newly created EquiContentBins XeBins XeObject
 * @param nS1 number of S1 slices
 * @param s1Mi Lower edge of S1 range
 * @param s1Ma upper esge of S1 range
 */
    EquiContentBins  *newEquiContentS1Bins(int ns1,double s1Mi,double s1Ma);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
/**
 * Print header of breakdown table
 */
    static void printSummaryHeader();

/**
 * get default color, for a given data type
 * @param  type Data type
 */
    static int  defaultColor(int type);

/**
 * tell if a data type comes from another one after cuts
 * @param  type Data type
 */
    static bool isAfterCuts(int type);

/**
 * given a data type , return the corresponding type after cuts
 * @param  type Data type
 */
    static int typeAfterCut(int type);

/**
 * given a data type , return the corresponding type before cuts
 * @param  type Data type
 */
    static int typeBeforeCut(int type);

/**
 * Is the type a background one?
 * @param  type Data type
 * @param warn print a warning statement in case of error
 */
  static bool isBackground(int type, bool warn=false);

/**
 * return type name
 * @param  type Data type
 */
    static string getTypeName(int type);

/**
 * return short type type name
 * @param  type Data type
 */
    static string getShortTypeName(int type);

/**
 * format count  for a table
 * @param  type Data type
 */
   static string formatContent(int type,double c,double w, bool total=false);

/**
 * Map S1shape enumeration to S1S2dataType enumeration
 * @param shape either SHAPE_SIGNAL,SHAPE_ER_BACKGROUND,SHAPE_ER_LEAKAGE,SHAPE_NR_BACKGROUND
 * @return the corresponding item S1S2 data type
 */
 static int mapS1Shape(int s1s2Shape);

/**
 * Update min, max, ...
 */
    bool  update();

/**
 * Internal method to help the constructor by setting vital members
 * @param run XeRun XeObject this data set belongs to
 * @param type Data type
 */
    void              setParameters(int type);
 

/**  
 *  print the data set
 *  @param level print level
 */
    bool     printIt(int level=1);


    static string     getTheName(S1S2Data* original, XeSetOfCuts* cuts);
    static string     getTheName(S1S2Data* original, DarkMatterCut* cut);

    bool              passes(int i, DarkMatterCut* cut);
    bool              passes(int i, XeSetOfCuts* cuts);

  protected:
    int               dataType;
    double            normalization;
    double            nEventsNormalized;
    double            s1Min;
    double            s1Max;
    double            s2Min;
    double            s2Max;
    double            s2overs1Min;
    double            s2overs1Max;
    S1S2Data*         cutData;
    vector<int>       bands;
    vector<int>       slices;

    float             S1;
    float             S2;
    float             logS2S1;
    float             flattened;
    int               band;
    int               slice;

      ClassDef(S1S2Data,1)

};


/**
   * A Slice (vertical cut in S1).
   * Defined by slice number, and characterized by S1Min, S1Max
*/

class S1Slice : virtual public DetectorComponent, virtual public S1S2Object
              , virtual public RunComponent  {

  public :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
 /**
    * Regulat constructor from an S1S2Bands XeObject
     * @param  w which slice
     * @param  s1s2b original S1S2Bands
 */
                    S1Slice(int w,S1S2Bands *s1s2b);
                    S1Slice();
                   ~S1Slice();

/**
 * Draw it, being and S1S2Object object
 */
    void  draw();

/**
 * Empty content
 */
    void  reset();
                    
/** 
 * Set S1 minimum value
 */
    void  setS1Min(double s1Min); 

/** 
 * Set S1 maximum value
 */
    void  setS1Max(double s1Max); 

 /**
     * fill one entry in a given band
     * Warning! Should not be called independently, to avoir mismatch 
     *          between bands and slices
     * @param s1  S1 value
     * @param s2  S2 value
     * @param  w   weight
 */
    void fill(double s1,double s2,double w=1.);
 /**
     * fill one entry in a given band
     * Warning! Should not be called independently, to avoir mismatch 
     *          between bands and slices
     * @param bin band index
     * @param  w   weight
*/
    void fill(int band,double w=1.);

 /**
     * Fill the S1overS2 table for study (i.e. define the bands)
     * @param s1  S1 value
     * @param s2  S2 value
 */
    void fillForStudy(double s1, double s2);

/**
 * Normalize (i.e. multiply by a factor)
 * @param ration multiplicatoin factor
 */
    void normalize(double ratio);

 /**
     * reweights  contents inside slices
     * @param   weights relative weight of each band (sum[weights]=1);
 */
    void reweightTheBands(double * weights);

/**
 * Return total count
 */
    double getContent(); 

/**
 * Return count in given band
 */
    double getContent(int band); 

/**
 * Return maximum count in all cells
 */
    double getMaxCount(); 

/**
 * Return minimum count in all cells
 */
    double getMinCount(); 

/**
 * Return central value of S1
 */
    double getCentralS1(); 

/**
 * Return minimum value of S1
 */
    double getS1Min(); 

/**
 * Return maximum value of S1
 */
    double getS1Max(); 

/**
 * Return width of S1 interval
 */
    double getS1Width(); 

/**
 * Return smallest S2overS1 bin
 */
    double getS2overS1Min(); 

/**
 * Return largest  S2overS1 bin
 */
    double          getS2overS1Max(); 

/**
 * Return area in s1*log(S2/S1) space
 */
    double          getArea(); 
    
/**
 * Return list of S2/S1 bins
 */
    vector<double>* getS2overS1Vector();
    
/**
 * Return  S2/S1 bins as an XeBins* XeObject
 */
    XeBins*         getS2overS1Bins();

/**
 * Return maximum 'X' value, as requested by being an S1S2Object object
 */
    double        maxX();

/**
 * Return maximum 'Y' value, as requested by being an S1S2Object object
 */
    double        maxY();

/**
 * Return minimum 'X' value, as requested by being an S1S2Object object
 */
    double        minX();

/**
 * Return minimum 'Y' value, as requested by being an S1S2Object object
 */
    double        minY();


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    void            setS2overS1Bins(XeBins* s2overS1Bins);
    static string   getTheName(int s,S1S2Bands *s1s2b);

/**  
 *  print the S1 Slice with minimum and maximum s1/s2 bands
 *  @param level print level
 */
    bool  printIt(int level=1);
/**  
 *  print the S1 Slices with the minimum and maximum LOG10(s2/s1) bands
 *  @param level print level
 */
    bool  printItLog(int level=1);

  protected :


    int             sequence;      /*!< slice sequence number, counting from 0*/
    int             nBands;                               /*!< number of bands*/
    double          S1Min;                            /*!< Minimum value of S1*/
    double          S1Max;                            /*!< Maximum value of S1*/
    double          count;                                    /*!< Total count*/

    S1S2Bands      *mother;            /*!< Pointer to S1S2Bands it belongs to*/
    XeBins         *S2overS1Bins;      /*!< Bands limits in S2/S1, as a XeBins*/
    vector<double>  S2overS1Values;            /*!<Values of S2/S1 for studies*/
    vector<double>  counts;                              /*!< getContents per band */


      ClassDef(S1Slice,1)

} ;

/**
   * A band (horizontal projection on S2)
*/

class S2Band : virtual public DetectorComponent, virtual public S1S2Object
             , virtual public RunComponent {

  public :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
 /**
     * Regular constructor from a S1S2Bands XeObkect
     * @param  w which ban
     * @param  s1s2b original S1S2Bands
 */
    S2Band(int w,S1S2Bands *s1s2b);

/**
 * Empty constructor for ROOT
 */
    S2Band();
    ~S2Band();

/**
 * Draw it, being and S1S2Object object
 */
    void  draw();

/**
 * Draw one of its edges
 * @param whichEdge  either UPPER_EDGE or LOWER_EDGE
 */
    void  draw(int whichEdge);

/**
 *  Draw contents of cells in color code
 *  @param cMin  content corresponding to lowest color
 *  @param cMax  content corresponding to highest color
 * @param zMode either LINEAR or LOG
 */
    void  drawCellContent(double cMin,double cMax,int zMode);

/**
 *  Draw band content in color code
 *  @param cMin  content corresponding to lowest color
 *  @param cMax  content corresponding to highest color
 * @param zMode either LINEAR or LOG
 */
    void  drawBandContent(double cMin,double cMax,int zMode);

/** 
 * Extend the lower and upper edges of S1 Range
 */
    void  extendS1Range(double s1Min, double s1Max);

 /**
     * fill one entry in a given S1 bin.
     * Warning! Should not be called independently, to avoir mismatch 
     *          between bands and slices
     * @param  s1  S1 value (which will be stored independently)
     * @param  w   weight
 */
    void  fill(double s1,double w=1.);

 /**
     * fill one entry in a given S1 bin.
     * Warning! Should not be called independently, to avoir mismatch 
     *          between bands and slices
     * @param  bin  s1 slice
     * @param  w   weight
 */
    void  fill(int bin,double w=1.);

/**
 * Empty content
 */
    void reset();
                    
/**
 * Print S1 values
 * @param maxEntries  number of entries, 0--> all
 */
    void printS1(int maxEntries=0);

/**
 * Normalize (i.e. multiply by a factor)
 * @param ration multiplicatoin factor
 */
    void normalize(double ratio);

/** 
 * fill a band uniformally
 * @param w total content of the band
 */
    void fillUniformly(double w);

/** 
 * fill a band according to a S1 distribution
 * @param w total content of the band
 * @param dist pointer to table of distributions
 */
    void fillAccordingToS1Dist(double w, double* dist);

/**
 * Set the vectors describing S2overS1 limits
 * @param s2overS1Min pointer to the vector of min s2overS1
 * @param s2overS1Max poaxter to the vector of max s2overS1
 */
    void setS2overS1Limits( vector<double>* s2overS1Min
                                   , vector<double>* s2overS1Max );
/**
 * Get the band sequence number
 * @return index running from zero
 */
    int  getSequence(); 
   
/** 
 * Get total content
 */
    double getContent(); 
   
/**
 * Get cell content, given the slice index
 * @param  slice  slice index, starting from zero
 */
    double getContent(int slice); 
   
/**
 * Get cell content, given the value of S1
 * @param  s1 value of s1
 */
    double getContent(double s1); 

/**
 * return pointer to counts
 */

   vector<double>* getContents();

/**
 * Get maximum cell count
 */
    double getMaxCount(); 

/**
 * Get minimum cell count
 */
    double getMinCount(); 

/**
 * Get minimum count in not empty cells
 */
    double getMinCountButZero();

/**
 * Get overall lowest S1 edge 
 */
    double getS1Min();

/**
 * Get overall lowest S1 edge 
 */
    double getS1Max();

/**
 * Return maximum 'X' value, as requested by being an S1S2Object object
 */
    double maxX();

/**
 * Return maximum 'Y' value, as requested by being an S1S2Object object
 */
    double maxY();

/**
 * Return minimum 'X' value, as requested by being an S1S2Object object
 */
    double minX();

/**
 * Return minimum 'Y' value, as requested by being an S1S2Object object
 */
    double minY();

/**
 * Get overall upper edge of S2/S1
 */
    double getUpperEdge();

/**
 * Get upper edge of S2/S1 for a given slice
 * @param slice slice index, starting from 0
 */
    double getUpperEdge(int slice);

/**
 * Set upper edge of S2/S1 for a given slice
 * @param slice slice index, starting from 0
 */
    void setUpperEdge(int slice, double UpEdge);

/**
 * Set lower edge of S2/S1 for a given slice
 * @param slice slice index, starting from 0
 */
    void setLowerEdge(int slice, double UpEdge);

/**
 * Get overall lower edge of S2/S1
 */
    double getLowerEdge();

/**
 * Get lower edge of S2/S1 for a given slice
 * @param slice slice index, starting from 0
 */
    double getLowerEdge(int slice);

/**
 * Get overall S1 width
 */
    double getWidth();

/**
 * Get  S1 width of a given slice
 * @param slice slice index, starting from 0
 */
    double getWidth(int slice);

/**
 *  return list of S2/S1 upper edges
 *  @return a pointer to the array
 */
    vector<double> getUpperEdges();

/**
 *  return list of S2/S1 lower edges
 *  @return a pointer to the array
 */
    vector<double> getLowerEdges();

/**
 * return list of stored S1 values
 *  @return a pointer to the array
 */
    vector<double>* getS1();

/**
 * return list of stored S1 bins, if stored from data returns the corresponding of "getS1" but 
 * in the relative S1 bin number and not S1 value.
 *  @return a pointer to the array
 */
    vector<int> getS1InBin();

 /**
     * Create the spectrum of Slice values
     * @return a pointer to the newly created XeSpectrum
     * @param  nb    number of bins (default: automatic)
     * @param  s1Min lower bound of Slice (default: automatic)
     * @param  s1Max upper bound of Slice (default: automatic)
 */
    XeSpectrum* newSliceXeSpectrum();
 /**
     * Create the distribution (i.e. normalized spectrum) of Slice values
     * @return a pointer to the newly created TabulatedDist
     * @param  nb    number of bins (default: automatic)
     * @param  s1Min lower bound of Slice (default: automatic)
     * @param  s1Max upper bound of Slice (default: automatic)
*/
     TabulatedDist* newSliceXeDist();

/**
 * Create a one-dim histogram of Slice Spectrum
 * @return pointer to newly created histogram
 * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
 */
  TH1F* newHistogramOfSliceSpectrum(int plot=NONE);

/**
 * Create a one-dim histogram of Slice Distribution
 * @return pointer to newly created histogram
 * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
 */
  TH1F* newHistogramOfSliceDistribution(int plot=NONE);

 /**
 * Create the spectrum of S1 values
 * @return a pointer to the newly created XeSpectrum
 * @param  nb    number of bins (default: automatic)
 * @param  s1Min lower bound of S1 (default: automatic)
 * @param  s1Max upper bound of S1 (default: automatic)
 */
    XeSpectrum* newS1XeSpectrum(int nb=AUTO
              ,double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);
 /**
 * Create the distribution (i.e. normalized spectrum) of S1 values
 * @return a pointer to the newly created TabulatedDist
 * @param  nb    number of bins (default: automatic)
 * @param  s1Min lower bound of S1 (default: automatic)
 * @param  s1Max upper bound of S1 (default: automatic)
*/
     TabulatedDist* newS1XeDist(int nb=AUTO
            ,double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);

/**
 * Create a one-dim histogram of S1 Spectrum
 * @return pointer to newly created histogram
 * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
 * @param  nb    number of bins (default: automatic)
 * @param  s1Min lower bound of S1 (default: automatic)
 * @param  s1Max upper bound of S1 (default: automatic)
 */
  TH1F* newHistogramOfS1Spectrum(int plot=NONE,int nb=AUTO
                      ,double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);

/**
 * Create a one-dim histogram of S1 Distribution
 * @return pointer to newly created histogram
 * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
 * @param  nb    number of bins (default: automatic)
 * @param  s1Min lower bound of S1 (default: automatic)
 * @param  s1Max upper bound of S1 (default: automatic)
 */
  TH1F* newHistogramOfS1Distribution(int plot=NONE, int nb=AUTO
                      ,double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);



    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   static string  getTheName(int s,S1S2Bands *s1s2b);

 /**  
 *  print the S2 Band
 *  @param level print level
 */
    bool printIt(int level=1);

  protected :


   S1S2Bands     *mother;           /*!<pointer to the S1S2Bands it belongs to*/
   XeBins        *s1Bins;                        /*!<Interval defining S1 bins*/
   double         count;                                       /*!<Total count*/
   int            sequence;               /*!<Sequence number , counted from 0*/
   int            nSlices;                                /*!<Number of slices*/
   vector<double> S1;     /*!<Stored values of S1 (does not necessarily exist)*/
   vector<double> counts;                                 /*!<counts per slice*/
   vector<double> S2overS1Min;                       /*!<lower limits if S2/S1*/
   vector<double> S2overS1Max;                       /*!<upper limits if S2/S1*/

      ClassDef(S2Band,1)


};

/**
   *  Bands and Slices maintained together
*/

class S1S2Bands : virtual public DetectorComponent, virtual public S1S2Object 
                , virtual public RunComponent {

  public :
   
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    

    ~S1S2Bands();
     S1S2Bands();

 /**
     * Constructor with equicontent S1 slices 
     * @param  nam   name of the object
     * @param  run run (actually XeRun) where it's taken from
     * @param  nBands number of bands 
     * @param  s1s2   original S1S2data
     * @param  nS1 nimber of slices in S1
 */
     S1S2Bands(string nam, XeRun *run,int nBands, S1S2Data* s1s2, int nS1);

 /**
     * Constructor with equidistant S1 slices 
     * @param  n   name of the object
     * @param  run run (actually XeRun) where it's taken from
     * @param  nBands number of bands 
     * @param  nS1  number of slices in S1
     * @param  s1Mi minimum S1
     * @param  s1Ma maximum S1
 */
     S1S2Bands(string n,XeRun *run,int nBands,int nS1,double s1Mi,double s1Ma);

 /**
     * Constructor giving description of s1 slices
     * @param  n   name of the object
     * @param  run run (actually XeRun) where it's taken from
     * @param  nBands number of bands 
     * @param  s1Bins bin interval (actually XeBins) describing slices
 */
     S1S2Bands(string n, XeRun *run,int nBands, XeBins* s1Bins);

 /**
     * Constructor copying from another S1S2Bands
     * @param  bands  Orignal bands to copy
     * @param  name   The name of the new S1S2Bands
     * @param  fillIt copy the contents
 */
     S1S2Bands(string name, S1S2Bands* bands,bool copyContent=false);

 /**
     * Constructor copying from another S1S2Bands
     * @param  bands  Orignal bands to copy
     * @param  name   The name of the new S1S2Bands
     * @param  fillIt copy the contents
 */
     S1S2Bands(string name, S1S2Bands bands,bool copyContent=false);
  
/**
 * Is such an object drawable in a given S1S2 Display mode ?
 * @param mode S2_VS_S1, S2_OVER_S1_VS_S1, FLATTENED_S2_VS_S1 ,BAND_VS_SLICE
 */
  virtual bool isDrawable(int mode);


 /**
     * Draw it (it's drawable in a S1S2 plot)
 */
     void draw();
 /**
  * Draw it (it's drawable in a S1S2 plot)
  * @param zMode either LINEAR or LOG
 */
     void  draw(int zMode);

 /**
     * Density of projected S1 (in 1/p.e.)
     * @param  s1 value at which it is evaluated 
 */
     double  getS1Density(double s1);

  /**
     * Fill an entry, (actually fill twice: one slice and one band) 
     * @param  s   slice number
     * @param  b   band number
     * @param  w   weight
 */
     void  fill(int s,int b,double w=1.);
 /**
     * Fill an entry, (actually fill twice: one slice and one band).
     * Also fills the table of S1 in the proper band
     * @param  s1  S1 Value at which it's filled
     * @param  s   slice number
     * @param  b   band number
     * @param  w   weight
 */
     void  fill(double s1, int s,int b,double w=1.);

 /**
     * Fill all bands equally for a given S1
     * @param  s1  Value at which it's filled
     * @param  w   total weight(i.e. the weight in each band is w/nBands)
 */
     void  dispatchS1InEqualBands(double s1, double w);

 /**
     * Fill all bands with density constant in S1, with variable weight
     * @param  weights  weight in each of the bands   
 */
     void  fillUniformly(vector<double>& weights);

 /** 
 * fill all bands according to the same S1 distribution
 * @param  weights  weight in each of the bands   
 * @param dist pointer to table of distributions
 */
    void fillAccordingToS1Dist(vector<double>& weights, double* dist);

 /**
     * Fill all bands with density exponentially decreasing in S1
     * with variable weight
     * @param  weights  weight in each of the bands   
     * @param  slope  parameter of the slope
 */
     void  fillExponential(vector<double>& weights,double slope);

  
/**
 * Return number of bands
 */
     int getNBands();

/**
 * Return number of Slices
 */
     int getNSlices();

/**
 * Return getContent	of a cell
 * @param s slice number (counted from 0)
 * @param b band number (counted from 0)
 */
     double getContent(int s, int b);

/**
 * return total count
 */
     double getTotalContent();

/**
 * return total count below a given s1
 * @param s1 upper limit
 */
     double getContentBelowS1(double s1);

/**
 * return count in a given band
 * @param b  band number (counted from 0)
 */
     double getBandContent(int b);

/**
 * return max count in a cell
 */
     double getMaxCount(); 

/**
 * return min count in a cell
 */
     double getMinCount(); 

/**
 * return max count in a band
 */
     double getMaxBandContent(); 

/**
 * return min count in a band
 */
     double getMinBandContent(); 

/**
 * return max X = S1 max (this method must be implemented)
 */
     double maxX();

/**
 * return min count in a band
 * return max Y = S2/S1 max (this method must be implemented)
 */
     double maxY();

/**
 * return min X = S1 min (this method must be implemented)
 */
     double minX();

/**
 * return min Y = S2/S1 min (this method must be implemented)
 */
     double minY();

/**
 * return total S1 range width
 */
     double getS1Width();

/**
 * return S1 central value in a given slice 
 * @param s slice number, counted from zero
 */
     double getS1Center(int s);

/**
 * return s1 bin width of a given slice (usesul if slices not equidistant!)
 * @param s slice number, counted from zero
 */
     double getS1Width(int s);

/**
 * return overall S1 max
 */
     double getS1UpperEdge();

/**
 * return  S1 max in a given slice 
 *  * @param s slice number, counted from zero
 */
     double getS1UpperEdge(int s);

/**
 * return overall S1 min
 */
     double getS1LowerEdge();

/**
 * return  S1 min in a given slice 
 * @param s slice number, counted from zero
 */
     double getS1LowerEdge(int s);

 /**
     * return  Total area in delta S1 * delta log(s2/s1)
 */
     double getArea(); 

 /**
     * return list of S1 upper edges, for all slices
     * @return pointer to the list 
 */
     double* getS1UpperEdges();

 /**
     * return list of S1 lower edges, for all slices
     * @return pointer to the list 
 */
     double* getS1LowerEdges();

/** 
 * return overall upper S1/S1 in a given band
 * @param b band number counted from zero
 * */
     double  getS2overS1UpperEdge(int b);

/** 
 * return upper S1/S1 in a given cell
 * @param b band number counted from zero
 * @param s slice number, counted from zero
 * */
     double getS2overS1UpperEdge(int b,int s);

/** 
 * return list of overall upper S1/S1 in a given band
 * @param b band number counted from zero
 * @return pointer to the list 
 * */
     vector<double> getS2overS1UpperEdges(int b);

/** 
 * return overall lower S1/S1 in a given band
 * @param b band number counted from zero
 * */
     double  getS2overS1LowerEdge(int b);

/** 
 * return lower S1/S1 in a given cell
 * @param b band number counted from zero
 *  * @param s slice number, counted from zero
 * */
     double getS2overS1LowerEdge(int b,int s);

/** 
 * return list of overall lower S1/S1 in a given band
 * @param b band number counted from zero
 * @return pointer to the list 
 * */
     vector<double> getS2overS1LowerEdges(int b);

/**
 * Return the S1 bin description
 * @return pointer to the XeBins XeObject
 */ 
     XeBins* getS1Bins();         

/**
 * Return a given band
 * @param b band number counted from zero
 * @return pointer to the S2Band XeObject
 */
     S2Band* getBand(int b); 

/**
 * Return a given slice given its index
 * @param s slice number, counted from zero
 * @return pointer to the S1Slice XeObject
 */
     S1Slice* getSlice(int s); 

/**
 * Return a given slice for a given S1
 * @param s1 value if s1
 * @return pointer to the S1Slice XeObject
 */
     S1Slice* getTheSlice(double s1); 


/**
 * Add the the given S1S2Bands to the caller
 * @param toBeAdded  bands to be added
 */

     void add(S1S2Bands *toBeAdded);

/**
 * Define the bands, given a reference data set.
 * Returns empty bands.
 * @param reference  pointer to reference S1S2Data data set 
 */
     void  compute(S1S2Data* reference);

/**
 * Define the bands, given a reference TH2F histogram.
 * The histogram is meant to be in log(S2/S1) VS S1 space.
 * The range in which bands are defined in log(S2/S1) is (default) [1.3,2.7]
 * Returns empty bands if problems occur.
 * @param reference histogram 
 */
     void  compute( TH2F reference, bool isLog = true);

/**
 * Print it band by band
 * @param level verbosity
 */ 
     void  printBands(int level=1);

/**
 * Print a few S1 values in each band
 * @param maxEntries  how many entries per band
 */ 
     void  printS1InBands(int maxEntries=0);

/**
 * Print it slice by slice
 * @param level verbosity
 */ 
     void  printSlices(int level=1);

/**
 * Draw it slice by slice
 */ 
     void  drawSlices();

/**
 * method to override drawFrame
 */
  void   drawTheFrame();

/**
 * Draw it band by band
 */ 
     void  drawBands();

/**
 * Draw the frame for cell by cell drawing
 * @param zMode either LINEAR or LOG
 */ 
     void  drawCellContentFrame(int zMode);

/**
 * Draw  cell by cell
 * @param zMode either LINEAR or LOG
 */ 
     void  drawCellContent(int zMode);


/**
 * Draw  cell by cell, with a frame
 * @param zMode either LINEAR or LOG
 */ 
     void  drawCellContentWithFrame(int zMode);


/**
 * Draw the frame for full band drawing
 * @param zMode either LINEAR or LOG
 */ 
     void  drawBandContentFrame(int zMode);

/**
 * Draw  full band content
 * @param zMode either LINEAR or LOG
 */ 
     void  drawBandContent(int zMode);

/**
 * Draw  full band content with fra,e
 * @param zMode either LINEAR or LOG
 */ 
     void  drawBandContentWithFrame(int zMode);

/**
 * produce a summary page
 */ 
     void produceGraphicOverview();

 /**
     * create  a 2d histogram whose content is follows slices and bands
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
     * @return  pointer to the newly  created  TH2F
 */
    TH2F* new2DHistogram(int plot=NONE);

 /**
 * create  a multigraph of Slice Spectrum
 * @return  pointer to the newly created  XeMultiGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     XeMultiGraph* newMultiGraphOfSliceSpectrum(int plot=NONE);

 /**
 * get a multigraph of Slice Distribution (normalized to 1. in each band)
 * @return  pointer to the newly created  XeMultiGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     XeMultiGraph* newMultiGraphOfSliceDistribution(int plot=NONE);



 /**
     * create a spectrum of Band content,
     * @return a pointer to  the newly created XeSpectrum
     * @param band   band index, starting at 0; default: ALL= cumulated
 */
    XeSpectrum* newBandXeSpectrum();

 /**
     * create a distribution of Band content 
     * @return a pointer to the newly created TabulatedDist
     * @param band   band index, starting at 0; default: ALL= cumulated
*/
     TabulatedDist* newBandXeDist();

 /**
     * Create an histogram of the spectrum of Band content
     * @return a pointer to the newly created histogram
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
  */
    TH1F* newHistogramOfBandSpectrum(int plot=NONE);

 /**
     * Create an histogram of the distribution(spectrum norm. to 1) of Band content
     * @return a pointer to the newly created histogram
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
  */
    TH1F* newHistogramOfBandDistribution(int plot=NONE);
 /**
   * simple implementation that returns histo filled with S1 distro for the integral of the specified bands
   * @param from_band  band index from which integrate
   * @param to_band    band index to which integrate, if not specified is the same band 
 */
    TH1F getHistogramOfS1Distro(int from_band, int to_band = 999);

 /**
    Returns TGraph of bands
  */ 
     vector<TGraph> getTGraphOfBands(bool isS2_vs_S1=true); 

 /**
     * create a spectrum of Slice content,
     * @return a pointer to  the newly created XeSpectrum
     * @param band   band index, starting at 0; default: ALL= cumulated
 */
    XeSpectrum* newSliceXeSpectrum(int band=ALL);

 /**
     * create a distribution of Slice content 
     * @return a pointer to the newly created TabulatedDist
     * @param band   band index, starting at 0; default: ALL= cumulated
*/
     TabulatedDist* newSliceXeDist(int band=ALL);

 /**
     * Create an histogram of the spectrum of Slice content
     * @return a pointer to the newly created histogram
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
  */
    TH1F* newHistogramOfSliceSpectrum(int band=ALL,int plot=NONE);

 /**
     * Create an histogram of the distribution(spectrum norm. to 1) of Slice content
     * @return a pointer to the newly created histogram
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
  */
    TH1F* newHistogramOfSliceDistribution(int band=ALL,int plot=NONE);

 /**
     * create a spectrum of S1 values ,
     * @return a pointer to  the newly created XeSpectrum
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param  nb    number of bins (default: automatic)
     * @param  s1Min lower bound of S1 (default: automatic)
     * @param  s1Max upper bound of S1 (default: automatic)
 */
    XeSpectrum* newS1XeSpectrum(int band=ALL, int nb=AUTO
              , double s1Min=AUTOMATIC,double S1Max=AUTOMATIC);
 /**
     * create a distribution of S1 values 
     * @return a pointer to the newly created TabulatedDist
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param  nb    number of bins (default: automatic)
     * @param  s1Min lower bound of S1 (default: automatic)
     * @param  s1Max upper bound of S1 (default: automatic)
*/
     TabulatedDist* newS1XeDist(int band=ALL, int nb=AUTO
                  , double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);

 /**
     * Create an histogram of the spectrum of S1 values
     * @return a pointer to the newly created histogram
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
     * @param  nb    number of bins (default: automatic)
     * @param  s1Min lower bound of S1 (default: automatic)
     * @param  s1Max upper bound of S1 (default: automatic)
  */
    TH1F* newHistogramOfS1Spectrum(int band=ALL,int plot=NONE, int nb=AUTO
                  , double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);

 /**
     * Create an histogram of the distribution(spectrum norm. to 1) of S1 values
     * @return a pointer to the newly created histogram
     * @param band   band index, starting at 0; default: ALL= cumulated
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
     * @param  nb    number of bins (default: automatic)
     * @param  s1Min lower bound of S1 (default: automatic)
     * @param  s1Max upper bound of S1 (default: automatic)
  */
    TH1F* newHistogramOfS1Distribution(int band=ALL,int plot=NONE, int nb=AUTO
                  , double s1Min=AUTOMATIC ,double s1Max=AUTOMATIC);




/**
 * Reset all contents
 */
     void reset();
              
/**
 * Extend the S1 range of
 * @param s1Min new and smaller lower S1 bound
 * @param s1Max new and larger upper S1 bound 
 */ 
     void extendS1Range(double s1Min, double s1Max);

/** 
 * Return the cumulated band count, from band 0 to (including) the specified one
 * @param band index, couting from zero
 */
     double getCumulatedBandContent(int b);

/**
 * Return the smallest overall non zero cell count
 */ 
     double getMinCountButZero(); 

/**
 * Return the smallest overall non zero band count
 */ 
     double getMinBandContentButZero(); 

 /**
     * Multiply everything by a constant factor 
     * @param   ratio  multiplication factor 
 */
     void normalize(double ratio);

 /**
     * normalize the overall content to "events"
     * @param   events the total sum in the bands after normalization
 */
     void normalizeToEvents(double events);

/**
 * Print 
 * @param level verbosity
 */ 
     bool  printIt(int level=1);

/**
 * Return band and slice numbers for a given s1,s2
 * @param s1  value of S1
 * @param s2 value of S2
 * @return a pair containing a band number and a slice number
 */
     pair<int,int>    whichSliceAndBand(double s1,double s2);

/** 
 * Fill with the cumulated (intgrated) content of an original S1S2Bands
 * @param original The original one
 */
     void  cumulate(S1S2Bands* original);

/** 
 * return overall minimum band size in s2/s1
  @param if isLogRequired = true returns the minimum in log(s2/s1)
 */
     double getMinimumBandSize( bool isLogRequired = false);


  protected :

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
      
              
/** 
 * Fill slices after bands were filled
 */
   void copyBandsToSlices();

/**
 * Establish the table of pointers to band counts vector
 */
   void establishPointersToCounts();

/**
 * Return a default name
 */
     static string    getTheName(int nBands,int nSlices);

 /**
     * reweights the whole bands and contents inside slices
     * @param   weights relative weight of each band (sum[weights]=1);
 */
     void reweightTheBands(double * weights);


/**
 * Check that band number either is valid or is equal to 'ALL'
 * @param band index to be checked, or ALL
 * @return a pair containing the fist and last indices, or NONE if error
 */
     pair<int,int>    whichBands(int band);

/**
 * return the name of a given (or all) bands
 * @param band index to be checked, or ALL
 */
   string bandName(int band);
 
 /**
 * Internal method: create a multigraph from nBands spectra 
 * @return pointer to newly created multigraph (actually XeMultiGraph)
 * @param  name    name of the multigraph to be created
 * @param  spectra  pointer to pointers of spectra
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     XeMultiGraph* newMultiGraphOfSlice(string name, XeSpectrum** spectra
                    , int plot=NONE);

 /**
 * Internal method: create a multigraph from nBands distributions 
 * @return pointer to newly created multigraph (actually XeMultiGraph)
 * @param  name    name of the multigraph to be created
 * @param  dists  pointer to pointers of disitribution
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     XeMultiGraph* newMultiGraphOfSlice(string name, TabulatedDist** dists
                    , int plot=NONE);
 
/** 
 * Internal method: instantiate bands and slices.
 * Assume that the s1Bins XeBins object already exists
 * @param name    Name to be passed
 * @param nBands Number of bands
 */
     void   instantiateBandsAndSlices(string name, int nBands);

/** 
 * Internal method: name the bands and the slices
 * @param name    Name to be passed
 */
     void   setNames(string name);

     int              nBands;                             /*!< number of bands*/
     int              nSlices;                           /*!< number of slices*/
     XeBins          *s1Bins;      /*!< pointer to description of the S1 bins */
     vector<S2Band*>  bands;                         /*!< pointer to the bands*/
     vector<S1Slice*> slices;                       /*!< pointer to the strips*/
     vector<vector<double>*> counts;      /*!< pointers to counts in each band*/

     TF1 *fBottom;
     
     ClassDef(S1S2Bands,1)

} ;

/**
   * Class for background model. Either neutron or electron
*/
class XeBackgroundModel : virtual public DetectorComponent
                        , virtual public RunComponent {

   public :


 /**
 * return  Bands for (regular) background
 */
     virtual  S1S2Bands* computeBands()=0;

/**  
 *  print the background model
 *  @param level print level
 */
     virtual ~XeBackgroundModel();

 /**
 * constructor
 * @param  name Name of the model
 * @param  run  XeRun on which the model is applied
 * @param  requested requested analysis mode: 
 * either NO_ANALYSIS, CUTS_ANALYSIS or PL_ANALYSIS
 */
     XeBackgroundModel(string name, XeRun* run, int requested);

/**
 * No arg constructor for root
 */
     XeBackgroundModel();

 /**
     * Check compatibility
     * @return flag telling that everything is ok
     * @param recomputeIfNeeded Recompute when thing aren't initialized
 */
    bool   checkIt(bool recomputeIfNeeded=true);

     void     printExpected();

 /**
     * return expected background
 */
     double   getExpected();

 /**
     * set value of expected backgound
     * @param exp expected background
 */
     void     setExpected(double exp);

   protected :

     int      requestedAnalysisMode;
     int      modelType;
     double   dayKg;
     double   expected;

     ClassDef(XeBackgroundModel,1)
};


/**
   * Class for neutron background model
*/
class NeutronBackgroundModel : public XeBackgroundModel {

  public:
  
     virtual ~NeutronBackgroundModel();

 /**
 * constructor
 * @param  name Name of the model
 * @param  run  XeRun on which the model is applied
 * @param  requested requested analysis mode: 
 */
     NeutronBackgroundModel(string name, XeRun* run,int requested);

/**
 * No arg constructor for root
 */
     NeutronBackgroundModel();
    

 /**
     * return  estimate in neutrons in Events/Kg/Day/Pe below
     * @param pe  upper equivalent number of photo electrons
 */
     double expected(double peMin,double peMax);

 /**
     * return  estimate in neutrons in Events/Kg/Day/Pe below
     * @param pe  upper equivalent number of photo electrons
 */
     virtual double expectedBelow(double pe)=0;

 /**
     * return a default model for a given run
     * @param  runNumber  The given run number
 */
    static NeutronBackgroundModel* newDefault(XeRun* run);

      ClassDef(NeutronBackgroundModel,1)
 /**
     * print details of background estimation
     * @param level  detail level
 */
     bool printIt(int level);

};


/**
   * Neutron background model Run10
*/
class NeutronBackgroundModelRun10 : public NeutronBackgroundModel {

  public:

     virtual ~NeutronBackgroundModelRun10();

 /**
     * constructor
     * @param  run  XeRun on which the model is applied
 */
     NeutronBackgroundModelRun10(XeRun* run);

/**
 * No arg constructor for root
 */
     NeutronBackgroundModelRun10();

 /**
     * return  Bands for (regular) background
 */
     S1S2Bands* computeBands();

 /**
     * return  estimate in neutrons in Events/Kg/Day/Pe below
     * @param pe  upper equivalent number of photo electrons
 */
     double expectedBelow(double pe);

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);


  protected :

    static double  _integrated[N_PE_NR_BACKGROUND];

      ClassDef(NeutronBackgroundModelRun10,1)

};


/**
   * Neutron background model Run8
*/
class NeutronBackgroundModelRun8 : public NeutronBackgroundModel {

  public:

     virtual ~NeutronBackgroundModelRun8();

 /**
     * constructor
     * @param  run  XeRun on which the model is applied
 */
     NeutronBackgroundModelRun8(XeRun* run);

/**
 * Empty constructor for ROOT
 */
     NeutronBackgroundModelRun8();

 /**
     * return  Bands for (regular) background
 */
     S1S2Bands* computeBands();


 /**
     * return  estimate in neutrons in Events/Kg/Day/Pe
     * @param pe  upper equivalent number of photo electrons
 */
     double expectedBelow(double pe);

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);


      ClassDef(NeutronBackgroundModelRun8,1)

};

/**
   * Class for electron background model
*/
class ElectronBackgroundModel : public XeBackgroundModel {


  public :

     virtual ~ElectronBackgroundModel();

 /**
 *  constructor
 * @param  name Name of the model
 * @param  run  XeRun on which the model is applied
 * @param  requested requested analysis mode: 
 */
     ElectronBackgroundModel(string name, XeRun* run,int requested);

/**
 * No arg constructor for root
 */
     ElectronBackgroundModel();

 /**
     * return  Bands for anomalous background
 */
     virtual  S1S2Bands* computeAnomalousBands()=0;

 /**
     * return a default model for a given run
     * @param  runNumber  The given run number
 */
    static ElectronBackgroundModel* newDefault(XeRun* run);

 /**
     * return expected gaussian background
 */
     double   getExpectedGaussian();

 /**
     * set value of expected gaussian  backgound
     * @param exp expected gaussian background
 */
     void     setExpectedGaussian(double exp);

 /**
     * return expected anomalous background
 */
     double   getExpectedAnomalous();

 /**
     * set value of expected anomalous  backgound
     * @param exp expected anomalous background
 */
     void     setExpectedAnomalous(double exp);

 /**
     * set value of exponential slope in S1 anomalous
     * @param slope expected anomalous background (<=0: uniform)
 */
     void     setAnomalousSlope(double slope);

 /**
     * return S1 slope of anomalous background
 */
     double   getAnomalousSlope();
/**
 * Set ER Gaussian S1 distribution be flat
 * @param flat do we want it flat?
 */
  void setERGaussianFlatInS1(bool flat);

/**
 * Is ER Gaussian S1 distribution flat ?
 */
  bool isERGaussianFlatInS1();

 /**
     * print details of background estimation
     * @param level  detail level
 */
     virtual bool printIt(int level=1);

  protected :


    bool   ERGaussianFlatInS1;    /*!<Do we want Er Gaussian S1 dist. be flat?*/

    double normCo60   ;                         /*!<Normalization of Co60 runs*/
    double meanCo60   ;                        /*!< Average of flattened S2/S1*/
    double sigmaCo60  ;                          /*!< Sigma of flattened S1/S1*/

    double gaussianExpected;                    /*!< Expected gaussian leakage*/
    double anomalousExpected;                  /*!< Expected anomalous leakage*/
    double anomalousSlope;                  /*!< S1 slope of anomalous leakage*/

    ClassDef(ElectronBackgroundModel,1)

};

/**
   * Simplistic Electron background model runs 8 an 10.
   * It uses hard coded values of average and sigma of gaussian backgroud;
   * Then it goes over all bands and computes the anomalous leakage for
   * each of them
*/
class SimplisticERBackground : public ElectronBackgroundModel {

  public :

     virtual ~SimplisticERBackground();

 /**
     *  Regular constructor
     * @param  name Name of the model
     * @param  run  XeRun on which the model is applied
 */
     SimplisticERBackground( XeRun* run);

/**
 * No arg constructor for root
 */
     SimplisticERBackground();

 /**
     * return  Bands for (regular) background
 */
     S1S2Bands* computeBands();

 /**
     * return  Bands for anomalous background
 */
     S1S2Bands* computeAnomalousBands();

 /**
     * Assume leakage is equally spread over the bands
     * @param spread  yes or no?
 */
    void spreadLeakageOverBands(bool spread);

 /**
     * return is Leakage is equally spread over the bands
 */
    bool isLeakageSpreadOverBands();
/**  
 *  print this simplistic background model
 *  @param level print level
 */
    bool printIt(int level=1);         

  protected :
  

    bool   doSpreadLeakageOverBands ;

      ClassDef(SimplisticERBackground,1)

};

/**
   * Simple Electron background model runs 8 an s10
   * It recomputes the gaussian parameters and establishes globally
   * the anomalous, assuming it constant in the flattened (log(s2/s1) variable
*/
class SimpleERBackground : public ElectronBackgroundModel {

  public :

    ~SimpleERBackground();

 /**
     *  Regular constructor
     * @param  run  XeRun on which the model is applied
 */
     SimpleERBackground(XeRun* run);

/**
 * No arg constructor for root
 */
     SimpleERBackground();
 /**
     * return  Bands for (regular) background
 */
    S1S2Bands* computeBands();

/**
 * Set the overall minimum limit in flattened s2/s1 space
 * @param fMin minimum value
 */
    void setFlattenedMin(double fMin=DEFAULT_FLATTENED_MIN);

/**
 * Set the overall maximum limit in flattened s2/s1 space
 * @param fMax maximum value
 */
    void setFlattenedMax(double fMax=DEFAULT_FLATTENED_MAX);

/**
 * Set the overall limits in flattened s2/s1 space, symmetric around zero
 * @param fSymmetric  range will be [-fSymmetric,+fSymmetric]]
 */
    void setFlattenedSymmetric(double fSymmetric);

/**
 * Set the  minimum limit in flattened s2/s1 space for the fit of regular ER
 * @param fMin minimum value
 */
    void setFlattenedFitMin(double fMin=DEFAULT_FLATTENED_FIT_MIN);

/**
 * Set the  maximum limit in flattened s2/s1 space for the fit of regular ER
 * @param fMax maximum value
 */
    void setFlattenedFitMax(double fMax=DEFAULT_FLATTENED_FIT_MAX);

/**
 * Set symmetric limits in flattened s2/s1 space for the fit of regular ER
 * @param fSymmetric  range will be [-fSymmetric,+fSymmetric]]
 */
    void setFlattenedFitSymmetric(double fSymmetric);

/**
 * Set the  minimum limit in flattened s2/s1 space for the anomalous leakage
 * @param fMin minimum value
 */
    void setFlattenedAnomalousMin(double fMin=DEFAULT_FLATTENED_ANOMALOUS_MIN);

/**
 * Set the  maximum limit in flattened s2/s1 space for the anomalous leakage
 * @param fMax maximum value
 */
    void setFlattenedAnomalousMax(double fMax=DEFAULT_FLATTENED_ANOMALOUS_MAX);

/**
 * Set symmetric limits in flattened s2/s1 space for the anomalous leakage
 * @param fSymmetric  range will be [-fSymmetric,+fSymmetric]]
 */
    void setFlattenedAnomalousSymmetric(double fSymmetric);

/**
 * Computes the exponetial factor "anomalousSlope" 
 * of anomalous bkg. and writes it into anomalousSlope member.
 * Based on Hagar's implementation, wiki: 
 * https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:run10_profile_likelihood
 */
    void computeAnomalousSlope();

 /**
     * return  Bands for anomalous background
 */
    S1S2Bands* computeAnomalousBands();
/**  
 *  print this simple background model
 *  @param level print level
 */
    bool printIt(int level=1);         

  protected :


    double flattenedMin          ;
    double flattenedMax          ;  
    double flattenedFitMin       ; 
    double flattenedFitMax       ;  
    double flattenedAnomalousMin ; 
    double flattenedAnomalousMax ;  

    ClassDef(SimpleERBackground,1)

};

/**
   * Class for ER background as published in PRL
   * Everything is tabulated 
*/
class PublishedElectronBackgroundRun10 : public ElectronBackgroundModel {

  public:

/**
 *  constructor
 *  @param  run  XeRun on which the model is applied
 */
     PublishedElectronBackgroundRun10(XeRun* run);

/**
 * No arg constructor for root
 */
     PublishedElectronBackgroundRun10();

    ~PublishedElectronBackgroundRun10();

 /**
    * return  Bands for anomalous background
 */
     S1S2Bands* computeBands();

 /**
 *  return  Bands for (regular) background
 */
     S1S2Bands* computeAnomalousBands();

/**
 * Set tabulated background
 * @param backgrounds a pointer to the array
 */
     void setBackground(double *backgrouds);

/**
 * Set tabulated anomalous background
 * @param backgrounds a pointer to the array
 */
     void setAnomalousBackground(double *backgrouds);
     
/**  
 *  print this background model
 *  @param level print level
 */
     bool printIt(int level=1);   

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);


  protected :
  
    vector<double> bkgs;
    vector<double> anomalousBkgs;

    ClassDef(PublishedElectronBackgroundRun10,1)

};


/**
   * Class for ER background as published in PRL for run 8
   * Everything is tabulated
   * Based on code from Ofer: ofer_bands.c and compute_single_band_content.c
   * Not yet verified with what reported in https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat3#final_background_prediction_x48kg0 
*/
class PublishedElectronBackgroundRun8 : public ElectronBackgroundModel {

  public:

/**
 *  constructor
 *  @param  run  XeRun on which the model is applied
 */
     PublishedElectronBackgroundRun8(XeRun* run);

/**
 * No arg constructor for root
 */
     PublishedElectronBackgroundRun8();

    ~PublishedElectronBackgroundRun8();

 /**
    * compute  Bands for anomalous background and returns them
 */
     S1S2Bands* computeBands();

 /**
 *  computes  Bands for (regular) background and returns them
 */
     S1S2Bands* computeAnomalousBands();

 /**
 *  return the Bands for (regular) background (runs the computation if needed)
 */
     S1S2Bands* getGausBands();

 /**
 *  return the Bands for Anomalous background (runs the computation if needed)
 */
     S1S2Bands* getLeakBands();
          
/**  
 *  print this background model
 *  @param level print level
 */
     bool printIt(int level=1);   

/** 
 * Implementation of the virtual method to define run compatibility
 * @param runNumber run number to be checked
 * @return ok or not
 */
  bool isRunCompatible(int runNumber);


  protected :
  

    float leak_RangeFraction; //Factor to correct from the 3-35pe range where the leakage is computed to the actual 4-30pe range
    TF1 *fsigmaco60;          //Sigma of the Co60 band as a function of S1
    int npco60;		      //Number of Co60 points passing all cuts between 3. and 35. pe in S1	

    S1S2Bands * gausBands;      // bands for gaussian component of ER bkg
    S1S2Bands * leakBands;      // bands for leakage component of ER bkg

    bool 	isLeakComputed;
    bool	isGausComputed;

    ClassDef(PublishedElectronBackgroundRun8,1)

};

/**
   * ModelImporter Class which holds general methods to import bkg and signal models from a TH2F and push it trough S1S2Bands.

   * Usage example:
     ModelImporter gauss_bkg(); 					// empty model
     gauss_bkg.setFile(f);      					// specify the file
     gauss_bkg.import("name_of_TH2F");  				// import the TH2F
     S1S2Bands *bandsFromTH2F = gauss_bkg.fillBands(referenceBands);    // get the S1S2Bands

   * Remarks:
     TH2F should have X-axis = S1 and Y-axis = log10(S2/S1) not flattened.
     The binning on X-axis has to be of 1 PE with center in X.5 were X is an integer. Example TH1D::[50, 0., 50.] (50 bins from 0. to 50.).
     The binning on Y-axis has to be smaller than the size of the smallest band, it will return empty bands otherwise.
     There is no limit on X and Y range, it will be cutted as defined in the reference bands.
     Keep the same TH2F binning for signal and bkg. 
*/
class ModelImporter
{

   public: 

/**
 * Costructor. Parameters:
 * @param file = TFile from which import the TH2F
 * Naming convention for TH2F:  
*/
	ModelImporter(TFile *file = NULL, TString name = "");

	~ModelImporter();	


/**
  * Set name of the histo you want to import and imports it if found in current file.
*/
	void setHistoName(TString name);	

/**
 * Sets the file
*/
	void setFile(TFile *file);

/**
 * Returns a copy of the imported TH2F
*/
	TH2F getModel();

/**
  * Returns the name of the current histogram.
*/
	TString getHistoName(){return histo_name;};

/**
 * Returns a copy of "reference_band" filled with the content of imported TH2F.
 * See remarks of this class for more details.
 * @param isLog has to be true if the histo is in log(s2/s1) VS s1 space, false if in s2 VS S1.
*/
	S1S2Bands *fillBands(S1S2Bands *reference_band, bool isLog = true); //should NOT return a pointer to local var FIXME

/**
 * Rebins the histogram's X axis.
 * @param Rebin rebin value as in TH1D::Rebin(rebin_val)
*/
    void rebinX(int Rebin);

/**
  * Scale, this is a non automatic feature, you'll have to apply yourself if needed
*/
   void scale(double norm) {if(model) model->Scale(norm);};

   private:

/**
 * Imports the TH2F you want to push into bands, the name has to be set.
*/
	void import();

	TString histo_name;	/*!<Name of histogram you want to import*/
	TH2F	*model;		/*!<Pointer to the stored histo */
	TFile   *f;		/*!<pointer to file from which one takes the histos*/

};



/**
   * Signal Model virtual class. It handles the generation of signal S1S2Bands for a fixed mass 
   * and for different value of Leff and Qy, it also generates reference bands.
*/
class SignalModel : public XeObject{
	
	public :

	  SignalModel(XeRun *run);
	  ~SignalModel(){};

	  /**
	    * Clones reference bands and fill it with computed signal
	    * @param "reference" bands (NOTE: parameter is passed by value)
	    * @param "T_Leff"  Leff T value
	    * @param "T_Qy"    Qy (charge yield) T value	    
	  */
	  virtual S1S2Bands computeSignalBands(double T_Leff, double T_Qy);


	  /**
	    * Returns the cross section multiplier, see Xsecion_multiplier for details.
	  */	
  	  double getXsecMultiplier();

 	 /**
	   * Set the scale factor to apply to the signal sample, Xsecion_multiplier.
	 */
	  void setXsecMultiplier(double factor); 

	  void setDefaultNorm(double factor)     {Default_XsectionNorm = factor;};
	
  	  double getDefaultNorm() 		{return Default_XsectionNorm;};

	  void normalizeToOne()                 {normOne = true;};

	  bool isNormOne()                      {return normOne;};
	
	  /**
	    * Returns the Xsection relative to the signal strenght parameter in the limits.
	      The limits are set in terms of "sigma" which goes from [0:100], to get the meaningful
	      cross section relative to that limit "sigma" has to be multiplied by this parameter.
	  */
	  double getRelativeXsec()		{return Default_XsectionNorm * Xsecion_multiplier;};

	  /**
	    * Allows to re-initialize the class with a new mass input.
	    * @param mass of the WIMP.
	  */	
 	  virtual bool setMass(double mass);


	  /**
	    * Compute the reference bands from signal template with nominal Leff and Qy.
	    * Default is bands from AmBe data points, but NOTE that this is a virtual method 
	    * and it is overwritten by daughters.
	    * Returns empty S1S2Bands.
	  */
	  virtual S1S2Bands computeReferenceBands();

	  /**
	    * Current mass in GeV.
	  */
	  double getCurrentMass(){return mass;};

	  /**
	    * Setting the printlevel for info.
	    * @param level if > 0 print additional info for debug, default is 0.
	  */
	  void setPrintLevel(int level) {printlevel = level;};
	
	protected:

	  XeRun   *currentRun; 		/*!<Pointer to the Run */
	  double  mass;			/*!<Current WIMP mass Value */
	  double  Xsection;		/*!<Final cross section corresponding to the actual number of events*/					
	  double  Default_XsectionNorm;     /*!<cross section with which the signal events are computed by default, example: this is the Xsec relative to the number of event in the histogram.*/
	  double  Xsecion_multiplier;     /*!<The limit on cross section will span several order of magnitude, Minuit "doesn't like" it, so the cross section for each mass is renormalized in order to keep the minuit parameter always within [0,100]. */ 
	  int     printlevel;

	  bool    normOne;




};



class SignalModel1D : virtual public SignalModel {

    public:

	~SignalModel1D(){};
	SignalModel1D(XeRun *run);
	
	void setSignalName(TString nameHisto);
	void setLeffSysNames(TString nameHisto_psigma, TString nameHisto_msigma);

 	double getLeffmultiplier(double t_val, int bin_number) ;

	void   normalizeToOne();

	TH1F *getDefaultSignal(){return Defaultsignal;};
	TH1F *getPsigmaSignal(){return PsigmaSignal;};
	TH1F *getMsigmaSignal(){return MsigmaSignal;};

	TH1F *Defaultsignal;

	TH1F *PsigmaSignal;  // plus 1 sigma signal FIXME:: this is supposed to be in a systematic uncertainty handler, but no time now....

	TH1F *MsigmaSignal;  // minus 1 sigma Leff

};


/**
  * Class that produces and fills bands by IMPORTING MC simulated signal histograms.
  * Many thanks to Hagar Landsman for her precious help building this class.
*/

class MCSignalModel : public SignalModel {

	public :

	/**
	  * Costructor  default rebinning is 10
        */  
	  MCSignalModel(XeRun *run, int Rebin =10, TString suffix="_");

	  ~MCSignalModel(){};

	/** 
	  * Allows to change wimp mass, checks whether the new mass input is 
	  * compatible with the one defined in the current run.
	  * @param mass_input wimp mass 
	*/
	 bool setMass(double mass_input);
	
	 
	/**
	  * Implementation of the virtual class from SignalModel,
	  * Retrieve signal bands from Histogram.
	  * @param T_Leff Leff t-Value 
	  * @param T_Qy Qy t-Value 
	*/
	 S1S2Bands computeSignalBands(double T_Leff, double T_Qy);

	/**
	  * Compute the reference bands from signal template 2D histogram with nominal  Leff and Qy.
	*/
	S1S2Bands computeReferenceBands();  

	/**
	  * Returns the loaded Signal Histo for the current wimp mass, Leff and Qy = 0.
	*/
	TH2F getHisto(); 

	private :
	/**
	  * return the name of the histogram to be loaded, based on MC signal naming convention.
	*/	
	  TString computeName(double Leff = 0., double Qy = 0. , double LCE = 1.,TString suffix="_");
	  
	  ModelImporter signalInput;  /*<! class that load, holds and push into bands 2D histograms. */
	  
	  int rebin;
	  TString suffisso; 

};


/**
   * The basis for an analysis : a Run with its (optional) data .
   * Contains S1S2Data for each of  AM_BE_DATA , E_GAMMA_DATA  , DM_DATA   
   * , AM_BE_CUT_DATA , E_GAMMA_CUT_DATA    , DM_CUT_DATA .
   * Contains S1S2Bands for the same data sets.
   * S1S2Bands of (nominal) background are created; they are called:
   * ER_BACKGROUND, ER_LEAKAGE, NR_BACKGROUND and ALL_BACKGROUNDS
*/
class XeRun : public integrable , public S1S2Object {

  public:
     

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     virtual      ~XeRun();


/**
 * No arg constructor for root
 */
     XeRun();

 /**
     * @param  run   runNumber
     * @param  datatype    data type (REAL_DATA or NO_DATA)
     * @param  file  name of root file containing graphs
     * @param  dmGraph     TGraph containing DM data
     * @param  erGraph     TGraph containing Electronic recoil data
     * @param  nrGraph     TGraph containing neutron recoil data
 */
  XeRun( int run              ,int datatype           ,string file="", string sFile="" 
     , string dmGraph="grNONE",string erGraph="grCo60",string nrGraph="grAmBe");

 /**
     * Set number of bands for analysis
     * @param   nb number of bands
 */
    static void    setNBands(int nb=DEFAULT_N_BANDS);
 /**
     * Set number of slices for analysis
     * @param   ns  Number of slices
 */
    static void    setNSlices(int ns=DEFAULT_N_SLICES);

/**
 * return short type type name, using value of cross section for graphs
 * @param  type Data type
 */
    string getShortTypeName(int type);



/**
  * return a 1D Graph from a TGraph2D
*/

   TGraph from2dto1dGraph(TGraph2D gr2d);

/**
 * Does the current run have bands
 */
   bool withBands();

/**
 * Is the current run without bands?
 */
   bool withoutBands();

 /**
     * Dump "candidates" (i.e. events in lowest bands)
     * @param   Maximum band number (number from 0)
     * @return went ok?
 */
    bool printCandidates(int maxBand);

 /**
     * Dump events of a give data source looking like "candidates"
     * (i.e. events in lowest bands)
     * @param   dt datatype`
     * @param   Maximum band number (number from 0)
     * @return went ok?
 */
    bool printEvents(int dt, int maxBand);

/**
 * Check whether a band index is correct and prints a warning if not
 * @param banf band index, from 0 to nBands-1, or "ALL"
 */
  bool checkBand(int band);


  MCSignalModel *getSignalHandler() 		{ return signalHandler; }; 

  double getSignalDefaultNorm()     		{ return signalHandler->getDefaultNorm(); };

  void   setSignalDefaultNorm(double norm)      { signalHandler->setDefaultNorm(norm); }

  double getSignalMultiplier()	    		{ return signalHandler->getXsecMultiplier(); };

  void   setSignalMultiplier(double factor)	{ signalHandler->setXsecMultiplier(factor); };


 /**
     * return data type
     * @return   NO_DATA, REAL_DATA, SIMULATED_DATA
 */
    int  getDataType();

/**
 * get the default file directory for a given run
 * @param rn run number
 */
    static string  getDefaultDirectory(int rn);
 
/**
 * get the default file name for a given run
 * @param rn run number
 */
    static string  getDefaultFileName(int rn);
   
/**
 * get the data type name for a given data type
 * @param dt data type (NO_DATA, REAL_DATA or SIMULATED_DATA)
 */
    static string  getDataTypeName(int dt);
   
/**
 * get the data type name for the current name
 */
    string  getDataTypeName();

/**
 * Method overriding the default one for XeObject.
 * It adds the data type
 */
    string  getName();
    
/**
  * Get number of bands defined for the current analysis.
*/
    int getNbands(){return nBands;};

/**
  * Get number of slices defined for the current analysis.
*/
    int getNslices(){return nSlices;};
               
/**
 * Draw it, being and S1S2Object object
 */
    void draw();
                   
/**
 * Draw the cuts
 */
    void drawCuts();

 /**
     * return Pointer to a default LEff class
 */
    LEff*  newDefaultLEff();

 /**
     * return Pointer to a default Set Of Selection Cuts 
 */
    XeSetOfSelectionCuts* newDefaultXeSetOfSelectionCuts();

/**
     * return Pointer to a default Set Of Dark Matter Cuts 
 */
    XeSetOfCuts*  newDefaultXeSetOfDarkMatterCuts();


 /**
     * return Pointer to a default Neutron Background Model 
 */
    NeutronBackgroundModel*  newDefaultNeutronBackgroundModel();

 /**
     * return pointer to  Neutron Background Model in use
 */
    NeutronBackgroundModel* getNeutronBackgroundModel();

 /**
     * return pointer to a default Electron Background Model 
 */
    ElectronBackgroundModel* newDefaultElectronBackgroundModel();

 /**
     * return pointer to  Electron Background Model in use
 */
    ElectronBackgroundModel* getElectronBackgroundModel();

 /**
     * Set selection cuts 
     * @param cuts  selection cuts to be set
 */
    void setSelectionCuts(XeSetOfSelectionCuts* cuts);

 /**
     * Set analysis (dark matter) cuts 
     * @param cuts  analysis (dark matter) cuts to be set
 */
    void setDarkMatterCuts(XeSetOfCuts* cuts);

 /**
     * Set model for computing electronic background (leakage)
     * @param model  The model to be used
 */
    void setElectronBackgroundModel(ElectronBackgroundModel* model);

 /**
     * Use the published background model
 */
    void setPublishedElectronBackgroundModel();

 /**
     * Set model for computing neutron background
     * @param model  The model to be used
 */
    void setNeutronBackgroundModel(NeutronBackgroundModel* model);

  
    void setMCSignalModel( MCSignalModel *signal) { signalHandler = signal;};



/** 
 * Set the interaction
 * @param inter a pointer the Interaction
 */
    void setInteraction(Interaction* inter);

/** 
 * Set the Galaxy model
 * @param gal pointer to the GalaxyModel
 */
    void setGalaxyModel(GalaxyModel* gal);

/** 
 * Set the WIMP Particle
 * @param wimp pointer to the Wimp
 */
    void setWimp(Wimp *wimp);

/** 
 * Set the WIMP mass
 * @param mass the requested mass
 */
    void setWimpMass(double mass);

/** 
 * Set the wimp-nucleon cross section
 * @param sigma the cross section
 */
    void setSigmaNucleon(double sigma);

/** 
 * Set the wimp-nucleon cross section for reports, graphs, tables...   
 * @param sigma the cross section
 * @param   unit  either SIGMA_UNIT or EVENT_UNIT
 * 
 */
    void setSigmaForReports(double sigma, int unit=SIGMA_UNIT);

/** 
 * Fet the wimp-nucleon cross section for reports, graphs, tables... 
 */
    double getSigmaForReports();


/**
 * return expected number of signal events for 1 cm2 cross section
 * @param lt L-Eff t-value
 */
    double nSignalPerCm2(double lt=0.);            

/** 
 * Set the target
 * @param target pointer to a Target XeObject
 */
    void setTarget(Target *target);

/** 
 * Set the  LEffective model
 * @param leff pointer to a LEff XeObject
 */
    void setLEff(LEff* leff);

/** 
 * Set the  LEffective "t-value"
 * @param tValue How many sigmas above or below central L-effective
 */
    void setLEffTValue(double tValue);

 /**
     * set hard cut of lowest energy for which LEff>0
     * @param e energy in KeV
 */
    void setLEffErMin(double e);

/** 
 * Set the  light yield
 * @param lightYield  the value
 */
    void setLightYield(double lightYield);

/**
 * Set the mode when s1 slices are equally filled
 * @param equal  true if yes, general case is no
*/
    void setEquallyFilledS1Slices(bool equal=false);

 /**
     * set the highest S2OverS1 limit in band 0
     * @param S2OverS1Max the maximal value
 */
    void setHighestS2OverS1(double S2OverS1Max);

 /**
     * set the lowest S2OverS1 limit in band 0
     * @param S2OverS1Min the minimal value
 */
    void setLowestS2OverS1(double S2OverS1Min);

 /**
     * set the lowest s1 value for band analysis
     * @param s1Min the minimal value
 */
    void setAnalysisS1Min(double s1Min);

 /**
     * set the highest s1 value for band analysis
     * @param s1Max the maximal value
 */
    void setAnalysisS1Max(double s1Max);

/** 
 * Set the  "t-value" for a selection cut
 * @param  i cut index starting at zero
 * @param tValue How many sigmas above or below central cut
 */
    void setSelectionCutTValue(int i, double t);

 /**
     * Fill the expected number of signal events for the sequence of selection cuts
     * @param  remaining  vector of remaing events after each selection cut
     * @param  LEfft  LEffective t-parameter
 */
    void fillSelectionCutsBreakdown(vector<double> &remaining,double LEfft=0.);

 /**
     * Dump of the expected number of signal events for the sequence of selection cuts
     * @param  mass  WIMP mass (0.: keep unchanged)
     * @param  sigma  cross section (0.: keep unchanged)
     * @param  LEfft  LEffective t-parameter
 */
    void printSelectionCutsBreakdown(double mass=0.,double sigma=0.
                                   ,double LEfft=0.);

 /**
     * simulate a run (i.e. fill DM_CUT_DATA data sets and bands)
     * @param      sigma  cross section being simulated
 */
    void simulate(double sigma=0.);

 /**
     *             fill simulated data, and corresponding band/slices
     * @return     true if s1 inside, false if not (this is an error condition)
     * @param      band  band number to be filled
     * @param      s1    simulated S1
 */
    bool fillSimulatedDataAndBands(int band,double s1);   

/**
 * Fill all data bands, after they have been defined
 */
    void fillDataBands();   

/**
 * Generate a DM data set according to the bkg model, fluctuating just stat. uncertainties.
   To retrieve the generated S1S2Band use the method XeRun::getS1S2Bands(DM_SIMULATED_DATA)
   @param seed, the seed with which you'd like to generate toys (TRandom3 is used)
   @param sigma, is the signal strenght that you want to inject.
 */
   void generateData(double seed = 1965, double sigma =0.);

/**
 * Generate a ASIMOV data set according to the bkg model + signal with signal strenght sigma,
   To retrieve the generated S1S2Band use the method XeRun::getS1S2Bands(ASIMOV_DATA)
 */
   void generateAsimovData(double sigma);

/** 
 * return run number
 */
    int  getNumber(); 

/** 
 * return number of selection cuts
 */
    int  getNSelectionCuts();

/**
 * This run was defined with data, i.e. it's not a mock-up
 */
    bool withData();

/**
 * This run was defined without data, i.e. it's a mock-up
 */
    bool withoutData();

 /**
     * Check the run attributes and recompute what is needed
     * @return flag telling that everything is ok
     * @param recomputeIfNeeded Recompute when thing aren't initialized
 */
    bool   checkIt(bool recomputeIfNeeded=true);
    
 /**
     * Force all run checks and band recomputing
     * @return flag telling that everything is ok
 */
    bool   recomputeEverything();

/**
 * Return current "t-value" for L-effective
 */
    double getLEffTValue();

 /**
     * get hard cut of lowest energy for which LEff>0
 */
    double getLEffErMin();

 /**
     * Return maximum Recoil energy for given mass
 */
    double ErMax();

/**
 * Return expected number of photons
 * @param Er Recoil Energy
 */
    double ErToNPhotons(double Er);
  
 /** 
 * Return Recoil energy, given expect number of photons
 * @param nP expected number of photons
 */
   double  NPhotonsToEr(double nP);

 /**
     * Compute event rates in events/KeV
     * @return  rate in events/KeV
     * @param   Er Recoil energy
 */
    double dRate(double Er);

 /**
     * get typical lower cut on recooil energy corresponding to cuts
 */
    double getTypicalErMin();

/**
 * return Wimp Mass
 */
    double getWimpMass();

/**
 * Return exposure (in days)
 */
    double getExposure();

/**
 * Return fiducial mass
 */
    double getFiducialMass();

/**
 *  Return wimp-nucleon cross section
 */
    double getSigmaNucleon();
   
/**
 * Return resolution of PMT
 */
    double getSigmaPMT();
   
/**
 * Return unsmeared S1 min for the relevant cuts
 */
    double getSelectionUS1Min();
   
/**
 * Return smeared S1 min for the relevant cuts
 */
    double getSelectionS1Min();
   
/**
 * Return smeared S1 max for the relevant cuts
 */
    double getSelectionS1Max();

/**
 * Get overall S1 min for the s1,s2 analysis
 */
    double getAnalysisS1Min();

/**
 * Get overall S1 max for the s1,s2 analysis
 */
    double getAnalysisS1Max();

/**
 * Get lowest overall s2/s1 for the s1,s2 analysis
*/
    double getLowestS2OverS1();

/**
 * Get highest overall s2/s1 for the s1,s2 analysis
*/
    double getHighestS2OverS1();

/**
 * Return light yield
 */
    double getLightYield();

/**
 * Return default E/gamma normalisation (i.e. expected number of bkg events)
 */
    double defaultEGammaNormalizationToEvents();

/**
 * Return default Neutron/gamma normalisation (i.e. expected number of bkg events)
 */
    double defaultAmBeNormalizationToEvents();

/**
 * Return default DarkMatter normalisation (i.e. 1)
 */
    double defaultDarkMatterNormalization();

/** 
 * Internal method: were thu cuts applied?
 */ 
   void    checkTheCuts();

 /**
     * print details of band contents for variable LEff
     * @param    leffT  t-value describing LEff
 */
    void printBandContent(double leffT=0.);

/**
 * Print the CLs limits while cumulating bands one to one
 * @param nb  number of bands to consider
 */
    void printBandCLsCounting(int nb);

 /**
     * Print a summary of the data bands
 */
    void printDataSummary();

/**
 * Return overall S1 max
 */
    double maxS1();

/**
 * Return overall S2 max
 */
    double maxS2();

/**
 * Return overall S1 min
 */
    double minS1();

/**
 * Return overall S2 min
 */
    double minS2();

/**
 * Return maximum 'X' value, as requested by being an S1S2Object object
 */
    double maxX();

/**
 * Return maximum 'y' value, as requested by being an S1S2Object object
 */
    double maxY();

/**
 * Return minimum 'X' value, as requested by being an S1S2Object object
 */
    double minX();

/**
 * Return minimum 'Y' value, as requested by being an S1S2Object object
 */
    double minY();
    
/**
 * return maximum S1 for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double maxS1(int type);
    
/**
 * return maximum S2 for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double maxS2(int type);
    
/**
 * return minimum S1 for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double minS1(int type);
    
/**
 * return minimum S2 for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double minS2(int type);
    
/**
 * return maximum 'X' for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double maxX(int type);
    
/**
 * return maximum 'Y' for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double maxY(int type);

/**
 * return minimum 'X' for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double minX(int type);

/**
 * return minimum 'Y' for a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    double minY(int type);

 /**
     * events rate for integration
     * @return The value for integration
     * @param  what  Sor far implemented RATE  -> dRate(x)
     * @param x      the point for the integrable is computed
 */
    double getValue(int what, double x);

 /**
     * compute flattened log(s1/s2) 
     * @return  flattened log(S1/S2)
     * @param   S1  S1 value
     * @param   S2  S2 value
 */
    double flatten(double S1, double S2);

/**
 * Return Se used for Light yield computation
 */
    double getSe();

/**
 * Return Sr used for Light yield computation
 */
    double getSr();

/**
 * Return L-effective model user for Light yield computation
 * @return pointer to LEff XeObject
 */
    LEff   *getLEff();
/**
 * Return Qy model user for Light yield computation
 * @return pointer to LEff XeObject
 */
    Qy   *getQy() {return qy;};

/**
 * Return target
 * @return pointer to Target XeObject
 */
    Target *getTarget();

/**
 * Return WIMP
 * @return pointer to Wimp XeObject
 */
    Wimp   *getWimp();
   
 /**
 * Return the pValue currently is use
 */
   PValue* getPValue();

 /**
 * Set the pValue currently is use
 * @param pV current PValue
 */
   void setPValue(PValue* pV);
 
/**
 * Return Galaxy model
 * @return pointer to Galaxy Model XeObject
 */
    GalaxyModel   *getGalaxyModel();

/**
 * Return object which flattens log(s2/s1)
 * @return pointer to RunFlattener XeObject
 */
    RunFlattener  *getFlattener();

/**
 * Return interaction
 * @return pointer to Interaction XeObject
 */
    Interaction   *getInteraction();

 /**
     * Get a tree corresponding to a given data file
     * @return pointer to a newly created TTree object
     * @param  type data type
 */
    TTree  *newTree(int type);

 /**
 * Get data set
 * @param type either AM_BE_DATA,E_GAMMA_DATA, or DM_CUT_DATA
 */
    S1S2Data *getS1S2Data(int type);

 /**
 * Get reference bands.
 * @return pointer to a S1S2Bands XeObject
 */
    S1S2Bands *getReferenceBands();

 /**
 * Get bands of  a given data type
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
    S1S2Bands *getS1S2Bands(int type);

 /**
     * create  a 2d histogram whose content is follows slices and bands
     * @return  pointer to the newly  created  TH2F
     * @param type data type
     * @param plot Do we want to plot the content? can be NONE,LINEAR,LOG
 */
    TH2F* new2DHistogram(int type, int plot=NONE);


/**
 * Return the dark matter cuts
 * @return pointer to a XeSetOfCuts XeObject
 */
    XeSetOfCuts *getDarkMatterCuts();

/**
 * Return all the cuts
 * @return pointer to a XeSetOfCuts XeObject
 */
    XeSetOfCuts *getAllCuts();

/**
 * Return a selection cut, given its sequence number
 * @param i index running from 0
 * @return a pointer to a SelectionCut XeObject
 */
    SelectionCut  *getSelectionCutBySequence(int i);

/**
 * Return all the selection cuts
 * @return pointer to a XeSetOfSelectionCuts XeObject
 */
    XeSetOfSelectionCuts *getSelectionCuts();

 /**
     *  Set SignalERecoil spectrum by brute force for (mainly for debugging)
     *  @param      Er  Spectrum, for 0 to N_ERS_POINTS*ER_STEP;
     *              If NULL, assumes that the spectrum is already filled in
 */
    void forceSignalERecoilSpectrum(double* Er=NULL);

 /**
 * Get the target in use for the run
 * @param a mass number 129/131/... or NATURAL,DEPLETED , XE100_MIXTURE
 */
    static Target* getTarget(int a);
 

 /**
     * return  the expected  number of events pass all cuts
     * @param leffT   LEff t-value
 */
    double expectedSignal(double leffT);

/**
 * Return theoretical upper poission sigma limit, assuming no background
 * @param cl confidence level
 * @param afterCuts evaluated after the cuts (alternative: before)
 */
    double getTheoreticalUpperSigma(double cl=DEFAULT_CL
                                           ,bool afterCuts=true);
/**
 * Return expected number of neutrons between two PE
 * @param peMin Lower bound
 * @param peMax Upper bound
 */
    double expectedNumberOfNeutrons(double peMin, double peMax);


 /**
     * Compute the expected number of signal events in an ERecoil window
     * @param E1 lower bound of the integral
     * @param E2 upper bound of the integral
 */
    double expectedSignalInErWindow(double E1=-1., double E2=-1.);

// =======================================================================
// ======== Methods to build bands and tabulate signal/backgrounds
// =======================================================================

/**
    * Tabulate S1S2Bands signal from histo. This is preliminary and has to be harmonized FIXME 
*/
    void tabulateSignalLeffQy();
/**
   * return the index corresponding to that Leff or Qy in the table. This is preliminary and has to be harmonized FIXME
*/
    unsigned int getTabulatedIndex(double Leff);

/**
  * returns interpolated bands between two values of Leff and Qy.
*/
   S1S2Bands interpolateLeffQy(double Leff, double Qy);
 
 /**
     * compute background bands from data, and estimate background
 */
    void computeBackgrounds();

 /**
     * compute the signal bands for a given LEff t-value
     * @param   LEFFt   t-value of LEff
 */
    void computeSignalBands(double LEFFt);

 /**
     * compute exoeceted Er and S1 spectra for signal, fill bands, etc...
 */
    bool computeSignal();

 /**
     *  Tabulate SignalERecoil spectrum, in the usual case it's not forced
 */
    void tabulateSignalERecoilSpectrum();


 /**
     *  Tabulate S1 spectrum
 */
    void tabulateSignalS1Spectrum();

 /**
     *  Tabulate S1 spectrum for one t value
     * @param  tIndex  Index (from 0 to 2*SIDE_LEFF_TABULATED)
 */
    void tabulateSignalS1Spectrum(int tIndex);

/**
    * Return the tabulated band close to the (always underestimating) t-values requested.
*/    
    S1S2Bands *getSignalBand(double t_Leff, double t_Qy);


// =======================================================================
// ========= Methods to access Band information ==============
// =======================================================================

/**
 * Return   content of a given band
 * @param band  which band  (ALL means the all together= the sum)
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
   double getBandContent(int type, int band=ALL);

/**
 * Return total content
 * @param type can be AM_BE_DATA,E_GAMMA_DATA,DM_DATA,ER_BACKGROUND,...
 */
   double getTotalContent(int type);


// =======================================================================
// ========= Methods to access Signal ERecoil information ================
// =======================================================================

  /**
     * Print the spectrum of E-Recoil
 */
    void printSignalERecoilSpectrum();
  
/**
 * Return the spectrum of Energy recoil, for current interaction, mass ,...
 */
    double *getSignalERecoilSpectrum();

/**
 * Return the distrbution of Energy recoil, for current interaction, mass ,...
 */
    double *getSignalERecoilDistribution();

 /**
     * get the spectrum of Signal ERecoil
     * @return  a pointer to th existing XeSpectrum
 */
    XeSpectrum* getSignalERecoilXeSpectrum();

 /**
     * get the distribution of Signal ERecoil
     * @return  a pointer to the existing TabultatedDist
 */
    TabulatedDist* getSignalERecoilXeDist();

 /**
     * Create a graph of Signal Energy Recoil
     * @return  Pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfSignalERecoilSpectrum(int plot=NONE);

 /**
     * Create a graph of Signal Energy Recoil, normalized to 1
     * @return  Pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfSignalERecoilDistribution(int plot=NONE);


 /**
     * Create an histogram of Energy Recoil
     * @return  Pointer to the newly created Xegraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F* newHistogramOfSignalERecoilSpectrum(int plot=NONE);

 /**
     * Create an histogram of Energy Recoil, normalized to 1
     * @return  Pointer to the newly created histogram
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F* newHistogramOfSignalERecoilDisribution(int plot=NONE);

 /**
     * Create a Multigraph of Energy Recoil
     * @return  Pointer to the Multigraph
     * @param Mr  Mass Range (XeRange) if NULL, use default one
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph* newMultiGraphOfSignalERecoilSpectrum(XeRange *mr=NULL
             ,int plot=NONE);

 /**
     * Create a Multigraph of Energy Recoil, normalized to 1
     * @return  Pointer to the Multigraph
     * @param Mr  Mass Range (XeRange) if NULL, use default one
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph* newMultiGraphOfSignalERecoilDistribution(XeRange *mr=NULL
         , int plot=NONE);

 /**
     * Create a 2D histogram of Er signal spectrum, after cuts
     * @param Mr  Mass Range (XeRange) if NULL, use default one
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH2F* new2DHistogramOfSignalERecoilSpectrum(XeRange *mr=NULL,int plot=NONE);
 /**
     * Create a 2D histogram of Energy Recoil
     * return Pointer to a new TH2F
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH2F *new2DHistogramOfSignalERecoilDistribution(XeRange *mr=NULL
   , int plot=NONE);



// =======================================================================
// ======== Methods to create graphs of physical quantities ==============
// =======================================================================

 /**
 * Create a graph of LEffective for the central value 
 * @return  Pointer to the newly created XeGraph
 * @param Er  Energy Range (ErRange) if NULL, use default one
 * @param t   LEff T-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfLEff(ErRange *er=NULL, double t=0., int plot=NONE);

 /**
 * Create a graph of LEffective for the current value 
 * @return  Pointer to the newly created XeGraph
 * @param Er  Energy Range (ErRange) if NULL, use default one
 * @param t   LEff T-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
     XeGraph* newGraphOfCurrentLEff(ErRange *er=NULL, int plot=NONE);

 /**
 * Create  a graph of the acceptance of the selection
 * @param   srange  S1 range (NULL : default one)
 * @param   smear  accpetance of smeared (true) or unsmeared (false)?
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfSelectionCutsAcceptance(S1Range* srange=NULL
                                               ,bool smear=true, int plot=NONE);

 /**
 * Create  a graph of the number of expected events in Er Window, before cuts
 * @param   mr  Mass range (NULL : default one)
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfExpectedSignalInErWindow(XeRange *mr=NULL,int plot=NONE);

 /**
 * Create  a graph of the number of signal events passing all cuts
 * @param   mr  Mass range (NULL : default one)
 * @param   tl  LEff t-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
 XeGraph* newGraphOfExpectedSignal(XeRange *mr=NULL,double tl=0.,int plot=NONE);

 /**
 * Create  a Multigraph of the number of signal events passing all cuts when LEff varies
 * @param   mr  Mass range (NULL : default one)
 * @param   tRange Tvalue Range (NULL= default one)
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph* newMultiGraphOfExpectedSignalByLEff(XeRange *mr=NULL
                                  , tValueRange* tRange=NULL, int plot=NONE);

/**
 * Create  a Multigraph of the number of signal events passing all cuts when VEsc varies
 * @param   mr  Mass range (NULL : default one)
 * @param   vRange VEscape Range (NULL= default one)
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph* newMultiGraphOfExpectedSignalByVEsc(XeRange *mr=NULL
                        , VEscRange* vRange=NULL, int plot=NONE);

 /**
 * Create  a graph of the theoretical upper sigma in the absence of bkg
 * @param   mr  Mass range (NULL : default one)
 * @param   cl confidence level
 * @param   after cuts   assuming a CLs anaylsis and consider all cuts
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfTheoreticalUpperSigma(XeRange* mr=NULL
                     ,double cl=DEFAULT_CL,bool afterCuts=true, int plot=NONE);


 /**
 * Return a Multigraph of the shape of S1 distribution for the signal,
 * neutron and electron background
 * @param plot either NONE, AUTO, LINEAR, or LOG
*/
    XeMultiGraph  *newMultiGraphOfS1Distributions( int plot=NONE);

// =======================================================================
// ============= Methods to access Signal S1 information =================
// =======================================================================

/**
 * Return the spectrum of Signal S1, for current interaction, mass ,..
 * @param tLEff "t-value" of L-effective
 */
    double *getSignalS1Spectrum(double tLEff);

/**
 * Return the spectrum of Signal S1, for current interaction, mass ,..
 * @param tindex index of the tabulation
 * */
    double *getSignalS1Spectrum(int tindex=SIDE_LEFF_TABULATED);


/**
 * Return the distribution of Signal S1, for current interaction, mass ,..
 * @param tLEff "t-value" of L-effective
 */
    double *getSignalS1Distribution(double tLEff=0.);

 /**
 * get the spectrum of Signal S1
 * @return  a pointer to a new XeSpectrum
 * @param tLEff "t-value" of L-effective
 */
    XeSpectrum* getSignalS1XeSpectrum(double lLEff);

 /**
 * get the spectrum of Signal S1
 * @return  a pointer to a new XeSpectrum
 * @param tindex index of the tabulation
 */
    XeSpectrum* getSignalS1XeSpectrum(int tindex=SIDE_LEFF_TABULATED);


 /**
 * get the distribution of Signal S1
 * @return  a pointer to a new Tabultated Dist
 * @param tLEff "t-value" of L-effective
 */
    TabulatedDist* getSignalS1XeDist(double tLEff=0.);

 /**
 * Create  a graph of the shape of S1 distribution for a signal
 * @param   tIndex  Index of LEff t-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph *newGraphOfSignalS1Distribution(
        int plot=NONE, int tIndex=SIDE_LEFF_TABULATED);

 /**
     * Create  an histogram  of the shape of S1 distribution for a signal
     * @return pointer to newly created Histogram
     * @param   tIndex  Index of LEff t-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F* newHistogramOfSignalS1Distribution(
        int plot=NONE, int tIndex=SIDE_LEFF_TABULATED);

 /**
 * Create  an histogram  of the shape of S1 distribution for a signal
 * @return pointer to newly created Histogram
 * @param   tIndex  Index of LEff t-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH2F* new2DHistogramOfSignalS1Distribution(XeRange *mr=NULL
              ,int plot=NONE, int tIndex=SIDE_LEFF_TABULATED);


 /**
 * Create  a Multigraph of the shape of S1 distribution for various masses
 * @return pointer to newly created XeMultiGraph
 * @param Mr  Mass Range (XeRange) if NULL, use default one
 * @param   tIndex  Index of LEff t-value
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph* newMultiGraphOfSignalS1Distribution(XeRange *mr=NULL
                     ,int plot=NONE, int tIndex=SIDE_LEFF_TABULATED);
 /**
 * create a graph of S1 spectrum for the expected signal
 * @return Pointer to a new XeGraph of S1
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph *newGraphOfSignalS1Spectrum(
           int plot=NONE, int tindex=SIDE_LEFF_TABULATED);

 /**
 * create an histogram of S1 spectrum for the expected signal
 * @return Pointer to a new XeGraph of S1
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F *newHistogramOfSignalS1Spectrum(
           int plot=NONE, int tindex=SIDE_LEFF_TABULATED);

 /**
 * create a 2D histogram of S1 spectrum for the expected signal
 * @return Pointer to a new TH2F
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH2F *new2DHistogramOfSignalS1Spectrum(XeRange *mr=NULL
                         ,int plot=NONE, int tindex=SIDE_LEFF_TABULATED);
 /**
 * Create a multigraph of S1 distribution for various LEff t-values
 * @return Pointer to a new XeMultiGraph of S1
 * @param s1Max  Max S1 scale
 * @param step distance between 2 t-index
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph*newMultiGraphOfSignalS1Spectrum(int step=1, int plot=NONE);

// =======================================================================
// ========= Methods to access Background(s) S1 information ==============
// =======================================================================

/**
 * Return the S1 distribution for a given type of background
 * @param background : ER_BACKGROUND, ER_LEAKAGE , NR_BACKGROUND 
 * or ALL_BACKGROUNDS
 */
    double *getBackgroundS1Distribution(int background);

 /**
 * get the S1 distribution (normalized to 1) for a given type of background
 * @return  a pointer to TabulatedDist
 * @param background : ER_BACKGROUND, ER_LEAKAGE , NR_BACKGROUND 
 * or ALL_BACKGROUNDS
 */
    TabulatedDist *getBackgroundS1XeDist(int background);

/**
 * Return the S1 spectrum for a given type of background
 * @param background : ER_BACKGROUND, ER_LEAKAGE , NR_BACKGROUND 
 * or ALL_BACKGROUNDS
 */
    double *getBackgroundS1Spectrum(int background);

 /**
     * get the S1 spectrum  for a given type of background
     * @return  a pointer to XeSpectrum
     * @param background : ER_BACKGROUND, ER_LEAKAGE , NR_BACKGROUND 
     * or ALL_BACKGROUNDS
 */
    XeSpectrum *getBackgroundS1XeSpectrum(int background);

 /**
 * Create a graph of the shape of S1 distribution for a given background
 * @return pointer to newly created XeGraph
 * @param background : ER_BACKGROUND, ER_LEAKAGE , NR_BACKGROUND 
 * or ALL_BACKGROUNDS
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfBackgroundS1Distribution(int background, int plot=NONE);

 /**
 * Return histogram of the shape of S1 distribution for a given background
 * @return pointer to newly created Histogram
 * @param background : ER_BACKGROUND, ER_LEAKAGE , NR_BACKGROUND 
 * or ALL_BACKGROUNDS
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F* newHistogramOfBackgroundS1Distribution(int background, int plot=NONE);

// =======================================================================
// ======= Methods to access Combined Backgrounds S1 information =========
// =======================================================================

/**
 * Return the spectrum of the background for a given band
 * @return pointer to Spectrum
 * @param  band  from 0 to nBands-1; if 'ALL', means all of them
 */
   XeSpectrum* getAllBackgroundsInBandsS1XeSpectrum(int band=ALL);

/**
 * Return the distribution of the background for a given band
 * @return pointer to TabulatedDist
 * @param  band  from 0 to nBands-1; if 'ALL', means all of them
 */
   TabulatedDist* getAllBackgroundsInBandsS1XeDist(int band=ALL);

/**
 * Return  array o fthe spectrum of the background for a given band
 * @return pointer to double
 * @param  band  from 0 to nBands-1; if 'ALL', means all of them
 */
   double* getAllBackgroundsInBandsS1Spectrum(int band=ALL);

/**
 * Return the array of distribution of the background for a given band
 * @return pointer to double 
 * @param  band  from 0 to nBands-1; if 'ALL', means all of them
 */
   double* getAllBackgroundsInBandsS1Distribution(int band=ALL);

/**
 * Produce a Graph of the S1 distribution of the background
 * @return A pointer to the XeGraph
 * @param  band  from 0 to nBands-1; if 'ALL', means all of them
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
   XeGraph* newGraphOfAllBackgroundsInBandsS1Distribution(int band=ALL
   , int plot=NONE);

/**
 * Produce a multi graph of the S1 distribution of the background
 * for all bands
 * @return A pointer to the XeMultiGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeMultiGraph* newMultiGraphOfAllBackgroundsInBandsS1Distribution(
       int plot=NONE);

/**
 * Produce a Graph of the S1 spectrum of the background
 * @return A pointer to the XeGraph
 * @param  band  from 0 to nBands-1; if 'ALL', means all of them
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
 XeGraph* newGraphOfAllBackgroundsInBandsS1Spectrum(int band=ALL,int plot=NONE);


/**
 * Produce a multi graph of the S1 spectrum of the background
 * for all bands
 * @return A pointer to the XeMultiGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
 XeMultiGraph* newMultiGraphOfAllBackgroundsInBandsS1Spectrum(int plot=NONE);

// =======================================================================
// ===== Methods to access Signal + all Backgrounds S1 information =======
// =======================================================================

/**
 * Create Spectrum of signal + background
 * @return Pointer to a newly created XeSpectrum
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 */
 XeSpectrum *newSpectrumOfSignalAndBackgroundS1(int tindex=SIDE_LEFF_TABULATED);

/**
 * Create Distribution of signal + background
 * @return Pointer to a newly created TabulatedDist
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 */
 TabulatedDist *newDistributionOfSignalAndBackgroundS1(
              int tindex=SIDE_LEFF_TABULATED);

 /**
 * create a graph of S1 spectrum for the expected signal and all backgrounds
 * @return Pointer to a newly created XeGraph 
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
  XeGraph *newGraphOfSignalAndBackgroundS1Spectrum(
          int plot=NONE, int tindex=SIDE_LEFF_TABULATED);

 /**
 * create a graph of S1 spectrum for the expected signal and all backgrounds
 * @return Pointer to a newly created XeGraph 
 * @param tindex  which of the tabulated LEff t values is to be plotted.
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
  XeGraph *newGraphOfSignalAndBackgroundS1Distribution(
         int plot=NONE, int tindex=SIDE_LEFF_TABULATED);

 /**
  * Decide if Bands are produced by AmBe or MC sample, true is MC, false AmBe.
  */
  void  setBandingFromMC(bool Do_mc_SignalBand) {DoMCSignalBand = Do_mc_SignalBand;};

 /**
   * Returns True if bands are computed from simulated MC signal, false if computed from AmBe
  */ 
  bool isBandFromMC() {return DoMCSignalBand;};
 
/**  
 *  print this run
 *  @param level print level
 */
    bool   printIt(int level=1);
/**
 * As requested by being a p-value, 
 * recompute what is needed when smthg has changed
 */
    virtual bool update();
 
/**
 * Get pointer to the data root file
*/
   TFile * getDataFile() {return dataFile; }
   TFile * getSignalFile() {return signalFile; }

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

 /**
     * print flags for tracing purpose
 */
    void traceTheFlags();


  protected :


 /**
     * check that all S1 limits (data, cuts, user analysis) are consistent
     * @return true if ok, false if not
 */
    bool checkS1LimitsConsistency();

/**
 * all operation flags reset to 'undone'
 */
    void resetAllFlags();
 
/**
 * set the smear PMT flag
 * @param smear  do we want to smear? (default is YES)
 */
    void smearPMT(bool smear=true);

/**
 * Instantiate signal bands, if needed  
 */
    void instantiateSignalBandsIfNeeded();

/**
 * Will pLots be in variable limits (i.e. according to min and max) ?
 * @param var true is automatic limits, false in fixed limits
 */
    void setAutomaticLimits(bool var=true);

/**
 * Set the relative normalizaition
 * @param egamma  e/gamma run/DM runs normalization 
 * @param ambe  AmBe run/DM runs normalization 
 * @param darkMatter  Dark Matter/DM runs normalization (usually 1.)
 */
    void setNormalization(double egamma=UNDEFINED , double ambe=UNDEFINED
                         ,double darkMatter=UNDEFINED ) ;

 
/**
 *  Return the standard run name
 *  @param runNumber run number
 */
    static string       getTheName(int runNumber);

/**
 * Initialize members to their default values
 */
    void initialize();
   
/**
 * delete S1S2Bands
 * @param first First S1S2Bands to be deleted
 * @param last Last S1S2Bands to be deleted
 */
    void deleteTheBands(int first=0,int last=N_ALL_DATA_TYPES-1);
   
/**
 * delete a given  S1S2Bands
 * @param data_type  data type (DM_DATA, ...)
 */
    void deleteTheBand(int data_type);

/**
 * Convert a L-Effective t-value into in index and a fraction ,for interplation
 * @param lefft L-Effective "t-value"
 * @return a pair consisting of an index and the fraction to the next index
 */ 
    pair<int,double>    LEffTBinAndFraction(double leffT);

    static int nBands;                    /*!< Number of bands in the analysis*/
    static int nSlices;                  /*!< Number of slices in the analysis*/

    int      runNumber;                                        /*!< Run number*/
    int      dataType;    /*!< Data type, NO_DATA, REAL_DATA or SIMULATED_DATA*/
    int      firstAnalysisS1Bin;           /*!< Lowest S1 bin for the analysis*/
    int      lastAnalysisS1Bin;           /*!< Highest S1 bin for the analysis*/

    bool     backgroundsComputed;  
             /*!< Have all background been computed ?*/
    bool     dataBandsAreBuilt;  
             /*!< Have Dark matter and calibration bands been filled?*/
    bool     signalTabulated;                        
             /*!< Is signal tabulated?*/
    bool     equallyFilledS1Slices; 
             /*!<Do we want S1 slices defined by equal content? (default:no)*/
    bool     DoMCSignalBand; 
             /*!<Bands from MC signal or from AmBe? default AmBe (default:false)*/
    bool     forcedSignalERecoilSpectrum;           
             /*!< Force spectrum of SignalERecoil*/
    bool     smearingPMT;   
             /*!< Are PMT signals to be smeared ? (default:yes)*/

    double   AmBeNormalizationToEvents;        /*!< Normalisation of AmBe runs*/
    double   AnalysisS1Max;                   /*!< Maximum S1 for the analysis*/
    double   AnalysisS1Min;                   /*!< Minimum S1 for the analysis*/
    double   DarkMatterNormalization;            /*!< Normalisation of DM (1.)*/
    double   EGammaNormalizationToEvents;        /*!< Normalisatoin of ER runs*/
    double   Exposure;                             /*!< Run exposure (in days)*/
    double   FiducialMass;                          /*!< Fiducial mass (in kg)*/
    double   LEffErMin;                                  /*!< Hard cut of LEff*/
    double   LightYield;                         /*!< Light yield  (in pe/Kev>*/
    double   S2OverS1Max;           /*!< upper limit of S2/S1 for the analysis*/
    double   S2OverS1Min;           /*!< Lower limit of S2/S1 for the analysis*/
    double   Se;                    /*!< Field quenching for Electronic recoil*/
    double   SelectionS1Max;         /*!< Maximum S1 in smeared selection cuts*/
    double   SelectionS1Min;         /*!< Minimum S1 in smeared selection cuts*/
    double   SelectionUS1Min;      /*!< Minimum S1 in unsmeared selection cuts*/
    double   SigmaPMT;                   /*!< Sigma when smearing PMT response*/
    double   Sr;                       /*!< Field quenching for Neutron recoil*/
    double   typicalErMin;  /*!< Ermin,for fast evaluation of number of events*/
    double   sigmaForReports;   /*!<Cross section used in plots,graphs,reports*/
    double   sigmaFactor;               /*!<sigmaForReports/getSigmaNucleon()*/


    vector<vector<S1S2Bands*>> 	    SignalBand;  
	 /*!< Holder of signal bands for different Leff and Qy, first vector runs on Leff and the second on Qy*/
    
    MCSignalModel *signalHandler;  //this is just  temporary for test.

    double        *SignalERecoilSpectrum;    
                   /*!< pointer to current the SignalERecoil spectrum*/
    double        *SignalERecoilDistribution; 
                   /*!<pointer to current SignalERecoil distribution*/
    XeSpectrum    *SignalERecoilXeSpectrum;              
                   /*!<Current SignalERecoil XeSpectrum*/
    TabulatedDist *SignalERecoilXeDist;
                   /*!<SignalERecoil TabulatedDist*/

    double        *SignalS1Spectrum[N_LEFF_TABULATED];  
                   /*!< pointers to expected signal S1 spect. for various Leff*/
    double        *SignalS1Distribution[N_LEFF_TABULATED]; 
                   /*!< pointers to expected signal S1 dist. for various Leff*/
    XeSpectrum    *SignalS1XeSpectrum[N_LEFF_TABULATED];   
                   /*!<p Signal S1 Spectrum, for various Leff*/
    TabulatedDist *SignalS1XeDist[N_LEFF_TABULATED];
                   /*!<Signal S1 distribution, for various Leff*/
    double         SignalS1Total[N_LEFF_TABULATED];     
                   /*!< Total events for various Leff*/

    double        *BackgroundS1Spectrum[N_BACKGROUNDS];        
                   /*!< pointer to S1 spectrum for various backgrounds*/
    double        *BackgroundS1SDistribution[N_BACKGROUNDS];        
                   /*!< pointer to  S1 distribution for various backgrounds*/
    XeSpectrum    *BackgroundS1XeSpectrum[N_BACKGROUNDS];
                   /*!<S1 spectrum for various backgrounds */
    TabulatedDist *BackgroundS1XeDist[N_BACKGROUNDS];       
                   /*!<S1 distributions for various backgrounds*/
    double         TotalContent[N_BACKGROUNDS];     
                   /*!< Total events for each background*/
    vector<double> BackgroundsInBands[N_BACKGROUNDS];
                   /*!<Nominal background count for each band*/

    vector<XeSpectrum*>  AllBackgroundsInBandsS1XeSpectrum;
                         /*!<Spectra of all backgrounds for each band*/
    vector<double*>      AllBackgroundsInBandsS1Spectrum;
                         /*!<Spectra of all backgrounds for each band*/
    vector<TabulatedDist*> AllBackgroundsInBandsS1XeDist;
                           /*!<Distribution of all backgrounfs for each band*/
    vector<double*>        AllBackgroundsInBandsS1Distribution;
                           /*!<Distribution of all backgrounfs for each band*/


    XeSetOfSelectionCuts *selectionCuts; 
                   /*!< pointer to set of selection cuts*/
    XeSetOfCuts   *darkMatterCuts; 
                   /*!< pointer to set of dark matter cuts*/
    XeSetOfCuts    allCuts;                    
                   /*!< Set of all the cuts*/

    PValue         *pValue;             /*!<pointer to PValue currently in use*/
    Target         *target;                       /*!<pointer to Target in use*/
    LEff           *leff;            /*!<pointer to description of L-Effective*/
    Qy             *qy;            /*!<pointer to description of Qy*/
    RunFlattener   *flattener;                   /*!<pointer to S1S2 flattener*/
    S1S2Data       *s1s2[N_CUT_DATA_TYPES];  /*!<pointer to all S1S2 data sets*/
    S1S2Bands      *referenceBands;        /*!<pointer to S1S2 reference bands*/
    S1S2Bands      *bands[N_ALL_DATA_TYPES];     /*!<pointer to all S1S2 bands*/

    ElectronBackgroundModel *ebm;     /*!<pointer to electron background model*/
    NeutronBackgroundModel  *nbm;      /*!<pointer to neutron background model*/
   
    TFile	*dataFile;     /*pointer to the stored data root file*/
    TFile	*signalFile;   /*pointer to the stored MC signal root file*/

      ClassDef(XeRun,1)
};

/**
   * virtual class  for calculating pvalues from S1S2 
*/
class S1S2PValue : virtual public PValue, virtual public RunComponent {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public :

     virtual       ~S1S2PValue();
 /**
     * return  New Instance of a virtual PValue
     * @param   name  its name
     * @param   mode  either CUTS_ANALYSIS or PL_ANALYSIS
     * @param   run   XeRun containing the data
 */
                    S1S2PValue(string name, int mode, XeRun* run=NULL);
     void           setInteraction(Interaction* inter);
     void           setWimpMass(double mass);
     void           updateSigmaScale();

/**
 * return expected number of signal events for 1 cm2 cross section
*/
     double         nSignalPerCm2();            

  protected :

     string         shortName;
     int            requestedMode;

     bool           checkEverything();

      ClassDef(S1S2PValue,1)

} ;

/**
   * pvalue for S1S2 CLs
*/

class S1S2CLs : public S1S2PValue {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   public :

                   ~S1S2CLs();
                    S1S2CLs(XeRun* run=NULL);
     double         pValueS(double nEvents);
     bool           checkPValue();
     void           printFlagsAndParameters();

     bool           update();
   protected:

     double         bkg;
     int            nEvents;

      ClassDef(S1S2CLs,1)


};

/**
   * pvalue for S1S2 Profile Likelihood
*/
class S1S2PL : virtual public ProfileLikelihood, virtual public RunComponent {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

  public :

    static constexpr  bool DEFAULT_PAR_STAT_BKG_FIT=false;
    static constexpr  bool DEFAULT_PAR_SYST_BKG_FIT=false;
    static constexpr  bool DEFAULT_PAR_LEFF_TVALUE_FIT=true;
    static constexpr  bool DEFAULT_PAR_QY_TVALUE_FIT=false;
    static constexpr  bool DEFAULT_S1_LIKELIHOOD=true;

   ~S1S2PL();
 /**
     * @param  run   The run (real or simulated) from which data are takem
 */
    S1S2PL(XeRun* run);

 /**
     * print the values of flags (fit LEff, fit Bkg, with S1 Likelihood...)
 */
    void printFlagsAndParameters();

 /**
     * initialize all members and set the flags (fit LEff, fit Bkg, with S1 Likelihood...) to their default values
 */
    void initialize();

/**
 * Set the wimp mass and trigger the new for a new S1 spectrum computation
 * @param mass the new WIMP mass, in GeV/c2
 */
    void setWimpMass(double mass);

    double getWimpMass(); 

 /**
     * set the flag to fit LEff t-value
     * @param  fit  Flag specifying whether LEff should be considered as a nuisance parameter
 */
    void fitLEffTValue(bool fit=DEFAULT_PAR_LEFF_TVALUE_FIT);

 /**
     * set the flag to fit Qy t-value
     * @param  fit  Flag specifying whether Qy should be considered as a nuisance parameter
 */
    void fitQyTValue(bool fit=DEFAULT_PAR_QY_TVALUE_FIT);

 /**
     * set the flag to fit the systematic background t-value
     * @param  fit  Flag specifying whether systematic Bkg should be considered as a nuisance parameter
 */
    void fitSystBkgTValue(bool fit=DEFAULT_PAR_SYST_BKG_FIT);

 /**
     * set the flag to fit the systematic background t-value
     * @param  fit  Flag specifying whether systematic Bkg should be considered as a nuisance parameter
 */
    void fitStatBkgTValue(bool fit=DEFAULT_PAR_STAT_BKG_FIT);

 /**
     * Set the flag "with or without" S1 shape likelihhod 
     * @param  s1  Flag specifying whether the shape of S1 distribution enters the likelihood
 */
    void   withS1Likelihood(bool s1=DEFAULT_S1_LIKELIHOOD);

    void   fitSignalAcc(bool doFit) {isSigAccFit = doFit;};

/**
 * return expected number of signal events for 1 cm2 cross section
 * @param lt L-Eff t-value
 */
    double nSignalPerCm2(double lt);            

/**
 * return expected number of signal events for 1 cm2 cross section
*/
    double nSignalPerCm2();            

 /**
     * Compute the Log Likelihood  
     * @return  Log Likelihood
     * @param  none  (all parameters are set thru LKParameter)
 */
    double computeTheLogLikelihood();

/**
  * SHOULD be called each time a change of condition happens, in practice is called each time 
   the wimp mass change. Fill internal tables from the run band content.
*/
    bool   update();

/**
  * Simulate an outcome
*/
   bool simulate(double sigma=0.);

 /**
 * get the systematic background t-value parameter
 */
    LKParameter* getSystBkgTValueParameter();

 /**
 * get the statistical background t-value parameters
 */
    vector<LKParameter*> getStatBkgTValueParameters();

 /**
 * get the LEff t-value parameter
 */
    LKParameter* getLEffTValueParameter();

 /**
 * get the Qy t-value parameter
 */
    LKParameter* getQyTValueParameter();

 /**
 * get the background sigma parameter
 */
    LKParameter* getSigmaParameter();


 /**
 * Things to be done Post exclusion finder 
 */
    void    exclusionComputed();

/**
 * get the cross section at exclusion
*/
    double  getSigmaForExclusion();

/**
 * get LEff t-value at exclusion
*/
    double  getLEffTValueForExclusion();

/**
 * get Delta LEff t-value at exclusion
*/
    double  getDLEffTValueForExclusion();

/**
 * get Qy t-value at exclusion
*/
    double  getQyTValueForExclusion();

/**
 * get Delta Qy t-value at exclusion
*/
    double  getDQyTValueForExclusion();

/**
 * get Background statistical t-value at exclusion
*/
    vector<double> getStatBkgTValueForExclusion();

/**
 * get Delta Background statistical t-value at exclusion
*/
    vector <double>  getDStatBkgTValueForExclusion();

/**
 * get Background systematic t-value at exclusion
*/
    double   getSystBkgTValueForExclusion();

/**
 * get Delta Background systematic t-value at exclusion
*/
    double  getDSystBkgTValueForExclusion();

/**
 * get effect of one systematic sigma in background
*/
    double  getSystBkgSigma();

/**
 * set effect of one systematic sigma in background
*/
    void    setSystBkgSigma(double sigma);

    void setData(int dataType);

    void generateAsimov(double mu_prime) ;

    void generateToyDataset(double seed, double mu_prime);

    double getSignalDefaultNorm() { return run->getSignalDefaultNorm(); };

    double getSignalMultiplier()  {return  run->getSignalMultiplier(); };
    void   setSignalMultiplier(double val) { run->setSignalMultiplier(val); };


/**
 * print current values
 * @param level print level
 */
    bool printIt(int level=1);

/**
 * Set if you'd like to have real data or extracted randomly from bkg, you need to previously generate them in this case
 * call XeRun::generateData() 
*/
   void setWithSimulatedata(bool flag) { withSimulatedData = flag;};

  protected :


    int     nBands;                                        /*!<Number of bands*/
    double  frac_Co_model;		   /*!<normalizzation factor between bkg model content and calibration sample*/
    double  N_bkg_tot;			   /*!<Total Bkg model content in all bands */
    bool    isLEffTValueFit;               /*!<Do we want to fit LEff t-value?*/
    bool    isQyTValueFit;                 /*!<Do we want to fit Qy   t-value?*/
    bool    isSystBkgTValueFit; /*!<Do we want to fit background syst t-value?*/
    bool    isSigAccFit; /*!<Do we want to fit background syst t-value?*/
    bool    isStatBkgTValueFit; /*!<Do we want to fit background stat t-value?*/
    bool    S1Likelihood;                       /*!<Fit shape of S1 spectrum ?*/
    double  SigmaForExclusion;         /*!<Value of fitted Sigma for exclusion*/
    double  LEffTValueForExclusion;   /*!<Value of  LEff t-value for exclusion*/
    double  DLEffTValueForExclusion;  /*!< Error on LEff t-value for exclusion*/
    double  QyTValueForExclusion;   /*!<Value of  Qy t-value for exclusion*/
    double  DQyTValueForExclusion;  /*!< Error on Qy t-value for exclusion*/
    double  SystBkgTValueForExclusion; 
                /*!<Value of normalization background systematic tvalue for exclusion*/
    double  DSystBkgTValueForExclusion;
                /*!<Error on background systematic tvalue for exclusion*/
    double  SystBkgSigma; 
                /*!<How much corresponds to one sigma in systematic background*/


    vector <vector<int>>     dataS1Bins;   /*!< Data S1 distro expressed in "S1Bins" (S1 bin number) for each band*/	

    vector <double>  StatBkgTValueForExclusion; 
                            /*!<Value of background stat. tvalue for exclusion*/
    vector <double>  DStatBkgTValueForExclusion;
                            /*!<Error on background stat. tvalue for exclusion*/
    vector <double>  StatBkgSigmas;
                 /*!<How much corresponds to one sigma in stat background band*/

    S1S2Bands 	    *BackgroundBand; 
        /*!< pointer to the current bkg band, NOTE: they are fixed for each signal mass*/

    S1S2Bands 	    *DataBand; 
        /*!< pointer to the current data band, NOTE: this can be set with setData()*/

    LKParameter     *le;             /*!< pointer to likelihood LEff parameter*/
    LKParameter     *lqy;	     /*!<pointer to likelihood Qy parameter*/
    LKParameter     *ls;             /*!<pointer to likelihood sigma parameter*/
    LKParameter     *lkSystBkg;      /*!< pointer to likelihood background parameter*/
    LKParameter     *lkacc;         /*!< pointer to likelihood signal acceptance parameter*/
    vector <LKParameter*>     lkStatBkgs;
              /*!< pointers to likelihood background in band t-value parameter*/
    bool    	    withSimulatedData;

    int 	    data_type;

    ClassDef(S1S2PL,1)

};


/**
   * Class for ER background imported from a TH2F.
   * Important:: the histograms has to be normalized to the expected number of events.	
*/
class ImportedElectronBackground : public ElectronBackgroundModel {

  public:

/**
 *  constructor
 *  @param  run  XeRun on which the model is applied
 *  By default the input file is the Run file and the default histos name are:
 *  RunXX_ERgauss (gauss bkg) and RunXX_ERleak (leakage).
 */
     ImportedElectronBackground(XeRun* run);

/**
 * No arg constructor for root
 */
     ImportedElectronBackground();

    ~ImportedElectronBackground();

 /**
    * Returns the Bands for Gauss bkg filled from a TH2F 
 */
     S1S2Bands* computeBands();

 /**
    * Returns the Bands for leakage bkg filled from a TH2F 
 */
     S1S2Bands* computeAnomalousBands();

          
/**  
 *  print this background model
 *  @param level print level
 */
     bool printIt(int level=1);   

/** 
 * Implementation of the virtual method to define run compatibility
 * just for historical reason... There is no way to implement it for this class 
 * Return always true.
 */
  bool isRunCompatible(int runNumber);

/**
  * Set name of histogram for gauss bkg.
*/
  void setHistoNameGauss(TString name) {gaussBkg.setHistoName(name);};

/**
  * Set name of histogram for leakage bkg.
*/
  void setHistoNameLeak(TString name)  {leakageBkg.setHistoName(name);};

  void scaleLeak(double norm)   {leakageBkg.scale(norm);};
  void scaleGauss(double norm)  {gaussBkg.scale(norm);};

  void rebinGauss(int b)        {gaussBkg.rebinX(b);};
  void rebinLeak(int b)        {leakageBkg.rebinX(b);};

  TH2F getHistoGauss(){return gaussBkg.getModel();};
  TH2F getHistoLeak(){return leakageBkg.getModel();};


  private:
 
    ModelImporter  gaussBkg;	/*<!gauss bkg importer*/
    ModelImporter  leakageBkg;  /*<!leakage bkg importer*/

    ClassDef(ImportedElectronBackground,1)

};


/**
   * Neutron background model Imported from TH2F
*/
class ImportedNeutronBackgroundModel : public NeutronBackgroundModel {

  public:

     ~ImportedNeutronBackgroundModel();

 /**
     * constructor
     * @param  run  XeRun on which the model is applied
     *  By default the input file is the Run file and the default histos name are:
     *  RunXX_Neutron.
 */
     ImportedNeutronBackgroundModel(XeRun* run);

/**
 * Empty constructor for ROOT
 */
     ImportedNeutronBackgroundModel();

 /**
     * Return S1S2Bands computed from TH2F
     * @return  Bands for (regular) background
 */
     S1S2Bands* computeBands();


 /**
     * return  estimate in neutrons in Events/Kg/Day/Pe
     * @param pe  upper equivalent number of photo electrons
 */
     double expectedBelow(double pe);

/** 
 * Implementation of the virtual method to define run compatibility
 * This cannot be implementd for this class. Returns always true.
 */
  bool isRunCompatible(int runNumber);

/**
  * Set name of histogram for leakage bkg.
*/
  void setHistoName(TString name)  {neutronBkg.setHistoName(name);};

  void scale(double norm)  {neutronBkg.scale(norm);};

  TH2F getHisto(){return neutronBkg.getModel();};

  private :   

    ModelImporter  neutronBkg;/*<!neutron bkg importer*/

      ClassDef(ImportedNeutronBackgroundModel,1)

};





#endif

