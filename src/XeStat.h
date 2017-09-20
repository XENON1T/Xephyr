#ifndef XeStat_h
#define XeStat_h

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

typedef double (*func_d)(double);

class DataSet;

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

    ClassDef(XeStat,1)
} ;


/**
   * extension to map<double,int>
*/
class XeValues : public map<double,int>, public XeObject {

  public : 

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    virtual ~XeValues();

     XeValues();
 /**
     * basic constructor 
     * @param  name  name of the XeObject
 */
     XeValues(string name);

 /**
     * Constructor from an array
     * @param  name  name of the XeObject
     * @param  d array of values
     * @parame n  number of entries
 */
    XeValues(string name,double *d,int n);

 /**
     * Constructor from a vector
     * @param  name  name of the XeObject
     * @param  vd vector of values
 */
    XeValues(string name,vector<double>& vd);


 /**
     * Constructor from a dataSet
     * @param  name  name of the XeObject
     * @param  dataSet the data set considered
     * @param  col    which column`
 */
    XeValues(string name,DataSet* dataSet,int col);

 /**
     * Return number of different values 
 */
    int  getNValues();

 /**
     * Return  number of original entries
 */
    int  getNEntries();

/**
  * print the object
  * @param level print level (0-> no print)
*/
    bool printIt(int level=1);

 /**
     * Return  largest gap between the values
 */
    double largestGap();

/**
   * Return quantile of the distrubtion
   * @return x value of the smallest values for a given fraction 
   * @param f fraction considered
*/
    double quantile(double f);

/**
   * Return x value above or below a certain number of standard deviations
   * @param sigma how many sigmas above/below (e.g. +1 ->.68)
*/
    double findSigma(double sigma);


    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    void reset();
                    
 /**
     * adds an array of values
     * @param values address of the array
     * @param n number of points
 */
    void add(double* values,int n); 

 /**
     * adds one value
     * @param value value to arr
 */
    void add(double value); 


  protected :

    int      nEntries;

     ClassDef(XeValues,1)

};

 const string SPECTRUM_HEADER="Spectrum: ";

/**
   * Class describing a spectrum. If normalized, it is called a XeDist
*/
class XeSpectrum : virtual public XeGraphics, virtual public XeStat{

  public :

   virtual ~XeSpectrum();

/**
 *  Empty constructor for ROOT
 */  
    XeSpectrum();

/**
 *  Minimal constructor
 */  
    XeSpectrum(string n);

/**
 * Return the  name with the header
 */
     static string getTheName(string nam);

/**
 * Return the name without the header
 */
     string getReducedName();


 /**
     * Constructor giving bin boundaries and possibly pointer to spectrum
     * @param nam name of the spectrum
     * @param nBins number of bins
     * @param min  Lower bound of the spectrum
     * @param max Upper bound of the spectrum
     * @param spec pointer to the spectrum to be copied (mode SPECTRUM)
     *            , or the values to be histogrammed (mode VALUES)
     * @param nValues number of values
     * @param mode whether spec refer to SPECTRUM or VALUES
 */
    XeSpectrum(string nam,int nBins, double min, double max
              , double* spec=NULL, int nValues=0, int mode=SPECTRUM);

/**
     * Constructor giving Linear Range and possibly pointer to spectrum
     * @param nam name of the spectrum
     * @param lr  Linear Range defining the boundaries
     * @param spec pointer to the spectrum to be copied (mode SPECTRUM)
     *            , or the values to be histogrammed (mode VALUES)
     * @param nValues number of values
     * @param mode whether spec refer to SPECTRUM or VALUES
 */
    XeSpectrum(string nam,LinearRange *lr, double *spec=NULL
              ,int nValues=0, int mode=SPECTRUM);

 /**
     * Constructor giving bin boundaries and passing vector of spectrum
     * @param nam name of the spectrum
     * @param nBins number of bins
     * @param min  Lower bound of the spectrum
     * @param max Upper bound of the spectrum
     * @param spec pointer to the spectrum to be copied (mode SPECTRUM)
     *            , or the values to be histogrammed (mode VALUES)
     * @param mode whether spec refer to SPECTRUM or VALUES
 */
    XeSpectrum(string nam,int nBins, double min, double max
              , vector<double> &spec, int mode=SPECTRUM);

/**
     * Constructor giving Linear Range and possibly pointer to spectrum
     * @param nam name of the spectrum
     * @param lr  Linear Range defining the boundaries
     * @param spec pointer to the spectrum to be copied (mode SPECTRUM)
     *            , or the values to be histogrammed (mode VALUES)
     * @param mode whether spec refer to SPECTRUM or VALUES
 */
    XeSpectrum(string nam, LinearRange *lr,vector<double> &spec
               , int mode=SPECTRUM);

 /**
     * Constructor from an XeValues of realizations
     * @param  nam  name of the original quantity
     * @param  nBins Number of bins
     * @param  min  lowest value of lowest bin
     * @param  max  highest value of highest bin
     * @param  values  Original table of Values (actually XeValues)
 */
    XeSpectrum(string nam,int nBins,double min,double max, XeValues* values);

/**
 * Fill content by a (bounded) exponential distribution
 * @param x0  1./slope
 * @param norm overall normalization
 * @param xmin lower cut
 * @param xmax upper cut
 */
    void fillExponent(double x0, double norm=1., double xmin=0.
                     ,double xmax=VERY_LARGE);

 /** 
 * add two Spectra
 */
   XeSpectrum operator+(XeSpectrum& other);

 /** 
 * multiply one double times one Spectrum
 */
   XeSpectrum operator*(double a);


/**
 * Add (cumulate) contents with those of another spectum
 * @param other  pointer to the other XeSpectrum
 */
    void add(XeSpectrum* other);

 /**
     * Fill one value with a given weight
     * @param   x value
     * @param   weight weight (default 1.)
 */
    void   fillDelta(double x, double weight=1.);

 /**
     * set content of a bin
     * @param   bin which bin, from 0 to nBins-1
     * @param   content what to set
 */
    void   setBin(int bin, double content);


/**
  * print the object
  * @param level print level (0-> no print)
*/
    bool   printIt(int level=1);

/**
 * print a detailed table of content
 * @param ignoreTail  don't print the uppermost empty part of the spectrum
 */
  void printDetailed(bool ignoreTail);

 /**
     * Return the sum of contents
 */
    double sum();


/**
    * Return the range (actually LinearRange) of x values
*/
    LinearRange* getRange();

/**
    * Return the number of bins
*/
    int          getNBins();

/**
    * Return the lowest edge
*/
    double       getMin();

/**
    * Return the highest edge
*/
    double       getMax();


/**
    * Return the  contents vector
    * @return pointer to the array of contents
*/
   vector<double>* getTheSpectrum();

/**
    * Return the  contents
    * @return pointer to the array of contents
*/
    double*      getSpectrum();

 /**
     * Return the bin number corresponding to a given value
     * @param  x x-value to be considered
 */
    int getBin(double x);

 /**
     * Return the bin number and fraction corresponding to a given value
     * @return a pair containing the bin and the fraction
     * @param  x x-value to be considered
 */
    pair<int,double> getBinAndFraction(double x);

 /**
     * computes the integral of the spectrum between two values 
     * @return sum of (fractional) bin contents
     * @param x1 lower bound of the integral
     * @param x2 upper bound of the integral
 */
    double integral(double x1, double x2);

 /**
     * @return A graph (actually an XeGraph) of the spectrum 
     * @param Xlab   label of x-axis
     * @param Ylb    label of y-axis
     * @param style  Style of the graph
     * @param legend Legend in case of a multi-graph
 */
    XeGraph* newGraph( string Xlab="", string Ylab="" , XeStyle* style=NULL
                     , string legend="");

 /**
     * Draw it as a graph
     * @param Xlab   label of x-axis
     * @param Ylb    label of y-axis
     * @param style  Style of the graph
     * @param options  TGraph->Draw() option
 */
    void drawGraph(string options="C",string Xlab="", string Ylab="" 
                   , XeStyle* style=NULL);
 /**
     * Draw it as a graph with its frame
     * @param Xlab   label of x-axis
     * @param Ylb    label of y-axis
     * @param style  Style of the graph
     * @param options  TGraph->Draw() option
 */
    void drawGraphWithFrame(string options="C",string Xlab="", string Ylab="" 
                   , XeStyle* style=NULL );

 /**
     * Draw it as an histogram
     * @param options  THist->Draw() option
 */
    void drawHistogram(string options="C",string Xlab="", string Ylab="");


 /**
 * @return an histogram
 * @param  name  Name of the histogram
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F* newHistogram(string name="",int plot=NONE);

 /**
 * @return an empty 2-D histogram form an array of spectra
 * @param  nSpectra  number of spectra
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH2F* newHistogram(int nSpectra,int plot=NONE);

 /**
     * fill one of the 2D histogram slices with the current spectrum
     * @param hist  The histogram to be filled
     * @param spectrum numbered from 0 to nSpectra-1
 */
    void fillHistogram(TH2F* hist, int spectrum);

/**
   * Reset the content
*/
    void reset();

  /**
   * Build the spectrum from an list of realizations
   * @param values array of values
   * @param nValues number of points
 */
   virtual void fillValues(double* values, int nValues);

                   
 /**
   * Copy the spectrum from an array
   * @param spec  values to be copied
 */
   virtual void  importSpectrum(double* spectrum);

 /**
     * copy the spectrum from a vector
     * @param spectrum  vector of values to be copied
 */
   virtual void  importSpectrum(vector<double>& spectrum);

 /**
     * add the spectrum from a vector
     * @param spectrum  vector of values to be added
 */
   virtual void  addSpectrum(vector<double>& spectrum);

/**
 * Internal method: set the range
 * @param: range pointer to LinearRange
 */
   void setTheRange(LinearRange *range);

/**
 * Internal method: set the range. Create a Linear Range only when needed
 * @param nB number of bins
 * @param min lower edge of range
 * @param max upper edge of range
 */
   void setTheRange(int nB, double min,double max);


  protected :

/** 
 * initialize members
 */
  void initialize();

    int            nBins;                                 /*!< number of bins */
    bool           rangeIsCreated;       /*!< a new XeRange had to be created */
    LinearRange *  range;   /*!< pointer to the XeRange (actually LinearRange)*/
    vector<double> spectrum;                                /*!< the values   */

    ClassDef(XeSpectrum,1)
};

class XeDist : virtual public XeStat  {

   public:

     virtual ~XeDist();

     XeDist();

 /**
     * Basic constructor
 */
              XeDist(string n);

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    

 /**
     * Normalize histograms the way I like
     * @param h      histogram to be normalized
     * @param toWhat requested sum of histogram after normalisation
 */
     static  void   normalize(TH1F* h,double toWhat=1.);

 /**
     * Normalize TGraph the way I like
     * @param graph     graph  to be normalized
     * @param toWhat requested sum of graph after normalisation 
     * (1.-> makes it a distribution!)
 */
     static  void   normalize(TGraph* graph,double toWhat=1.);

 /**
   * Return quantile of an histogram
   * @return value of the smallest values for a given fraction 
   * @param histogram to be considered
   * @param q quantile (fraction of distribution considered)
 */
     static  double findQuantile(TH1* h, double q);

 /**
   * Return value for n sigmas above/below in a histogram
   * @return value of the smallest values for a given fraction 
   * @param histogram to be considered
   * @param sigma quantile expressed in +/- n sigmas (e.g - 0> .5, +1 -> 0.68)
 */
     static  double findSigma(TH1* h, double sigma);

/**
 * Wraper to random generator
*/
     static  double rndm();

 /**
     * generates a pseudo random variable according to the distribution  
     * @return A pseudo random variable generated according to the distribution
 */
     virtual double generate()=0;

 /**
     * generates a vector of  pseudo random variables according to the distribution  
     * @param  n number of elements
     * @return A vector of size n
 */
     vector<double> generateVector(int n);

 /**
     * returns probability density function  
     * @return  Probability density function
     * @param   x Point at which is it evaluated
 */
     virtual double pdf(double x)=0;

 /**
     * returns logs of probability density function 
     * @return  Log of probability density function
     * @param   x Point at which is it evaluated
 */
     virtual double logPdf(double x);

 /**
     * returns cumulated probability density function  
     * @return  Cumulated probability density function
     * @param   x Point at which is it evaluated
 */
     virtual double cdf(double x)=0;

 /**
     * computes integrated probability density function between bounds 
     * @return Integrated pdf inside interval
     * @param x1  Lower bound of integration
     * @param x2  Upper bound of integration
 */
     double cdf(double x1,double x2);

/**
 * Return the corrected name
 */
     static string getTheName(string nam);

  protected :
 
     static TRandom3 random;                            /*!< random generator */

     ClassDef(XeDist,1)

};

/**
   * A Generic Tabulated distribution.
   * This is a virtual class. It is a untity normalized XeSpectrum
*/

class TabulatedDist : public XeSpectrum, public XeDist {

  public:
  
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    virtual ~TabulatedDist();  

  /**
    * Empty Constructor for ROOT
 */
    TabulatedDist();  

  /**
    * Constructor with name only
 */
    TabulatedDist(string n);  

 /**
     * Constructor giving Linear Range
     * @param  nam   name of the original quantity
     * @param  lr    LinearRange describing the bins
 */

    TabulatedDist(string nam,LinearRange *lr, double *spec=NULL
                 , int nValues=0, int mode=SPECTRUM);

    double cdf(double x);
    double pdf(double x);
 
 /**
    * Constructor from an XeSpectrum
 */
    TabulatedDist(XeSpectrum *sp);

 /**
    * Constructor from an array of either the spectrum of realizations
     * @param  nam  name of the original quantity
     * @param  nBins Number of bins 
     * @param  min  lowest value of lowest bin
     * @param  max  upper value of upper bin
     * @param spec pointer to the spectrum to be copied (mode SPECTRUM)
     *            , or the values to be histogrammed (mode VALUES)
     *            if NULL, nothing is done
     * @param nValues number of values
     * @param mode whether spec refer to SPECTRUM or VALUES
 */
    TabulatedDist(string nam,int nbins, double min,double max,double *spec=NULL
                 , int nValues=0, int mode=SPECTRUM);

 /**
     * Constructor from an XeValue of realizations
     * @param  nam  name of the original quantity
     * @param  nBins Number of bins 
     * @param  min  lowest value of lowest bin
     * @param  max  upper value of upper bin
     * @param  values  Original vector of values
 */
    TabulatedDist(string nam,int nBins,double min,double max,XeValues* values);


 /**
     * Copy content of distribution from an array, normalize and cumulate
     * @param array  pointer to array of values to be copied
 */
   void importSpectrum(double* array);

 /**
     * Copy content of distribution from a vector, normalize and cumulate
     * @param spectrum  vector of values to be copied
 */
   void importSpectrum(vector<double>& spectrum);

 /**
     * Copy content of distribution from a spectrum, normalize and cumulate
     * @param spectrum  vector of values to be copied
 */
   void importSpectrum(XeSpectrum* spectrum);

  /**
   * Build the distribution from an list of realizations
   * @param values array of values
   * @param nValues number of points
 */
   virtual void fillValues(double* values, int nValues);


 /**
     * Set the spectrum of the distribution from a distribution, 
     * normalize and cumulate. 
     * Note: the original distribution need not be tabulated
     * @param dist mother distribution
 */
   void importDistribution(XeDist* dist);

  /**
     * add one delta in the distribution, normalize and cumulate
     * @param x where
     * @param w weight
 */
   void addDelta(double x, double w=1.);

 /**
 * @return the cumulated histogram
 * @param  name  Name of the histogram
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    TH1F* newCumulatedHistogram(string name="",int plot=NONE);

/**
     * @return A graph (actually an XeGraph) of the cumulated distribution
     * @param Xlab   label of x-axis
     * @param Ylb    label of y-axis
     * @param style  Style of the graph
     * @param legend Legend in case of a multi-graph
 */
    XeGraph* newCumulatedGraph( string Xlab="", string Ylab="" 
                                 , XeStyle* style=NULL, string legend="");

 /**
     * Draw the cumulated distribution as a graph
     * @param Xlab   label of x-axis
     * @param Ylb    label of y-axis
     * @param style  Style of the graph
     * @param options  TGraph->Draw() option
 */
    void drawCumulatedGraph(string options="C",string Xlab="", string Ylab="" 
                   , XeStyle* style=NULL);
 /**
     * Draw the cumulatedDistribution as a graph with its frame
     * @param Xlab   label of x-axis
     * @param Ylb    label of y-axis
     * @param style  Style of the graph
     * @param options  TGraph->Draw() option
 */
    void drawCumulatedGraphWithFrame(string options="C",string Xlab=""
                                    , string Ylab="" , XeStyle* style=NULL );

 /**
     * Draw the cumulatedDistribution as an histogram
     * @param options  THist->Draw() option
 */
    void drawCumulatedHist(string options="C",string Xlab="", string Ylab="");


 /**
     * generate a variable according to the distribution
 */
    double generate();

 /**
     * Convert cdf into x (to find the top 10%, etc...)
     * @return the x value corresponding to the Cdf value x
     * @param c  cdf value (i.e .5=median, .9=top 10%, etc...)
 */
    double quantileX(double c);

 /**
     * Convert cdf into a bin number (to find the top 10%, etc...)
     * @return the x value corresponding to the Cdf value x
     * @param c  cdf value (i.e .5=median, .9=top 10%, etc...)
 */
    int quantileBin(double c);

/**
  * print the distribution
  * @param level print level (if >=3, also print cumulated distribution)
*/
    bool printIt(int level=1);

 /**
     * print the cumulated distribution
     * @param level of print 
*/
    bool printCumulated(int level=1);

 /**
     * Return cumulated values
     * @return pointer to the values
*/
    double* getCumulated();

/**
     * Normalized to 1 then Cumulate the tabulated distribution
 */
    
    void   normalizeAndCumulate();

/*
 * return pointer to the "cumulated values-> bin" map
 */
    map<double,int>* getCdfMap();

/*
 * return pointer to the cumulated values vector
 */
    vector<double>*  getCdfVector();


  protected:
    
    vector<double>  cdfVector;  /*!< List of cumulated probabilities */
    map<double,int> cdfMap;   /*!< Used for generating according to this dist.*/

    ClassDef(TabulatedDist,1)

};

/**
   * Uniform distribution
*/

class UniformDist : public XeDist {
 
  public :
  
/**
 * Default constructor, uniform between 0 and 1 
 */
    UniformDist();
  
/**
 * Regular constructor, uniform between 0 and 1 
 * @param a lower bound
 * @param b upper band
 */
    UniformDist(double a, double b);
   ~UniformDist();

 /**
 * Reset the limits
 * @param a lower bound
 * @param b upper band
 */
    void   setLimits(double a, double b);

/**
 * return the cumulated distribution function
 * @param x the random variable
*/
    double cdf(double x);

/**
 * return the normalized  probability distribution function
 * @param x the random variable
 */
    double pdf(double x);

/**
 * generate a pseudo random variable according to this distribution
 */
    double generate();
    static double generate(double a, double b);

  protected :
   
    double x0;
    double x1;
    double onew;
     
     ClassDef(UniformDist,1)

};

/**
   * Exponential distribution
*/

class ExponentialDist : public XeDist {
 
  public :
  
    ExponentialDist(double x0=1.);
   ~ExponentialDist();
 
    void   setX0(double x0);

/**
 * return the cumulated distribution function
 * @param x the random variable
*/
    double cdf(double x);

/**
 * return the normalized  probability distribution function
 * @param x the random variable
 */
    double pdf(double x);

 /**
     * returns logs of probability density function 
     * @return  Log of probability density function
     * @param   x Point at which is it evaluated
 */
    double logPdf(double x);

/**
 * generate a pseudo random variable according to this distribution
 */
    double generate();

/**
 * generate a pseudo random variable according to exponential distribution
 * @param x0 1/slope
 */
    static double generate(double x0);

  protected :
   
    double x0; 
     
     
     ClassDef(ExponentialDist,1)

};

/**
   * Chi2 distribution
*/

class Chi2Dist : public XeDist {
 
  public :
  
/**
 * Empty constructor for ROOT
 */
    Chi2Dist();
  
/**
 * Regular  constructor 
 * @param ndof number of degrees of freedom
 */
    Chi2Dist(int dof);
   ~Chi2Dist();
 
/**
 * set number of degrees of freedom 
 * @param ndof number of degrees of freedom
 */
    void   setNdof(int ndof);

/**
 * return the cumulated distribution function
 * @param x the random variable
*/
    double cdf(double x);

/**
 * return the normalized  probability distribution function
 * @param x the random variable
 */
    double pdf(double x);

/**
 * generate a pseudo random variable according to this distribution
 */
    double generate();

/**
 * Compute the pdf (static method)
 * @param x random variable
 * @param dof number of degrees of freedom
 */
    static  double pdf(double x,int dof);

/**
 * Compute the integreal above a given chi2
 * @param chi2 tested chi2
 * @param dof number of degrees of freedom
 */
    static  double above(double chi2, int ndof);

/**
 * Compute the integreal below a given chi2
 * @param chi2 tested chi2
 * @param dof number of degrees of freedom
 */
    static  double below(double chi2, int ndof);

/**
 * Compute the integral between two chi2
 * @param chi2_1 lower bound
 * @param chi2_2 upper bound
 * @param dof number of degrees of freedom
 */
    static  double between(double chi2_1, double chi2_2, int ndof);

/**
 * generate a pseduo ranodm variable distributed accoring to a chi2 dist.
 * @param dof number of degrees of freedom
 */
    static  double generate(int dof);

  protected :
   
    int    ndof;
     
     ClassDef(Chi2Dist,1)

};

/**
   * Gaussian distribution
*/

class GaussianDist : public XeDist {
 
  public :
  
  
/**
 * Empty constructor for ROOT
 */
           GaussianDist();

/**
 * Regular constructor 
 * @param mu mean
 * @param sigma sigma
 */
           GaussianDist(double mu, double sigma);
          ~GaussianDist();
 
/**
 * Set the aveage
 */
           void   setMu(double mu);
 
/**
 * Set the standard deviation
 */
           void   setSigma(double sigma);

/**
 * return the cumulated distribution function
 * @param x the random variable
*/
           double cdf(double x);

/**
 * return the normalized  probability distribution function
 * @param x the random variable
 */
           double pdf(double x);

 /**
     * returns logs of probability density function 
     * @return  Log of probability density function
     * @param   x Point at which is it evaluated
 */
           double logPdf(double x);

/**
 * generate a pseudo random variable according to this distribution
 */
           double generate();

/**
 * Compute the fraction of the distribution above n Sigmas
 */
    static double above(double nSigmas);

/**
 * Compute the fraction of the distribution below n Sigmas
 */
    static double below(double nSigmas);

/**
 * Compute the fraction of the distribution between Xmin and Xmax
 * @param Xmin lower bound of the integral
 * @param Xmax upper bound of the integral
 * @param sigma std of the distribution (default is 1.)
 */
    static double inside(double Xmin, double Xmax, double sigma=1.);
    static double logPdf(double x, double mean, double sigma);
    static double generate(double m, double sig);

  protected :
   
    double mu;
    double sigma;
     
     
     ClassDef(GaussianDist,1)

};

/**
   * Poisson distribution
*/

class PoissonDist : public XeDist {

  public :
  
           PoissonDist(double mu);
          ~PoissonDist();
 
           void   setMu(double mu);

/**
 * return the cumulated distribution function
 * @param x the random variable
*/
           double cdf(double x);

/**
 * return the normalized  probability distribution function
 * @param x the random variable
 */
           double pdf(double x);

 /**
     * returns logs of probability density function 
     * @return  Log of probability density function
     * @param   x Point at which is it evaluated
 */
           double logPdf(double x);

/**
 * generate a pseudo random variable according to this distribution
 */
    double generate();

    static double generate(double m);
    static double pdf(int n, double m);
    static double pdf(double n, double m);
    static double logPdf(int n, double m);
    static double logPdf(double n, double m);
    static double aboveEqual(int n, double m);
    static double above(int n, double m);
    static double belowEqual(int n, double m);
    static double below(int n, double m);
    static double between(int n1, int n2, double m);
    static void   fillTable(XeTable& table, double m, double eps);

  protected :

    double mu;


     ClassDef(PoissonDist,1)

};

/**
   * Describing  a cross section, in terms of both cross section and equivalent
   * number of event;
*/
class XeSigma : public XeObject {

  public :

   ~XeSigma();
    XeSigma();

/**
 * Partial constructor
 * @param mass WIMP mass it refers to
 * @param name name (which will be appended with the mass info.)
 */
      XeSigma(double mass,string name="");

/**
 * Full constructor
 * @param mass WIMP mass it refers to
 * @param sigma cross section, in units of cm^2
 * @param events cross section, in units of expected number of signal events
 * @param name name (which will be appended with the mass info.)
 */
    XeSigma(double mass,double sigma, double events,string name="");    

 /**
     * Return the mass 
 */
    double getMass();

 /**
     * Return the cross section, in units of cm^2
 */
    double getSigma();

 /**
     * Return the cross section, in units of expected number of signal events
 */
    double getEvents();

 /**
     * Return the cross section, either in units of cm^2 
     * or expected number of signal events
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
    double getValue(int unit);

/**
  * print the object
  * @param level print level (0-> no print)
*/
    bool   printIt(int level=1);

/**
 * Internal method : set the members
 * @param mass WIMP mass it refers to
 * @param sigma cross section, in units of cm^2
 * @param events cross section, in units of expected number of signal events
 * @param name name (which will be appended with the mass info.)
 */
    void set(double mass,double sigma, double events,string name="");    

/**
 * Internal method : set the values
 * @param sigma cross section, in units of cm^2
 * @param events cross section, in units of expected number of signal events
 */
    void setValues(double sigma, double events);    

/**
 * Internal method : set mass and name
 * @param mass WIMP mass it refers to
 * @param name name (which will be appended with the mass info.)
 */
    void setMassAndName(double mass,string name="");   


  protected :

    double mass;            /*!< Mass for which the cross section is evaluated*/
    double sigma;                                   /*!< Cross section in cm^2*/
    double events;            /*!< Cross section in evaluated number of events*/

/**
 * Return the name 
 * @param mass WIMP mass For for which the cross section is evaluated
 * @param eTitle  title (to which the mass wil be appended)
 */ 
    static string getTheName(double mass,string eTitle);

    ClassDef(XeSigma,1)
};

/**
   * Defining a limit
*/
class XeLimit : public XeStat {

  public :

   ~XeLimit();
    XeLimit();

/**
 * Partial constructor
 * @param mass WIMP mass it refers to
 * @param name name (which will be appended with the mass info.)
 */
 XeLimit(double mass,string name="");

/**
 * Full constructor
 * @param mass WIMP mass it refers to
 * @param estimatedSigma estimated cross section, in units of cm^2
 * @param estimatedEvents estimated cross section, 
 * in units of expected number of signal events
 * @param sigmaLimits estimated cross section, in units of cm^2
 * @param eventsLimits estimated cross section, 
 * in units of expected number of signal events
 * @param name name (which will be appended with the mass info.)
 */
  XeLimit(double mass,double estimatedSigma, double estimatedEvents
             , double sigmaLimit, double eventsLimit,string name="");    

/**
 * Set the mass and the name
 * @param mass WIMP mass it refers to
 * @param name name (which will be appended with the mass info.)
 */
  void setMassAndName(double mass,string name="");    

/**
 * Set the limits
 * @param estimatedSigma estimated cross section, in units of cm^2
 * @param estimatedEvents estimated cross section, 
 * in units of expected number of signal events
 * @param sigmaLimits estimated cross section, in units of cm^2
 * @param eventsLimits estimated cross section, 
 * in units of expected number of signal events
 */
  void setLimits(double estimatedSigma, double estimatedEvents
                ,double sigmaLimit, double eventsLimit);    

 /**
     * Return the mass 
 */
    double getMass();

 /**
     * Return the Limit cross section, in units of cm^2
 */
    double getUpperSigma();

 /**
     * Return the Limit cross section, in units of expected number of signal events
 */
    double getUpperEvents();

 /**
     * Return the Limit cross section, either in units of cm^2 
     * or expected number of signal events
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
    double getUpperLimit(int unit);

 /**
     * Return the Estimated cross section, in units of cm^2
 */
    double getEstimatedSigma();

 /**
     * Return the Estimated cross section, in units of expected number of signal events
 */
    double getEstimatedEvents();

 /**
     * Return the Estimated cross section, either in units of cm^2 
     * or expected number of signal events
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
    double getEstimated(int unit);

 /**
     * Return the Estimated or Limit cross section, either in units of cm^2 
     * or expected number of signal events
     * @param what either ESTIMATED or UPPER_LIMIT
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
    double getValue(int what,int unit);


/**
  * print the object
  * @param level print level (0-> no print)
*/
    bool   printIt(int level=1);

  protected :

/**
 * Return the name 
 * @param mass WIMP mass For for which the limit is evaluated
 * @param eTitle  title (to which the mass wil be appended)
 */ 
    static string getTheName(double mass,string eTitle);

    XeSigma estimated;                            /*!< Estimated cross section*/
    XeSigma limit;                              /*!< Upper limit cross section*/

    ClassDef(XeLimit,1)
};


typedef map<double, XeSigma*>::iterator SigmaIterator;
/**
   * Describing  a set of cross sections
*/
class XeSigmas : public map<double, XeSigma*> , public XeStat  {

   public :

     ~XeSigmas();
/**
 * Simple constructor
 */
      XeSigmas();

/**
 * Constructor
 * @param title name of the XeObject
 */
      XeSigmas(string title);

/**
  * print the object
  * @param level print level (0-> no print)
*/
      bool  printIt(int level=1);

/**
 * Read from a flat file
 * @param fName name of the file (prepended by the results directory name
 */
      void  read(string fName);

/**
 * Add one entry
 * @param sigma  pointer to  XeSigma value
 */
      void  add(XeSigma *sigma);

/**
 * Add one entry
 * @param mass WIMP mass it refers to
 * @param sigma cross section, in units of cm^2
 * @param events cross section, in units of expected number of signal events
 */
      void          add(double mass, double sigma, double events);

/**
 * Create a graph in units of cm2 
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @return a pointer to the newly created XeGraph
 */
      XeGraph*      newSigmaGraph(int plot);

/**
 * Create a graph in its of expected number of signal events
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @return a pointer to the newly created XeGraph
 */
      XeGraph*      newEventsGraph(int plot);

/**    
 * Create a graph in either in units of cm^2 or expected number of signal events
 * @param unit either SIGMA_UNIT or EVENT_UNIT
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @return a pointer to the newly created XeGraph
 */
      XeGraph*      newGraph(int unit,int plot=NONE);

    ClassDef(XeSigmas,1)
 
};


typedef map<double, XeLimit*>::iterator UpperLimitIterator;

/**
   * Describing  a set of cross sections
*/
class XeLimits : public map<double, XeLimit*>,  public XeStat  {

   public :

     ~XeLimits();
/**
 * Simple constructor
 */
      XeLimits();

/**
 * Constructor
 * @param title name of the XeObject
 */
      XeLimits(string title);


/**
  * print the object
  * @param level print level (0-> no print)
*/
      bool printIt(int level=1);

/**
 * Read from a flat file
 * @param fName name of the file (prepended by the results directory name
 */
      void read(string fName);

/**
 * Add one entry
 * @param sigma  the XeLimit value
 */
      void add(XeLimit *limit);

/**
 * Add one entry
 * @param mass WIMP mass it refers to
 * @param estimatedsigma estimated cross section, in units of cm^2
 * @param estimatedevents estimated cross section,
 * in units of expected number of signal events
 * @param sigmaLimit limit cross section, in units of cm^2
 * @param eventsLimit limit cross section,
 * in units of expected number of signal events
 */
 void  add(double mass, double estimatedSigma, double estimatedEvents
          , double sigmaLimit, double eventsLimit );

/**
 * Create a graph of estimated cross section in units of cm2 
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
  XeGraph* newGraphOfEstimatedSigma(int plot=NONE);

/**
 * Create a graph of estimated cross section in units
 * of expected number of signal events
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
  XeGraph* newGraphOfEstimatedEvents(int plot=NONE);

/**
 * Create a graph of estimated cross section in units of cm2 
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
  XeGraph* newGraphOfUpperSigma(int plot=NONE);

/**
 * Create a graph of estimated cross section in units
 * of expected number of signal events
 * @return a pointer to the newly created XeGraph
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
  XeGraph* newGraphOfEventLimit(int plot=NONE);


/**    
 * Create a graph in either in units of cm^2 or expected number
 * of signal events
 * @param mode either ESTIMATED or UPPER_LIMIT
 * @param unit either SIGMA_UNIT or EVENT_UNIT
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @return a pointer to the newly created XeGraph
 */
   XeGraph* newGraph(int mode, int unit, int plot=NONE);

    ClassDef(XeLimits,1)
 
};


/**
 * Class describing one set of sensitivity: median, +/- 1 and +/- 2 sigma upper limit
*/
class XeSensitivity :  public XeStat {

   public :

    ~XeSensitivity();
/**
 * Empty constuctor for ROOT
 */ 
     XeSensitivity();

/**
 * Constructor without the limits
 * @param mass  mass at which the sensitivity is measured
 * @param signalPerCm2  Number of expected events per cm^2
 * @param eTitle  Name of the sensitibity (the mass will be added!)
 */
     XeSensitivity(double mass,double signalPerCm2,string eTitle="");

/**
 * Constructor with the limits
 * @param mass  mass at which the sensitivity is measured
 * @param signalPerCm2  Number of expected events per cm^2
 * @param limits array of upper limits for -2,1,0,+1,+2 sigmas
 * @param eTitle  Name of the sensitibity (the mass will be added!)
 */
     XeSensitivity(double mass,double signalPerCm2, double * limits
                  ,string eTitle="");

/**
  * print the object
  * @param level print level (0-> no print)
*/
     bool    printIt(int level=1);

/**
 * Set one limit
 * @param more MINUS_TWO_SIGMAS , MINUS_ONE_SIGMA , MEDIAN 
 *        , PLUS_ONE_SIGMA , PLUS_TWO_SIGMAS
 * @param limit Upper cross section incm^2
 */
     void    setLimit(int mode, double limit);

/**
 * Set all limits
 * @param limits pointer limits for for -2,1,0,+1,+2 sigmas
 */
     void    setLimits(double* limits);

 /**
     * get the sensitivity
     * @param which MINUS_TWO_SIGMAS , MINUS_ONE_SIGMA , MEDIAN 
     *        , PLUS_ONE_SIGMA , PLUS_TWO_SIGMAS
     * @param unit wither SIGMA_UNIT or EVENT_UNIT
 */
     double  getLimit(int which,int unit=SIGMA_UNIT);

/**
 * get the limits in SIGMA_UNIT (cm^2)
 */
     double* getLimits();

/**
 * Return the mass for which the limits are computed
 */   
     double  getMass();
  
/**
 * Return the conversion factor cm^2->events (SIGMA_UNIT -> EVENT_UNIT)
 */ 
    double  getSignalPerCm2();

   protected :

/**
 * Compute the name , appending the WIMP mass
 * @param  mass mass
 * @param eTitle simple title
 */
     static string getTheName(double mass,string eTitle);

     double mass;             /*!< Mass at which the sensitivity is evaluated */
     double signalPerCm2;                 /*!<conversion factor cm^2->events  */
     double limits[N_SENSITIVITY_MODES];      /*!< limits in SIGMA_UNIT (cm^2)*/

     ClassDef(XeSensitivity,1)
};

typedef map<double, XeSensitivity*>::iterator SensitivityIterator;

/**
   * Describing the median, +/- 1 and +/- 2 sigma upper limits bands
*/
class SensitivityBands : public XeStat
           , public map<double, XeSensitivity*> ,public XeStylized {

   public :

     ~SensitivityBands();

/**
 * Empty constructor for ROOT
 */
      SensitivityBands();

/**
 * Regular constructor
 * @param runName Name of sensitivity bands
 */
      SensitivityBands(string runName);

/**
 * Return line style
 * @param which MINUS_TWO_SIGMAS , MINUS_ONE_SIGMA , MEDIAN 
 *        , PLUS_ONE_SIGMA , PLUS_TWO_SIGMAS
 */
      static int lineStyle(int which);

/**
 * Return number of sigmas
 * @param which MINUS_TWO_SIGMAS , MINUS_ONE_SIGMA , MEDIAN 
 *        , PLUS_ONE_SIGMA , PLUS_TWO_SIGMAS
 */
      static double numberOfSigmas(int which);

/**
 * Return band name
 * @param which MINUS_TWO_SIGMAS , MINUS_ONE_SIGMA , MEDIAN 
 *        , PLUS_ONE_SIGMA , PLUS_TWO_SIGMAS
 */
      static string bandName(int which);

/**
 * Read for a flat file
 */
      void read(string fName);

 /**
     * draw the sensitivity bands in either sigma or events limits
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
      void   draw(string options="C",int unit=SIGMA_UNIT);

/**
 * Return max X
 */
      double getMaxX();

/**
 * Return min X
 */
      double getMinX();

/**
 * Return max Y
 */
      double getMaxY();

/**
 * Return min Y
 */
      double getMinY();

/**
 * Return min Y when different from 0, for plots in log scale
 */
      double getMinYNotZero();

/**
  * print the object
  * @param level print level (0-> no print)
*/
      bool printIt(int level=1);

/**
 * Add a sensitivity entry
 * @param sensitivity Pointer to sensitivity
 */
      void add(XeSensitivity *sensitivity);

/**
 * Add an entry
 * @param mass WIMP mass for which it was evaluated
 * @param signalPerCm2 conversion factor cm^2->events
 * @param limits pointer to table of limits
 */
      void add(double mass, double signalPerCm2, double* limits);

 /**
     * get the graph of one sensitivity curve
     * @param which which band MINUS_RWO_SIGMAS,MINUS_ONE_SIGMA,CENTRAL,
     *                         PLUS_ONE_SIGMA,PLUS_TWO_SIGMAS
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
      XeGraph*      newGraph(int which,int unit=SIGMA_UNIT);

 /**
     * get the multigraph of all sensitivity curves
     * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
      XeMultiGraph* newMultiGraph(int unit=SIGMA_UNIT);

/**
 * Draw a part of the brazilian flag
 * @param low index of the lower part of the part
 * @param high index of the upper part of the part
 * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
      TPolyLine*    newPolygon(int low, int high, int unit);

      ClassDef(SensitivityBands,1)
};

/**
   * Data set : a collection of nCol colmuns *  nEntries entries
*/
class DataSet : public XeStat   {

   public :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
/** 
 * Regular constructor
 * @param name  DataSet name
 * @param nCol  number of columns
 */
     DataSet(string name,int nCol);

/** 
 * Short constructor
 * @param nCol  number of columns
 */
     DataSet(int nCol);


/**
 * Empty constructor for ROOT
 */
     DataSet();
    ~DataSet();

/**
 * Add an entry
 * @param x pointer to the ncolumns values
 */
     void             addEntry(double* x);

/**
  * print the object
  * @param level print level (0-> no print)
*/
     bool             printIt(int level=1);

/**
 * Clear everything
 */
     void             clear();
   
/**
 * Get an entry
 * @param entry  entry number, starting at 0
 * @param values pointer to receive the array of values
 */
     void             getEntry(int entry,double* values); 

/** 
 * Return number of columns ( values per entry)
 */
     int              getNColumns();
 
/**
 * Return number of events, or entries  
 */
     int              getNEvents();
   
/** 
 * Return one value in a given entry
 * @param entry entry number, starting at zero
 * @param col value index, starting at zero
 */
     double           getValue(int entry,int col); 

/**
 * Return maximum value in a column
 * @param col value index, starting at zero
 */
     double           getMax(int col);

/**
 * Return minimum value in a column
 * @param col value index, starting at zero
 */
     double           getMin(int col);

/**
 * Return one column as an array
 * @param col value index, starting at zero
 * @return pointer to first element
 */
     double*          getVector(int col);

/**
 * Return one column as a vector
 * @param col value index, starting at zero
 * @return pointer to vector
 */
     vector<double>*  getColumn(int col);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
/**
 * Virtual method to be called when smthg has changed
 */
     virtual bool     update();

/** 
 * Initialization
 * @param nCol  number of data columns
 */
  void setColumns(int nCol);

   protected :
   
     vector<vector<double>* > entries;  /*!< vector of columns*/
     int                      nColumns;  /*!< Number of columns*/
     int                      nEvents;  /*!< Number of entries=events*/
     
      ClassDef(DataSet,1)
};

/**
   * Simulated Data set 
*/

class SimulatedDataSet : public DataSet {

   public :
    
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     SimulatedDataSet(string name, int nDim);
     SimulatedDataSet();

     virtual void   initialize();
     virtual void   generateSignal(double *x)=0;
     virtual void   generateBackground(double *x)=0;
             void   simulate(double muSignal, double muBackground);
             void   generate(int nSignal, int nBackground);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    ~SimulatedDataSet();

   protected :

     vector<double> transient;
     
      ClassDef(SimulatedDataSet,1)

};

/**
   * One dimensional Simulated Data set 
*/

class OneDimSimulatedDataSet : public SimulatedDataSet {

  public :
   
   /* -------------------------------------------------------------
    *                     Basic methods 
    * ------------------------------------------------------------*/
                   
    OneDimSimulatedDataSet();

    void    generateSignal(double *x);
    void    generateBackground(double *x);
    void    setSignal(XeDist* s);
    void    setBackground(XeDist* b);

   /* -------------------------------------------------------------
    *                     Advanced methods 
    * ------------------------------------------------------------*/
                    
    XeDist* getSignal();
    XeDist* getBackground();    

   /* -------------------------------------------------------------
    *                Internal methods (not for user)
    * ------------------------------------------------------------*/
                   
   ~OneDimSimulatedDataSet();
    void    initialize();

  protected :

    XeDist* signal;
    XeDist* background;    
 
      ClassDef(OneDimSimulatedDataSet,1)

};


/**
   * For MC studies: exopnential over flat background
*/

class SimulatedExponentAndBackground : public OneDimSimulatedDataSet {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    SimulatedExponentAndBackground(double x0=DEFAULT_SIMULATED_X0
                                  ,double xMaxB=DEFAULT_SIMULATED_XMAX_B);
     
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    void   setX0(double x0);
    void   setXmaxB(double xm);
                    
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

          ~SimulatedExponentAndBackground();

  protected :

    ExponentialDist *exp;
    UniformDist     *uni;    
   
      ClassDef(SimulatedExponentAndBackground,1)


};

/**
   * For MC studies: gaussian over flat background
*/

class SimulatedGaussianAndBackground : public OneDimSimulatedDataSet {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    SimulatedGaussianAndBackground( double mu=DEFAULT_SIMULATED_MU_GAUSS
                                  , double sigma=DEFAULT_SIMULATED_SIGMA
                                  , double xMaxB=DEFAULT_SIMULATED_XMAX_B
                                  ) ;

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    void   setXmaxB(double xm);
    void   setMu(double mu);
    void   setSigma(double sigma);

                    
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

          ~SimulatedGaussianAndBackground();

  protected :

    GaussianDist *gau;
    UniformDist  *uni;    

      ClassDef(SimulatedGaussianAndBackground,1)
};

/**
   * The basic class for CLs and PL : a signal over background.
*/


class SignalAndBackground :  virtual public XeStat {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
  public : 

    virtual ~SignalAndBackground();
             SignalAndBackground(double sigToE=1., double bkg=0., double dbkg=0.
                        ,double eff=1., double deff=0.);
            // sigToE: how many events per units of sigma and exposure
            // bkg   : how many background events per unit of exposure
            // eff   : efficiency affecting signal only

    void     combine(SignalAndBackground *sb1, SignalAndBackground *sb2);
    void     combine(vector<SignalAndBackground*> &vsb);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    void     printSBParameters();
    void     setBackground(double b=0.);
    void     setDBackground(double db=0.);
    void     setExposure(double exposure=1.);
    void     setSigToEvents(double v=1.);
    void     setNObserved(int observed=0);
    void     setNObservedAsBackground();
    void     setEfficiency(double e=1.);
    void     setDEfficiency(double de=0.);
    void     fillPoissonTable(XeTable& table, double eps=1.e-5, double sig=0.);
    int      getNObserved();
    double   getBackground();
    double   getDBackground();
    double   getSigToEvents();
    double   getExposure();
    double   getEfficiency();
    double   getDEfficiency();
    double   nSignalEvents(double sigma);
    double   nBackgroundEvents();
    double   nTotalEvents(double sigma);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    void     reset();
    double   probability(int observed=AUTO,double sig=0.);

  protected:

    
    int      nObserved;
    double   sigToEvents;
    double   exposure;
    double   efficiency;
    double   dEfficiency;
    double   background;
    double   dBackground;

      ClassDef(SignalAndBackground,1)
};


/**
   * Simple Confidence Interval.
   * This is a virtual class
*/


class CI : virtual public solvable, public XeStat {  

  public: 

                CI(int mode=UNDEFINED_INT);
    virtual    ~CI();
    bool        isInside(double x);
    bool        withLowerLimit();
    bool        withUpperLimit();
    double      getLowerLimit();
    double      getUpperLimit();
    void        applyCLs(bool doIt);
    bool        isCLsApplied();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

  
  protected:

    int         mode;
    bool        withUpper;
    bool        withLower;
    bool        doApplyCLs;
    double      LowerLimit;
    double      UpperLimit;
    static bool withCLs(int mode);
    static bool withLowerLimit(int mode);
    static bool withUpperLimit(int mode);
    static bool isSingleLowerLimit(int mode);
    static bool isSingleUpperLimit(int mode);
};

/**
   * Poisson Confidence Interval
*/


class PoissonCI : public CI {



    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/

    public :

               PoissonCI(SignalAndBackground *sb,int mode=UNDEFINED_INT);
               PoissonCI(double sigToE=1.,double bkg=0.,int mode=UNDEFINED_INT);
      virtual ~PoissonCI();
      void     computeLimits(int nObserved, double CL=DEFAULT_CL);
      void     computeLimits(double CL=DEFAULT_CL);
      bool     expectsBackground();
      double   computeUpperSigma(double CL=DEFAULT_CL);
      double   coverage(double sigma,double CL=DEFAULT_CL);

      static   bool   expectsBackground(int mode);
      static   double upperLimit(int m, int nObs, double bkg, double CL);
      static   double lowerLimit(int m, int nObs, double bkg, double CL);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

      double   getValue(double sigma);
      double   oneLimit(int currentMode, double CL);
      void     setBackground(double b=0.);
      void     setExposure(double exposure=1.);
      void     setSigToEvents(double v=1.);
      void     setNObserved(int observed);
      void     setNObservedAsBackground();
      void     setEfficiency(double e=1.);
      void     setSignalAndBackground(SignalAndBackground *sb);

      SignalAndBackground* getSignalAndBackground();

    protected :

      bool                 deleteSB;
      int                  currentMode;
      SignalAndBackground* SB;
      ClassDef(PoissonCI,1)
};


/**
   * Rejection modes (internal and technical)
*/

class Rejection : virtual  public XeStat {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public: 

    Rejection();
   ~Rejection();

    static string   getModeName(double e);
    static string   getModeName(int mode);
    static int      getMode(double eff);
    static double   getRejection(int mode);

  protected :

    static double rejections[N_REJECTION_MODES];
    ClassDef(Rejection,1)

};



/**
   * Basic class for calculating P values.
   * This is a virtual class.
   * In this class, cross sections are the real ones (1E-45 or so); 
   * No scaling is performed.
*/
class PValue : virtual public XeStat {

 
  public :

    /* ------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
/**
 * void constructor for ROOT
 */
  PValue(); 

 /**
 * Regular constructor
 * @param name      name of the p-value calculator
 * @param analMode  either PL_ANALYSIS,CUTS_ANALYSIS,NO_ANALYSIS
 */
    PValue(string name, int analMode=NO_ANALYSIS);

 /**
     * compute p-value for a given cross section
     * @return p-value for a given cross section
     * @param  sigma true value of the cross section
 */
   virtual double  pValueS(double sigma)=0;   

 /**
     * Virtual method returning exclusion q value for a tested sigma
     * @return  q(sigma)
     * @param  sigma     cross sections being testes
 */
   virtual double qExclusion(double sigma);

 /**
     * Virtual method returning current estimated sigma
     * @return  current value of estimated sigma
 */
   virtual double getSigmaHut(); //ALE FIXME this has to be reomeved, actually the whole class has to be removed

/**
     * Typical number of expected signal events per cm2 of cross sections
 */
   virtual double  nSignalPerCm2()=0;

 /**
     * update after a change in conditions (mainly sigma test) 
     * @return bool was update OK ?
 */
   virtual bool    update()=0;

 /**
     * print All flags and parameters
 */
   virtual void    printFlagsAndParameters()=0;

 /**
     * retrieve parameters at end of exclusion finder, and print them if needed
 */
   virtual void    exclusionComputed();

 /**
    * Simulate an outcome
    * @param  sigma  cross section simulated
 */
   virtual bool    simulate(double sigma=0.);

/**
 * Virutal method: pass the wimp mass
 */
   virtual void    setWimpMass(double mass);

/** 
 * Virtual method: update cross section, in cm^2, corresponding to 1 event
 */
   virtual void    updateSigmaUnit();
  
/**
 * virtual method: estimate the cross section (actually for P.L.)
 */
   virtual void    estimateCrossSection();
  


/**
 * Virtual method saying is CLs can be used
 * @return DO_NOT_FORCE, FORCE_TRUE or FORCE_FALSE
 */
  virtual int     forceCLs();
           
 /**
     * @param  scale   sigma unit to be passed to Minuit
 */
   void            setSigmaUnit(double scale);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

/**
 * Intitialize at instantiation.
 * Method used as a workaround to diamond hierarchy problems
 * @param mode either NO_ANALYSIS, CUTS_ANALYSIS or PL_ANALYSIS
 */
   void            initializeIt(int mode=NO_ANALYSIS);

   virtual        ~PValue();

/**
 * Check that the mode is compatible with the current PL or Cut Based anaylsis
 */
   bool            isAnalysisModeOK();

/**
 * Is all initialized ok?
 */   
   bool            isInitializedOK();

/**
 * Is there an error ?
 */   
   bool            isInError();

/**
 * Checks (once) that all is iniatlized ok
 */
   bool            initialize();  
  
/**
 * Commpute p-value for the background
 */
   double          pValueB();

/**
 * Return upper limit of number of events in the fit
 */
   virtual double  eventsUpperLimit();

/**
 * Return lower limit of number of events in the fit
 */
   virtual double  eventsLowerLimit();

/**
 * Set upper limit of number of events in the fit
 */
   virtual void    setEventsUpperLimit(double l);

/**
 * Set lower limit of number of events in the fit
 */
   virtual void    setEventsLowerLimit(double l);
                    
  protected :

/**
 * check that pValue is OK
 */
   virtual bool    checkPValue()=0;

   int             requestedAnalysisMode;  /*!<Requested mode (PL or cuts)*/
   bool            initializedOK;  /*!<Is all ok?*/
   bool            inError;  /*!<Has an error been detected?*/
   double          sigmaUnit;  /*!<Sigma corresponding to 1 events*/
   double          lowerLimit;  /*!<Lower number of events in fit*/
   double          upperLimit;  /*!<Upper number of events in fit*/

   ClassDef(PValue,1)
};


/**
   * Pvalue for poisson distrbution.
   * It must inherit virtually from  SignalAndBackground
*/
class PVcountingSB :public PValue, public SignalAndBackground {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
  public :

   ~PVcountingSB();
    PVcountingSB(double sigToE=1., double bkg=0. ,double eff=1.);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

 /**
     * print All flags and parameters
 */
    void    printFlagsAndParameters();

 /**
     * print All flags and parameters
 */
    void    setEventsUpperLimit(double lim);
    bool    update();
    bool    checkPValue();
    double  pValueS(double sigma);   

/**
 * return expected number of signal events for 1 cm2 cross section
*/
    double  nSignalPerCm2();
    double  eventsUpperLimit();
    double  getSigmaHut();   
 
    ClassDef(PVcountingSB,1)

};


/**
   * Exclusion a la Yellin
*/


class YellinPValue : public PValue {
 
  /* virtual class */

  public :


    /* ------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

    virtual   ~YellinPValue();
               YellinPValue(DataSet* dataSet, XeDist* signalDistribution
                           , double sToE=1., int col=0);
               YellinPValue(OneDimSimulatedDataSet* dataSet,double sToE=1.);
    double     pValueS(double sigma);   

/**
 * return expected number of signal events for 1 cm2 cross section
*/
    double     nSignalPerCm2();      
    double     eventsUpperLimit();
    bool       update();
    bool       checkPValue();

 /**
     * print All flags and parameters
 */
    void       printFlagsAndParameters();
    int        forceCLs();
   
    DataSet*   getDataSet();

  protected :

    int        column;
    double     sigToEvents;
    double     maxGapForOne;
    DataSet*   dataSet;
    XeDist*    signal;

    ClassDef(YellinPValue,1)
     
};


/**
   * Basic machinery for exclusion
*/

class Exclusion : public XeStat, public solvable {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    

   ~Exclusion();
/**
 * Empty constructor for ROOT
 */
    Exclusion();

 /**
     * Constructor with one p-value (PValue) calculator
     * @param pValue  a p-value (PValue) calculator
 */
    Exclusion(PValue *pValue);

 /**
     * Constrictor with several p-value (PValue) calculators
     * @param pValueS vectors of p-values calculator; 
     * the Fishcher correction will be used
 */
    Exclusion(PValue** pValueS, int nP);
  

 /**
     * Compute cross section
     * @return pointer to a newly created XeSigma
 */
    XeSigma* computeCrossSection();
 
 /**
     * Compute sigma upper limit
     * @param  cl  Confidence Level
 */
    void computeLimit(double cl=DEFAULT_CL);

 /**
     * Compute and return sigma upper limit
     * @return limit in cross section
     * @param  cl  Confidence Level
 */
    double computeUpperSigma(double cl=DEFAULT_CL);

/**
 * Return lastly computed upper sigma
 */
   double getUpperSigma(); 

/**
 * Return lastly computed estimated sigma
 */
   double getEstimatedSigma(); 

/**
 * Return lastly computed upper number of events
 */
   double getUpperEvents(); 

/**
 * Return lastly computed estimated number of events
 */
   double getEstimatedEvents(); 

/**
 * Return lastly computed XeLimit object
 */
   XeLimit* getLimit(); 



 /**
    * Simulate an outcome
    * @param  sigma  cross section simulated
 */
    void simulate(double sigma);

 /**
     * @return  q(sigma) according to PL recipe
     * @param  sigma     cross sections being testes
 */
    double qExclusion(double sig);

 /**
     * Simulate n runs and return the n upperLimits 
     * @param n   Number of simulations
     * @return sigmas the list of upper sigmas (actually XeValues)
     * @return events the list of upper events (actually XeValues)
     * @param  cl  Confidence Level
 */
    void simulateNUpperLimits(int n, XeValues *sigmas=NULL
                             , XeValues *events=NULL,double cl=DEFAULT_CL);


 /**
     * Simulate n runs and return the n qValues 
     * @param n   Number of simulations
     * @return qs list of q values        (actually XeValues)
     * @return sHats list of sigma values after fit (actually XeValues)
     * @return eHats list of Nevents values after fit (actually XeValues)
     * @param sigSim Cross section for the simulation
     * @param sigTest Cross section to be tested
 */
  void simulateNQValues(int n, XeValues *qs=NULL, XeValues* sHats=NULL
                  ,XeValues* eHats=NULL,double sigSim=0.,double sigTest=0.);

 /**
 * @return  an XeGraph of Qexclusion for a sigma range
 * @param sr Sigma range
 * @param unit either SIGMA_UNIT or EVENT_UNIT
 */
    XeGraph* newGraphOfQExclusion(XeRange* sr, int unit=EVENT_UNIT);

 /**
     * Return one point in the sensitivity bands
     * @return  pointer to a newly created sensitivity point
     * @param mass WIMP mass
     * @param n   Number of simulations
     * @param  cl  Confidence Level
 */
    XeSensitivity* newSensitivity(double mass, int n, double cl=DEFAULT_CL);


 /**
     * Return one point in the sensitivity bands
     * @return  pointer to a newly created sensivitity band
     * @param n   Number of simulations
     * @param  cl  Confidence Level
 */
    SensitivityBands* newSensitivityBands(int n, XeRange* mr=NULL
                                        ,double cl=DEFAULT_CL);

 /**
     * Return estimated and limit sigmas, in terns of both sigma and events
     * @return  pointer to a newly created XeLimits
     * @param mr  mass range (NULL means default ine)
     * @param  cl  Confidence Level
 */
     XeLimits* newUpperLimits(XeRange* mr=NULL,double cl=DEFAULT_CL);

 /**
 * Return graph of p-values for various test cross sectoins 
 * @return the XeGraph
 * @param sr an XeRange describing the Sigma range
 * @param plot either NONE, AUTO, LINEAR, or LOG
 */
    XeGraph* newGraphOfPValue(XeRange* sr,int plot=NONE);

  /**
 * Return the XeGraph of upper sigma 
 * @return  An XeGraph of upper limit for a given mass range
 * @param   mr Mass Range for the WIMP particle
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @param cl Confidence level
 */                   
    XeGraph* newGraphOfUpperSigma(XeRange* mr, int plot=NONE
                                 ,double cl=DEFAULT_CL);

  /**
 * Return the XeGraph of upper events 
 * @return  An XeGraph of upper limit for a given mass range
 * @param   mr Mass Range for the WIMP particle
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @param cl Confidence level
 */                   
    XeGraph* newGraphOfUpperEvents(XeRange* mr,int plot=NONE
                                  ,double cl=DEFAULT_CL);


  /**
 * Return the XeGraph of upper limit(sigma or event )
 * @return  An XeGraph of upper limit for a given mass range
 * @param   unit  either SIGMA_UNIT or EVENT_UNIT
 * @param   mr Mass Range for the WIMP particle
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @param cl Confidence level
 */                   
    XeGraph* newGraphOfUpperLimit(int unit,XeRange* mr,int plot=NONE
                                , double cl=DEFAULT_CL);

 /**
 * Reset everything
 */
    void     reset();

 /**
     * set current wimp mass
     * @param mass requested WIMP mass
 */
    void     setWimpMass(double mass);

 /**
     * return current WIMP mass
 */
    double   getWimpMass();

 /**
     * Set the scales and limits, to deal with numbers of order unit
     * @param sigmaUnit  Sigma corresponding to 1 expected event after cuts
     * @param eMin       Lower limit of events sigma_up search (usually 0)
     * @param eMax       Upper limit of events sigma_up search (usually 20-30)
 */
    void     setSigmaUnitAndLimits(double sigmaUnit, double eMin, double eMax);

 /**
     * Set the "CLs" correction flag
     * @param doIt apply or don't apply the CLs correction
 */
    void     applyCLs(bool doIt);

 /**
     * Set the "CLs" correction flag and remember the current setting
     * @param doIt apply or don't apply the CLs correction
 */
    void     applyCLsAndSave(bool doIt);

 /**
     * Save the current setting of the "CLs" correction flag 
*/
    void     saveCLs();

 /**
     * Restore the current setting if the "CLs" correction flag 
*/
    void     restoreCLs();

 /**
     * Set the Fisher correction flag when more than one p-value is checked
     * @param doIt apply or don't apply the Fisher correction
 */
    void     applyFisherCorrection(bool doIt);

/**
 * Is Cls correction to be applied ?
 */
    bool     isCLsApplied();
    
 /**
     * For rescaling, gives a typical upper limit for an observed number of events
     * @param nev number of events
 */
    static  double  typicalUpperLimit(int nev);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

 /**
     * Preparatory computations for the exclusion limits.  
     * To be called each time a new simulated run is created
     * @return  bool went ok or not
 */
    bool     prepareForLimit();
               
 /**
     * Update stuff when a parameter (e.g. test mass) changed
 */
    bool     update();
               
 /**
     * Make sure everything is OK, in place
 */
    bool     initialize();

 /**
     * Compute (CLs corrected) pValue
     * @return p-value between 0 and 1
     * @param sigma test cross section
     * @param usedSaved use the already computed and save p-value for bkg only
 */
    double   correctedPValue(double sigma,bool useSaved=false);

 /**
     * Technical wrapped for equation solving
     * @param sigma  test cross section
 */
    double   getValue(double sigma);
   
/**
 * Return the limit object
 * @return pointer to the limit
 */

  protected :
   
    int      nPValues;                      /*!<Number of p values calculators*/
    bool     initializedOK;                     /*!<Was everything analyzed ok*/
    bool     doApplyCLs;                          /*!<Do we want to apply CLs?*/
    bool     savedApplyCLs;                 /*!<Temporary save of 'doApplyCLs'*/
    bool     fisherCorrection;      /*!<Apply Fisher correction if >1 p-value?*/
    double   wimpMass;                                  /*!< Current WIMP mass*/
    double   savedPBkg;            /*!< P value at background used for the CLS*/
    double   signalPerCm2;      /*!<Translation from cross section to n Events*/
    double   eventMin;       /*!<Lower bound of number of events for exclusion*/
    double   eventMax;       /*!<Upper bound of number of events for exclusion*/
    XeLimit  limit;                                          /*<!Current limit*/
    vector<PValue*>  pValues;                   /*!<P-value calculators in use*/

/**
 * Add p-Values to the exclusion
 * @param pV pointer to array to pointers of PValue
 * @param l number of pValues
 */
    void     addPValues(PValue** pV,int l);

/**
 * Add one p-Value to the exclusion
 * @param pV pointer to  PValue
 * @param l number of pValues
 */
    void     addPValue(PValue* pV);

/**
 * Compute p-value for background only hypothesis
 */
    double   computePB();

    ClassDef(Exclusion,1)
     
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

    ClassDef(LKParameter,1)
     
} ;

/**
   * Predifined parameter: sigma
*/
class SigmaParameter : public LKParameter {
  public :
    SigmaParameter();
   ~SigmaParameter();

    ClassDef(SigmaParameter,1)
     
} ;

/**
   * Predefined parameter: t-value for systematic background uncertainty
*/
class TSystBkgParameter : public LKParameter {
  public :
    TSystBkgParameter(int run);
   ~TSystBkgParameter();

  static string getTheName(int run);

    ClassDef(TSystBkgParameter,1)
} ;

class TGaussParameter : public LKParameter {
  public :
    TGaussParameter(int id, int run);
   ~TGaussParameter();

    void   setUncertainty(double unc) { uncertainty = unc;};
    double getUncertainty()           { return uncertainty; };

    static string getTheName(int b, int run);

    double uncertainty;

    ClassDef(TGaussParameter,1)
} ;

/**
   * Predefined parameter: t-value for statistical background uncertainty
*/
class TStatBkgParameter : public LKParameter {

  public :
    TStatBkgParameter();
    TStatBkgParameter(int band, int run);
   ~TStatBkgParameter();

    ClassDef(TStatBkgParameter,1)

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

    ClassDef(TLEffParameter,1)

} ;

/**
   * Predefined parameter: t-Qy
*/
class TQyParameter : public LKParameter {
  public :
    TQyParameter();
   ~TQyParameter();

    ClassDef(TQyParameter,1)

} ;


/**
   * Predefined parameter: t-efficiency
*/
class TEfficiencyParameter : public LKParameter {
  public :
    TEfficiencyParameter();
   ~TEfficiencyParameter();

    ClassDef(TEfficiencyParameter,1)

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
     * return the shape LL that a set of values comes from a parent distribution
     * @return Log Likelihood
     * @param values  values
     * @param  dist parent distribution
 */
     static   double shapeLikelihood(vector<double>* values, XeDist* dist);

 /**
     * return the shape LL that a set of values comes from a 
     * mixture of parent distributions.
     * @return Log Likelihood
     * @param  values  values
     * @param  dists parent distributions
     * @param  weights relative weight of distrubtions (needs not be normalized to 1)
 */
     static   double shapeLikelihood(vector<double>* values
                                    ,vector<XeDist*>& dists
                                    ,vector<double>& weights);
 
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

    ClassDef(Likelihood,1)

} ;

typedef map<int,LKParameter*>::iterator ParameterIterator;

#define TRAVERSE_PARAMETERS(it) \
  for(ParameterIterator it=parameters.begin(); it!=parameters.end(); it++) 


/**
   * Likelihood from a data set.
   * This is a virtual class
*/

class LikelihoodFromDataSet :  public Likelihood {

  /* virtual class */

  public : 
                    
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    LikelihoodFromDataSet(DataSet* data);
    virtual        ~LikelihoodFromDataSet();
    virtual double  computeIndividualLog(double *values)=0;
            double  computeTheLogLikelihood();
    DataSet*        getDataSet(); 
   
  protected :

    DataSet*        data;
    vector<double>  values;
    static string   getTheName(DataSet* data);

    ClassDef(LikelihoodFromDataSet,1)

} ;


/**
   *  Profile likelihood calculator
*/

class ProfileLikelihood : virtual public Likelihood, virtual public PValue {
    
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

 /**
     * @return  q(sigma) according to PL recipe
     * @param  sigma     cross sections being testes
 */
     double   qExclusion(double sigma);
    
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
  
     void     updateSigmaUnit();
 
 /**
     * Compute the p-value for a given cross section
     * @return  double pValue for a given cross section
     * @param   sigma  Cross section being tested
 */
    double   pValueS(double sigma);
    virtual ~ProfileLikelihood();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

 /** 
 * check that p-values are all ok for limit computations
 */
    virtual  bool  checkPValue();

 /**
     * print All flags and parameters
 */
    void    printFlagsAndParameters();


 /**
     *  estime cross section and save corresponding likelihood 
 */
     void    estimateCrossSection();
 

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

    ClassDef(ProfileLikelihood,1)

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
    void               setWimpMass(double);
    double 	       getWimpMass();

    void setData(int dataType);

    void generateAsimov(double mu_prime) ;

    void generateToyDataset(double seed, double mu_prime);

    double getSignalDefaultNorm();
    double getSignalMultiplier();
    void   setSignalMultiplier(double val); 
/**
 * return expected number of signal events for 1 cm2 cross section
*/
    double             nSignalPerCm2();      


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
             
    bool    update();
    

  protected:

    bool    checkPValue();   
    map<int,ProfileLikelihood* > exps;
    int                          nCommon;
    double                       sigToEvents;

    ClassDef(CombinedProfileLikelihood,1)

};

typedef map<int,ProfileLikelihood*>::iterator expIterator;

#define TRAVERSE_EXPERIMENTS(it) \
 for(expIterator it=exps.begin(); it!=exps.end(); it++) 


/**
   * Simple profile likelihood based on Poisson S+B
*/


class PLcountingSB : virtual public SignalAndBackground, 
                     virtual public ProfileLikelihood{

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
  public :

   ~PLcountingSB();
    PLcountingSB(double sigToE=1., double bkg=0., double dbkg=0.
                ,double eff=1., double deff=0.);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    double  computeTheLogLikelihood();

 /**
     * print All flags and parameters
 */
    void    printFlagsAndParameters();
    void    setWimpMass(double mass);
    bool    update();

/**
 * return expected number of signal events for 1 cm2 cross section
*/
    double  nSignalPerCm2();
    double  eventsUpperLimit();


    ClassDef(PLcountingSB,1)

};
#endif



