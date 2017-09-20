#ifndef XeMath_h
#define XeMath_h
#include <set>
#include "Math/Functor.h"
#include "Math/BrentRootFinder.h"
#include "XeCore.h"

/**
 * Base class for math operations
 */

class XeMath : virtual public XeCore {
  
public:


/**
 * Destructor
 */
 ~XeMath();

/**
 * Void constructor for ROOT
 */
  XeMath();

/**
 * equalize a vector, set all its value equal, conserving the sum
 * @param v pointer to array
 * @param n number of elements
 */
  static bool    equalizeVector(double *v, int n);

/**
 * equalize a vector, set all its value equal, conserving the sum
 * @param v vector
 */
  static bool    equalizeVector(vector<double>&v);

/**
 * normalize a vector so that its sum is 1.
 * @param vin pointer to input array
 * @param n number of elements
 * @param vout pointer to output array; if NULL, will be replaced in place
 */
  static bool    normalizeVector(double *vin, int n, double *vout=NULL);

/**
 * normalize a vector so that its sum is 1.
 * @param v vector (both input and output)
 */
  static bool    normalizeVector(vector<double>&v);

/**
 * normalize a vector so that its sum is 1.
 * @param vin input vector
 * @param vout output vector
 */
  static bool    normalizeVector(vector<double>&vin,vector<double>&vout);

/** 
 * compute the sum of a vector
 * @param v pointer to input array
 * @param n number of elements
 */
  static double  vectorSum(double *v, int n);

/**
 * print content of a vector of double with optional header
 */
  static bool    printVector(double *v, int n, string header="");

/**
 * print the log10 of content of a vector of double with optional header
 */
  static bool    printLogVector(double *v, int n, string header="");

/**
 * print content of a vector of double with optional header
 */
  static bool    printVector(vector<double>&v, string header="");

/**
 * print content of a vector of int with optional header
 */
  static bool    printVector(int *v, int n, string header="");

/**
 * print content of a vector of int with optional header
 */
  static bool    printVector(vector<int>&v, string header="");


/** 
 * compute the sum of a vector
 * @param v vector 
 */
  static double  vectorSum(vector<double>&v);
 
/**
 * Return the index of  value in a interval
 * @param v Value to be checked
 * @param vMin Lower edge of the interval
 * @param vMax Upper edge of the interval
 * @param nInterval number of intervals
 * @param log are we in log scale?
 */
  static int    indexInterval(double v, double vMin, double vMax
                             ,int nInter, bool log=false);

/**
 * Compute polynomail in a efficient way
 * @param x value at which it's evaluated
 * @param n degree of a polynomial
 * @param coef array of n+1 coefficients
 */
  static double computePolynomial(double x, int n, double *coef);

 /**
     * return name of a t-value 
     * @param t "t-value"
 */
    static string tValueName(double t);

  ClassDef(XeMath,1)
};

/**
   * Tolerance checker
*/

class XeTolerance: public XeMath {

  public :


/**
 * Destructor
 */
   ~XeTolerance();

 /**
 * Constructor 
 * @param  epsilon Accuracy for tolerance check
 */
    XeTolerance(double epsilon);

/**
 * Check that 2 double are equal within tolerance
 * @param d1 first value
 * @param d2 second value
 * @param print do we print a message in case of disagreement?
 * @return true if within tolerance
 */
    bool check(double d1,double d2,bool print=true);  

/**
 * Check that 3 double are equal within tolerance
 * @param d1 first value
 * @param d2 second value
 * @param d3 third value
 * @param print do we print a message in case of disagreement?
 * @return true if within tolerance
 */
    bool check(double d1,double d2,double d3,bool print=true);  

/**
 * Check that 4 double are equal within tolerance
 * @param d1 first value
 * @param d2 second value
 * @param d3 third value
 * @param d4 fourth value
 * @param print do we print a message in case of disagreement?
 * @return true if within tolerance
 */
    bool check(double d1,double d2,double d3, double d4, bool print=true);  

/**
 * Check that n double are equal within tolerance
 * @param d  values
 * @param n number of values
 * @param print do we print a message in case of disagreement?
 * @return true if within tolerance
 */
    bool check(double* d, int n,bool print=true);
/**
 * Check that n double inside a vector are equal within tolerance
 * @param d vector of values
 * @param print do we print a message in case of disagreement?
 * @return true if within tolerance
 */
    bool check(vector<double> &d, bool print=true);

  protected :
   
    double epsilon; /*!< absolute (not relative)  tolerance */

};

/**
   * Class for linear/parabolic interpolation
*/

class XeInterpolation: public XeMath {

  public : 


/**
 * Void constructor for ROOT
 */
    XeInterpolation();

/**
 * Destructor
 */
   ~XeInterpolation();


 /**
     * general interpolation procedure  
     * @return Interpolated quantity
     * @param  mode  LINEAR or PARABOLIC
     * @param  x table of x's
     * @param  y table of y's
     * @param  n size of input points
     * @param  value value at which interpolation is performed
 */
    static double interpolate(int mode,double *x,double *y,int n,double value);

/**
     * Interpolation from an internal fraction 
     * @return Interpolated quantity
     * @param  mode  LINEAR or EXPONENTIAL
     * @param  fraction where between x0 (faction=0) and x1 (fraction=1)
     * @param  y0 value at x0
     * @param  y1 value at x1
*/
    static double interpolate(int mode, double fraction, double y0, double y1);
 
 /**
     * linear interpolation
     * @return Interpolated quantity
     * @param  x table of x's
     * @param  y table of y's
     * @param  n size of input points
     * @param  value value at which interpolation is performed
 */
    static double linear(double *x, double *y, int n, double value);


 /**
     * linear interpolation between x0 and x1
     * @return Interpolated quantity
     * @param  x table of x0 and x1
     * @param  y table of y0 and x1
     * @param  value value at which interpolation is performed
     */
    static double linear(double *x, double *y, double value);

/**
     * Linear Interpolation from an internal fraction 
     * @return Interpolated quantity
     * @param  fraction where between x0 (faction=0) and x1 (fraction=1)
     * @param  y0 value at x0
     * @param  y1 value at x1
*/
    static double linear(double fraction, double y0, double y1);

 /**
     * parabolic interpolation
     * @return Interpolated quantity
     * @param  x table of x's
     * @param  y table of y's
     * @param  n size of input points
     * @param  value value at which interpolation is performed
 */
    static double parabolic(double *x, double *y, int n, double value);

  /**
     * parabolic interpolation between x0,x1,x2 
     * @return Interpolated quantity
     * @param  x table of x0, x1 and x2
     * @param  y table of y0, y1, and y2	
     * @param  n size of input points
     * @param  value value at which interpolation is performed
 */
    static double parabolic(double *x, double *y, double value);

/**
     * Exponential Interpolation from an internal fraction 
     * @return Interpolated quantity
     * @param  fraction where between x0 (faction=0) and x1 (fraction=1)
     * @param  y0 value at x0
     * @param  y1 value at x1
*/
    static double exponential(double fraction, double y0, double y1);
};


/**
   * General class to handle bins and contents.
   * It is described by nBins+1 vector edges and isn't optimised for special
   * cases of bin widths
*/
class XeBins : public XeObject {

  public:

/**
 * Void constructor for ROOT
 */
                     XeBins();

/**
 * Simple constructor
 * @param name  name of the object
 */
                     XeBins(string name);

/**
 * Full constructor
 * @param name  name of the object
 * @param edges nbins+1 values
 */
                     XeBins(string name, vector<double> &edges);

/**
 * Shorter method to get the bin corresponding to a value, when it can be done
 * @param x value to be tested
 * @return bin number, starting at 0
 */
    virtual int      getBin(double x);

    int              getNBins();
    double           getCenter(int i);
    double           getLowerEdge();
    double           getUpperEdge();
    double           getLowerEdge(int i);
    double           getUpperEdge(int i);
    double*          getLowerEdges();
    double*          getUpperEdges();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

/**
 * Get bin and fraction of bin for a given value
 * @param x value to be tested
 * @return a pair containing the bin and the fractional bin;
 */
    pair<int,double> getBinAndFraction(double x);

/**
 * Get a list of bins and weights for an interval.
 * The distribution is assumed flat. 
 * the sum of the weights is one only is the XeBins object fully covers it.
 * @param x1 lower edge of the interval
 * @param x2 upper edge of the interval
 * @return a map of bins to weight
 */
    map<int,double> getBinsAndWeights(double y1, double y2);


/**
 * Destructor
 */
    virtual         ~XeBins();

/**
 * print it, as requested because it's an XeObject 
 */
    bool             printIt(int level=1);
    bool             contains(double value);
                    
/**
 * Create a new Xebins, copy of this one, whose edges are log10(original)
 */
   XeBins* newLog10();



    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    void             extend(double xmin, double xmax);
    bool             isOk();

  protected:

    int              nBins;
    double           vMin;
    double           vMax;
    bool             ok;
    vector<double>   edges;
    map<double,int>  bins;   

    void             establishTheMap();

    ClassDef(XeBins,1)

};

/**
   * XeBins when the bins are equidistant
*/

class EquidistantBins : public XeBins {


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public: 

/**
 * Void constructor for ROOT
 */
   EquidistantBins();
   EquidistantBins(string what, int nBins,double a, double b);

/**
 * Destructor
 */
   virtual ~EquidistantBins();
   int getBin(double value);

  protected :

    double          step;
    static string   getEquidistantName(string w,int n, double a, double b);

    ClassDef(EquidistantBins,1)

};

/**
   * Bins when the contents are all equal 
*/

class EquiContentBins :public XeBins {


    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public:


/**
 * Void constructor for ROOT
 */
   EquiContentBins();
   EquiContentBins(string what, int nBins, vector<double>* data
                   ,double vmin=VERY_SMALL, double vmax=VERY_LARGE);

/**
 * Destructor
 */
   virtual  ~EquiContentBins();
  
  protected :

    static string   getEquiContentName(string w,int n);
    bool            fillFromMap(map<double,int>* data, int nBins, int n=0);

    ClassDef(EquiContentBins,1)

};


//------------- Linear and Logarithmic Ranges -----------------------------

/**
   * General class to handle intervals divided into bins.
   * Can be of type LINEAR, LOG, GENERAL.
   * Virtual class.
*/

class ErRange;
class YRange;
class S1Range;
class tValueRange;
class VEscRange;
class MassRange;
class SigmaRange;
class EventRange;
class LinearRange;
class LogRange;


/**
   * General class to handle range of points
   * It is described by nBins bins and isn't optimised for special
   * cases of bin widths
*/

class XeRange : public XeMath, public XeObject {


 public:


    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
/**
 * Constructor
 * @param what name of the quantity
 * @param n number of points
 * @param a lowest value
 * @param b hightest value
 * @param mode either LINEAR, LOG, or GENERAL
 */
   XeRange(string what, int n,double a, double b,int mode);
   
/**
 * Empty constructor for ROOT
 */
   XeRange();

/**
 * Destructor
 */
   virtual ~XeRange();
   
/** 
 *  get the mode
 * @return  LINEAR, LOG, or GENERAL
 */
   int getMode();

/**
 * Get the number of points
 */
   int getNPoints();  

/**
 * Get the number of bins
 */
   int getNBins();  

/**
    * return the index of the value closest to default 
*/ 
   int getReferenceIndex();

/**
    * return the lowest value
*/ 
   double  getMin();

/**
    * return the highest value
*/ 
   double  getMax();

/**
    * return the (logarithmic) step, whenever defined
*/ 
   double getStep();

 /**
     * return x-value for a given index
 */
   double getValue(int index);

 /**
     * return x-values 
     * @return         pointer to the x-values
 */
   double* getValues();

 /**
     * return vector of the x-values 
     * @return         pointer to the vector of x-values
 */
   vector<double>* getTheValues();

 /**
     * return center of "bin" between indexAndNextOne
 */
   double centerOfBin(int index);

 /**
     * return x-value for a index after given index
 */
   double getNextValue(int index);

 /**
     * return limits of a given interval
     * @return         pair of double
 */
   pair<double,double> getInterval(int i);

/**
 * print it, as requested because it's an XeObject 
 */
   bool printIt(int level=1);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

/**
 * Get the default S1 XeRange, instantiate it if necessary
 */
   static S1Range*     getDefaultS1Range();

/**
 * Get the default ERecoil XeRange, instantiate it if necessary
 */
   static ErRange*     getDefaultErRange();

/**
 * Get the default Y XeRange, instantiate it if necessary
 */
   static YRange*      getDefaultYRange();

/**
 * Get the default Mass XeRange, instantiate it if necessary
 */
   static MassRange*   getDefaultMassRange();

/**
 * Get the default t-value XeRange, instantiate it if necessary
 */
   static tValueRange* getDefaultTValueRange();

/**
 * Get the default V-escape XeRange, instantiate it if necessary
 */
   static VEscRange*   getDefaultVEscRange();

/**
 * Get the default Sigma XeRange, instantiate it if necessary
 */
   static SigmaRange*  getDefaultSigmaRange();

/**
 * Get the default Number of event XeRange, instantiate it if necessary
 */
   static EventRange*  getDefaultEventRange();

/**
 * check whether this range is equal to an already defined default (log or lin)
 * @return  The default range if found, NULL if not
*/
   XeRange* findDefaultRange();

/**
    * set the index of the value closest to reference
*/ 
   void  setReferenceIndex(double def);
              
/**
    * get the index of the value closest to a test valie
    * @param test the test value
*/ 
   int getClosestIndex(double test);
              
 /**
     * return index corresponding to a value
     * @param x input value
 */
   virtual int getIndex(double x)=0;

/** 
     * compute value corresponding to partial bin number
     * @param  fractIndex fractional index
*/
   virtual double computeValue(double fractIndex)=0;

/** 
     * compute value corresponding to index
     * @param  index index
*/
   virtual double computeValue(int index)=0;

/**
 * Is it linear?
 */
   bool           isLinear();

/**
 * Is it logarithmic?
 */
   bool           isLogarithmic();

/**
 * Is it general?
 */
   bool           isGeneral();

/**
 * Comparison operator
 */
   bool operator==(const XeRange& other);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
/**
 * Fill internal members, especially the vector of values
 */
   void           completeRange();

/**
 * Compute the name of the range
 */
   static string  getTheName(string name,int n,double v0, double v1);

 protected:

   int            mode;           /*!< either LINEAR, LOG, or GENERAL*/
   int            nPoints;        /*!< number of points */
   int            nBins;          /*!< number of "bins" = nPoints-1*/
   int            referenceIndex; /*!< index of value closest to default value*/
   double         lowest;         /*!< smallest value */
   double         highest;        /*!< largest value*/
   double         step;           /*!< step, whenever defined */
   vector<double> xs;             /*!< all values*/

   static S1Range      *defaultS1Range;     /*!< default S1(p.e.p) XeRange*/
   static ErRange      *defaultErRange;     /*!< default ERecoil XeRange*/
   static YRange       *defaultYRange;      /*!< default Y XeRange*/ 
   static MassRange    *defaultMassRange;   /*!< default Mass XeRange*/
   static tValueRange  *defaultTValueRange; /*!< default t-value XeRange*/
   static VEscRange    *defaultVEscRange;   /*!< default V-escape XeRange*/
   static SigmaRange   *defaultSigmaRange;  /*!<default Sigma XeRange */
   static EventRange   *defaultEventRange;  /*!< default nEventsXeRange*/

   ClassDef(XeRange,1)

 };

/**
   * The most general interval description
*/
class GeneralRange: public XeRange {

 public:
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
/**
 * Constructor
 * @param what is begind described
 * @param n number of values
 * @param values list of values
 */
   GeneralRange(string what, int n,double* values);
   
/**
 * Empty constructor for ROOT
 */
   GeneralRange();

/**
 * Destructor
 */
   virtual ~GeneralRange();

 /**
     * return index corresponding to a value
     * @param x input value
 */
   int     getIndex(double x);
   double  computeValue(int i);
   double  computeValue(double i);
   
/**
 * set a value
 * @param index index
 * @param value value
 */
   void    setValue(int index, double value);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
 protected: 

   map<double,int>  bins;   

   ClassDef(GeneralRange,1)

};

/**
   * Log-equidistant intervals
*/

class LogRange : public XeRange {

 public:

/**
 * Regular constructor
 * @param what quantity in the range
 * @param n nimber of values
 * @param a lowest value
 * @param b hightest value
 */
   LogRange(string what, int n,double a, double b);

/**
 * Void constructor for ROOT
 */
   LogRange();

/**
 * Destructor
 */
  ~LogRange();

 /**
     * return index corresponding to a value
     * @param x input value
 */
   int    getIndex(double x);
   double computeValue(int i);
   double computeValue(double i);

   ClassDef(LogRange,1)

};

/**
   * Linear-equidistant intervals
*/

class LinearRange : public XeRange {

 public:

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    

/**
 * Regular constructor
 * @param what quantity in the range
 * @param n nimber of values
 * @param a lowest value
 * @param b hightest value
 */
   LinearRange(string what, int n,double a, double b);

/**
 * Void constructor for ROOT
 */
   LinearRange();

/**
 * Destructor
 */
   virtual ~LinearRange();

 /**
     * return index corresponding to a value
     * @param x input value
 */
   int    getIndex(double x);
   double computeValue(int i);
   double computeValue(double i);
 
/**
 * Converts a value to a bin number and fraction
 * returns bin -1 if outside the range
 * @return pair representing the bin and the fraction
 * @param  x value to be converted
*/
    pair<int,double>    getIndexAndFraction(double x);

/**
 * Converts a value to a bin number and fraction
 * Returns lowest of highest possible value if outside the range
 * @return pair representing the bin and the fraction
 * @param  x value to be converted
*/
    pair<int,double>    getTolerantIndexAndFraction(double x);

   ClassDef(LinearRange,1)

};


/**
   * PRL Mass Range
*/

class PRLMassRange : public  GeneralRange {

   public :
  
     static constexpr  int    N_PRL_MASSES=38;

/**
 * Destructor
 */
    ~PRLMassRange();

/**
 * Void constructor for ROOT
 */
     PRLMassRange();
   
   protected:

     static double prlMasses[PRLMassRange::N_PRL_MASSES];

   ClassDef(PRLMassRange,1)
};


/**
   * Log Mass Range
*/

class MassRange : public  LogRange {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   static constexpr   int    DEFAULT_NM   =   25;
   static constexpr   double DEFAULT_MMIN =    6;
   static constexpr   double DEFAULT_MASS =   50;
   static constexpr   double DEFAULT_MMAX = 1000;
   public :
     MassRange(int n=DEFAULT_NM,double m0=DEFAULT_MMIN,double m1=DEFAULT_MMAX);

/**
 * Destructor
 */
    ~MassRange();

   ClassDef(MassRange,1)

};

/**
   * PE range
*/

class S1Range : public LinearRange {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   static constexpr   int    DEFAULT_N_PE   = 201;
   static constexpr   double DEFAULT_PE_MIN =  0.;
   static constexpr   double DEFAULT_PE_MAX = 50.;
   public :
     S1Range(int n=DEFAULT_N_PE,double p0=DEFAULT_PE_MIN,double p1=DEFAULT_PE_MAX);

/**
 * Destructor
 */
    ~S1Range();

   ClassDef(S1Range,1)

};

/**
   * Escape Velocity range
*/
class VEscRange : public LinearRange {

   static constexpr   int    DEFAULT_NVESC   =  11;
   static constexpr   double DEFAULT_VESCMIN = 500.;
   static constexpr   double DEFAULT_VESC    = 544.;
   static constexpr   double DEFAULT_VESCMAX=  600.;
  
 public :

     VEscRange(int n=DEFAULT_NVESC,double e0=DEFAULT_VESCMIN
              ,double e1=DEFAULT_VESCMAX);

/**
 * Destructor
 */
     ~VEscRange();

   ClassDef(VEscRange,1)

} ;

/**
   * "t-value" range
*/
class tValueRange : public LinearRange {

   static constexpr   int    DEFAULT_NTVAL   =   9 ;
   static constexpr   double DEFAULT_TVALMIN = -2.0;
   static constexpr   double DEFAULT_TVAL    =   0.;
   static constexpr   double DEFAULT_TVALMAX=   2.0;
  
 public :

     tValueRange(int n=DEFAULT_NTVAL,double e0=DEFAULT_TVALMIN
              ,double e1=DEFAULT_TVALMAX);

/**
 * Destructor
 */
     ~tValueRange();

   ClassDef(tValueRange,1)

} ;



/**
   * Enery recoil range
*/
class ErRange : public LinearRange {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   static constexpr   int    DEFAULT_NER   =  201;
   static constexpr   double DEFAULT_ERMIN =   0.;
   static constexpr   double DEFAULT_ERMAX=  100.;
   public :
     ErRange(int n=DEFAULT_NER,double e0=DEFAULT_ERMIN,double e1=DEFAULT_ERMAX);

/**
 * Destructor
 */
    ~ErRange();

   ClassDef(ErRange,1)

};

/**
   * Sigma intervals (log)
*/

class SigmaRange : public LogRange {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   static constexpr   int    DEFAULT_NS    =     61;
   static constexpr   double DEFAULT_SMIN  = 1.E-46;
   static constexpr   double DEFAULT_SIGMA = 1.E-45;
   static constexpr   double DEFAULT_SMAX  = 1.E-40;
   public :
     SigmaRange(int n=DEFAULT_NS,double s0=DEFAULT_SMIN,double s1=DEFAULT_SMAX);

/**
 * Destructor
 */
    ~SigmaRange();

   ClassDef(SigmaRange,1)

};

/**
   * Sigma intervals (linear)
*/

class EventRange : public LinearRange {

    /* -------------------------------------------------------------
     *                        Basic methods 
     * -----------------------------------------------------------*/
                    
   static constexpr   int    DEFAULT_NS   =   51;
   static constexpr   double DEFAULT_SMIN =   0.;
   static constexpr   double DEFAULT_SMAX =  10.;
   public :
     EventRange(int n=DEFAULT_NS,double s0=DEFAULT_SMIN
                     ,double s1=DEFAULT_SMAX);

/**
 * Destructor
 */
    ~EventRange();

   ClassDef(EventRange,1)

};



/**
   * Range for the "Y" variable in Form Factor computations
*/

class YRange : public LinearRange {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   static constexpr   int    DEFAULT_NY   =  51;
   static constexpr   double DEFAULT_YMIN =  0.;
   static constexpr   double DEFAULT_YMAX=  10.;
   public :
     YRange(int n=DEFAULT_NY,double y0=DEFAULT_YMIN,double y1=DEFAULT_YMAX);

/**
 * Destructor
 */
    ~YRange();

   ClassDef(YRange,1)

};


/**
   * Wrapper for integration
*/


class integrable : virtual public XeMath {
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public:

/**
 * Void constructor for ROOT
 */
                   integrable();

/**
 * Destructor
 */
    virtual       ~integrable();

 /**
     * @return  sum {f(x) dx}
     * @param   what Flag telling what to integrate
     * @param   x0 lower bound of integration
     * @param   x1 upper bound of integration
 */
    double         integrate(int what ,double x0, double x1) ;

 /**
     * @return value of the function to be integrated
     * @param  what  Flag telling what to compute
     * @param  x     Where the function is evaluated
 */
    virtual double getValue(int what,double x)=0;
};


/**
   * Wrapper for one dim equation solver.
   * It is automatically scaled in order to find a solution between 0 and 1. 
*/

class solvable: virtual public XeMath {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public:

/**
 * Void constructor
 */
                   solvable();

/**
 * Constructor
 * @param xmin lower bound of the solution
 * @param xmax upper bound of the solution
 */
                   solvable(double xmin,double xmax);

/**
 * Destructor
 */
    virtual       ~solvable();
/** 
 * 'f' to be implemented in the equation f(x)=0
 */
    virtual double getValue(double x)=0;

/**
 * Interface with the TF1 function
 * @param x pointer to parameter between 0 and 1
 * @apram p unused
 */
    double Evaluate(double *x, double *p);

 /**
     * return x such f(x)=y
     * @param  y  value for which solution is searched  
 */
    double solve(double y);

 /**
     * @param  xmi minimum value of solution
     * @param  xma maximum value of solution
 */
  void   setSolverLimits(double xmi,double xma);

  /**
 *  Set min and max
     * @param  xmi minimum value of solution
     * @param  xma maximum value of solution
 */
   void setLimits(double xmi,double xma);

  protected:

   double xmin;                               /*!< lower bound of the solution*/
   double xmax;                               /*!< Upper bound of the solution*/
   double range;                                                 /*!<xmax-xmin*/
   TF1*   theFunction;           /*!<Pointer to built-in function to be solved*/
         
};

#endif
