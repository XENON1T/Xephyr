#ifndef XeCore_h
#define XeCore_h

#include "XeVersion.h"

#include <iomanip>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <math.h>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TGLayout.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGWindow.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TPolyLine.h"
#include "TPolyMarker.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector.h"

#define deleteWithPointer(ptr)  {if(ptr!=NULL){ delete ptr;ptr=NULL;}}

#define RESET_PRECISION resetiosflags(ios_base::floatfield )<<setprecision(6)

enum bandType  { NO_BAND = -1 
               , LOWER_EDGE  
               , BAND_CENTER  
               , UPPER_EDGE   
               , N_EDGES      
               } ; 

enum symmetryType { SYMMETRIC
                  , ASYMMETRIC
                  } ;

                   
enum statusType { UNCHECKED, IN_ERROR, CHECKED_OK};

enum displayType { UNKNOWN_DISPLAY    = -1
                 , S2_VS_S1           
                 , S2_OVER_S1_VS_S1  
                 , FLATTENED_S2_VS_S1
                 , BAND_VS_SLICE
                 , N_DISPLAYS       
                 , CURRENT_DISPLAY    =  N_DISPLAYS
                 } ;

static const  double VERY_LARGE =                  1.E19
                  , TINY            =             1.E-99
                  , VERY_SMALL      =        -VERY_LARGE
                  , VERY_SMALL_LOG  =            -10000.
                  , UNDEFINED       =        -9999
                  , AUTOMATIC       =       UNDEFINED*10.
                  , CRAZY           =  98765432123456789.
                  ;


enum XeFonts      { DEFAULT_FONT=62
                  , LATEX_FONT=132
                  } ;

enum myColors     { darkGreen=kGreen+1 
                  };

enum trackingPoint { ENTRY
                   , EXIT
                   } ;

enum lineStyle    { PLAIN            = 1
                  , DASHED           = 2
                  , DOTTED           = 3
                  , DASH_DOTTED      = 4
                  , LONG_DASH_DOTTED = 5
                  , LONG_DASHED      = 7
                  } ; 

enum styleDefault { DEFAULT_LINE_COLOR   = kBlue
                  , DEFAULT_LINE_WIDTH   = 2
                  , DEFAULT_LINE_STYLE   = PLAIN
                  , DEFAULT_MARKER_COLOR = kBlack
                  , DEFAULT_MARKER_STYLE = 7       
                  , DEFAULT_FILL_COLOR   = 0
                  , DEFAULT_FILL_STYLE   = 0      
                  , CUT_LINE_COLOR       = kSpring
                  , CUT_FILL_WIDTH       = 102
                  , CUT_FILL_COLOR       = kSpring
                  , CUT_FILL_STYLE       = 3001
                  } ;

static const  double DEFAULT_MARKER_SIZE=1.; // no effect on markers style 1,6,7   

enum legendMode     {WITH_LEGEND, NO_LEGEND};

enum flagType     { VERY_LARGE_INT =  9999999
                  , VERY_SMALL_INT = -VERY_LARGE_INT
                  , UNDEFINED_INT  
                  , SAME     
                  , NEXT
                  , AUTO        
                  , ALL
                  , NONE 
                  , GENERAL
                  , LINEAR
                  , LOG
                  , EXPONENTIAL
                  , LOGX_LINY
                  , PARABOLIC
                  } ;

const string UNDEFINED_STRING="???";
const string XENON_100_REFERENCE
 ="https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis";

const string LABEL_BAND            = "Band number"
                  , LABEL_ER              = "E_{r} (KeV)"
                  , LABEL_EVT_AFTER_SEL   = "Events after selection"
                  , LABEL_EVT_DAY         = "Events (day^{-1})"
                  , LABEL_EVT             = "Events"
                  , LABEL_EVT_MAX         = "Events_{max}"
                  , LABEL_EVT_KEV         = "Events (KeV^{-1})"
                  , LABEL_EVT_KG_DAY      = "Events (Kg^{-1} day^{-1})"
                  , LABEL_EVT_KG_DAY_KEV  = "Events (Kg^{-1} day^{-1} KeV^{-1})"
                  , LABEL_EVT_PE          = "Events (p.e^{-1})"
                  , LABEL_FF              = "F^{2}"
                  , LABEL_FLATTENED       = "Flattened log_{10}(S2/S1)"
                  , LABEL_LEFF            = "L_{eff}"
                  , LABEL_LEFF_TVALUE     = "t-Value for L_{eff}"
                  , LABEL_MASS            = "Mass (GeV/c^{2})"
                  , LABEL_NORMTO1         = "Distribution norm. to 1"
                  , LABEL_LOG_LIKELIHOOD  = "-log(Likelihood)"
                  , LABEL_PVALUE          = "P-value"
                  , LABEL_Q_EXCLUSION     = "q_{excl.}"
                  , LABEL_S1              = "S1 (p.e.)"
                  , LABEL_S2              = "S2"
                  , LABEL_S2_OVER_S1      = "S2/S1"
                  , LABEL_SA              = "SA(y)"
                  , LABEL_SIGMA           = "\\sigma (cm^{2})"
                  , LABEL_SIGMA_MAX       = "\\sigma_{max} (cm^{2})"
                  , LABEL_SLICE           = "Slice Number"
                  , LABEL_VESC            = "V_{esc}"
                  , LABEL_Y_UOVER2        = "y=u/2"
                  ;

static const  int  DEFAULT_N_BANDS         = 12  
               ,  DEFAULT_N_SLICES        =  27
               ;

class stopper;

/**
   *  A collection of static methods for general handling 
*/
class XeCore {

   public  :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                       
     virtual ~XeCore();
     XeCore();

 /**
     * set debug level
     * @param l  level
 */
     static void    setDebugLevel(int l);

 /**
     * set trace level
     * @param l  level. (0=none, 1=partial, 2=full)
 */
     static void    setTraceLevel(int l);


 /**
     * Switch the "show references" flag on/off
 */
     static void    showTheReferences(bool show=true);

 /**
     * Suppress warnings
 */
     static void    suppressWarnings(bool suppress=true);

 /**
     * Re-enables warnings
 */
     static void    showWarnings(bool show=true);

 /**
     * Return debug level
 */
     static int     getDebugLevel(); 

 /**
     * Return trace level
 */
     static int     getTraceLevel(); 

 /** 
 * Fill an histogram with a vector of double
 */
     static void fillHist(TH1* h, vector<double> &v);

 /**
     * Cumulates a 1-dim histogram, bin by bin
     * @return Pointer to the cumulated histogram 
     * @param  h original histogram
 */
     static TH1F*   getCumulated(TH1F* h);

/**
   * Save all Xephyr objects to a new file (still not working)
   * @param fileName Name of the file
*/
    static void saveObjects(string fileName);

/**
   * Restore all Xephyr objects from a file (still not working)
   * @param fileName Name of the file
*/
    static void restoreObjects(string fileName);


    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     static bool    doWarn();
     static void    resetPrecision();
     static void    restoreWarnings();
     static void    saveWarnings();
     static void    displayFlag(string title,int tab,bool value,bool def);
     static void    displayFlag(string title,int tab,int value,int def);
     static void    displayFlag(string tit,int tab,double val,double def);
     static void    displayFlag(string title,int tab,string what,bool def);
     static bool    isALaTeXName(string s);
     static bool    isALaTeXChar(char c);
     static bool    fileExists(string path);
     static string  addASpace(string s);
     static string  makeItAFileName(string s);
     static string  makeItALaTeXName(string s);
     static string  upperCase(string s);
     static string  lowerCase(string s);
     static string  trim(string s);
     static string  justify(string s, int w, bool left,bool trim);
     static string  leftJustify(string s,int w, bool trim=false);
     static string  rightJustify(string s, int w, bool trim=false);
     static string  format0I(int v, int w=1);
     static string  formatI(int v, int w=1, int trailer=0);
     static string  formatLI(int v, int w=1, int p=1);
     static string  formatRI(int v, int w=1, int p=1);
     static string  formatF(double v, int w=1, int p=1);
     static TString  formatR(double v, int w=1, int p=1);
     static string  formatLF(double v, int w=1, int p=1);
     static string  formatRF(double v, int w=1, int p=1);
     static string  formatG(double v, int w=1, int p=1);
     static string  formatLG(double v, int w=1, int p=1);
     static string  formatRG(double v, int w=1, int p=1);
     static string  formatLatex(double v, int p=1);
     static string  yesOrNo(bool b);
     static string  doOrDont(bool b);
     static string  removeHeader(string original, string header);
     static string  replaceSlashes(string s);
     static string  replaceBlanks(string s);
     static string  replaceParentheses(string s);
     static string  positiveFlag(int flag);
     static vector<double> decodeVector(string text);
     

     static bool    initializeCore();
     static string getEdgeName(int e);

                  
   protected :

     static int      debugLevel;
     static int      traceLevel;
     static bool     coreInitialized;
     static bool     withWarnings;
     static bool     savedWarnings;
     static stopper *mainStopper;



} ;


/**
   *  A collection of static methods for gaphics
   *
*/
class XeGraphics : virtual public XeCore {

  public:
                    
     XeGraphics();
     virtual ~XeGraphics();

 /**
     * Draw a frame with limits and legends
     * @param title  frame title
     * @param xmi  lower x limit
     * @param xma  upper x limit
     * @param ymi  lower y limit
     * @param yma  upper y limit
     * @param tx    x title
     * @param ty    y title
 */
     static void   newFrame( string title,double xmi,double xma
                             , double ymi,double yma,string tx="",string ty=""
                             , int xMode=LINEAR, int yMode=LINEAR);
 /**
     * Draw a frame for 3 D flots with limits and legends
     * @param title  frame title
     * @param xmi  lower x limit
     * @param xma  upper x limit
     * @param ymi  lower y limit
     * @param yma  upper y limit
     * @param zmi  lower z limit
     * @param zma  upper z limit
     * @param zscale z scale, LOG or LINEAR?
     * @param tx    x title
     * @param ty    y title
 */
     static void   newFrameZ( string title,double xmi,double xma
                               , double ymi,double yma,double zmi,double zma
                               , int zscale=LINEAR,string tx="",string ty="");

/**
 * Create a canvas and open it
 * @param name  canvas name
 * @param nx divide it in nx columns of sub canvases
 * @param ny divide it in ny rows of subcanvas
 * @param square make it square
 */
   static void newCanvas(string name="", int nx=1 ,int ny=1, bool square=false);
/**
 * Save current canvas in pnf format, looking for ./pictures  or pictures dir.
 * If such directories don't exist, save in current directory
 * @param header  string to be added to the name
 */  
     static void   saveCanvas(string header="");
  
 /**
 *  Define a sub canvas, and increment  window;
 *  Called with no argument will keep divison an increment the sub canvas
 *  @param  nx Number of vertical sub-canvas; SAME means no change
 *  @param  ny Number of horizontal sub-canvas; SAME means no change
 *  @param  window where to plot, starting at zero; can be SAME or NEXT
 */
     static void  subCanvas(int nx,int ny, int window=0);

/** 
 * Switches to next subCanvas
 */
     static void nextSubCanvas();

/** 
 * Switches to given subCanvas
 *  @param  window where to plot, starting at zero; can be SAME or NEXT
 */
     static void subCanvas(int window);

 /**
 *  Define a floating sub canvas
 * @param xmin  Lower x-edge in overall canvas, from 0. to 1.
 * @param xmin  Upper y-edge in overall canvas, from 0. to 1.
 * @param ymin  Lower y-edge in overall canvas, from 0. to 1.
 * @param ymin  Upper y-edge in overall canvas, from 0. to 1.
 */
     static void   subCanvas(double xmin,double xmax,double ymin,double ymax);

/**
 * Update the graphics = force output
 */
     static void   updateGraphics();

/**
 * Draw a 1D histogram with reasonable titles and options
 * @param h pointer to histogram to be drawn
 * @param tx Title for x-axis 
 * @param ty Title for y-axis 
 * @param opt option to be passed to h->Draw();
*/
     static void   drawHist(TH1* h,string tx="",string ty="",string opt="");

/**
 * Opotionnaly Draw a 1D histogram with reasonable titles and options
 * @param h pointer to histogram to be drawn
 * @param ysc yScale LINEAR, LOG, or NONE (do nothing)
 * @param tx Title for x-axis 
 * @param ty Title for y-axis 
 * @param opt option to be passed to h->Draw();
*/
 static void  drawHist(TH1* h,int ysc, string tx="",string ty="",string opt="");


/**
 * Draw a 2D histogram in a Canvas, with reasonable titles and options
 * @param h pointer to scatter plot to be drawn
 * @param tx Title for x-axis 
 * @param ty Title for y-axis 
 * @param zscale LINEAR, LOG, or NONE (do nothing)
*/
     static void drawHist(TH2* h,string tx="",string ty="",int zscale=LINEAR);


/**
 * Draw the cumulated content of an histogram
 * @param color which color to be drawn
 */
     static TH1F*  drawCumulated(TH1F* h, int color=kRed);

/**
 * Set the font the cumulated content of an histogram
 * @param font  which font  to be used
 */
     static void   setFont(int font);

/**
 * Get the font 
 */
     static int    getFont();

   /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
 
/**
 * Check that a flag is either LOG or LINEAR
 * Optionally Print a warning message if not
 * @param flag value to be tested
 * @param warn print a warning in case of error
 */
     static bool    isScaleOK(int flag,bool warn=true);


/** 
 * Return 0 or 1 for ROOT, depending on given scale
 * @return 1 if scale==LOG, 0 otherwise
 * @param scale graphic mode being checked 
 */
     static int    isLog(int scale);

/**
 * Return error if flag isn't either LOG or LINEAR
 * Optionally Print a warning message if not
 * @param flag value to be tested
 * @param warn print a warning in case of error
 */
     static bool    isScaleError(int flag,bool warn=true);

  protected:

     static int        font;                                  /*!<Current font*/
     static TCanvas*   canvas;                              /*!<Current canvas*/
     static TPad*      pad;                                    /*!<Current pad*/
     static int        xCanvas;     /*!<Number of vertical subCanvas in canvas*/
     static int        yCanvas;    /*!<Number of horizontal subCancs in canvas*/
     static int        wCanvas;                          /*!<current subcanvas*/

}; 

/**
 * Repository of development tools
 */
 
class XeTool : virtual public XeCore {
    
    public:

      XeTool();
      virtual ~XeTool();
};

/**
   * Simpler wrapper to stop watch
*/

class stopper : public XeTool, public TStopwatch {

   public :

     stopper(string name);
    ~stopper();

 /**
     * restart the stopper;
 */
     void restart();

 /**
     * stop the stopper and print the CPU time
 */
     void print();

/**
 * print lapse time
 */
     void lapse();


   protected :

     string name; 
};

    

/**
   * A static class to count method invocation 
*/

class MethodCounter : public XeTool {
  
  public :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
               ~MethodCounter();
                MethodCounter();

 /**
     * count the number of times a method was invoked
     * @param method name of the method
 */
    static void count(string method);

 /**
     * list the number of times all methods were invoked
 */
    static void list();

 /**
     * reset all counters
 */
    static void reset();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  protected :

     static map<string,int>  counts;

} ;

class XeTopObject;

/**
   * Basic Xephyr Object.
   * Later on, will extend TObject
*/
class XeObject : virtual public XeCore, virtual public TObject {

  
  public:


    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/

/**
 * Virtual method which checks the object, fill its status 
 * @return if checked OK
 * @param recomputeIfNeeded Should we compute if needed?
 */
    virtual  bool checkIt(bool recomputeIfNeeded=true);

/**
 * Virtual method to trace its internal flags
  * This is a workaround for invoking "traceTheFlags" in interactive ROOT
 */
    virtual  void traceFlags();

/**
 * Desctructor
 */
    virtual ~XeObject();

 /**
 * Empty constucutor for ROOT
 */ 
   XeObject();

 /**
 * Regular constructor 
 * @param name object name
 */ 
   XeObject(string name);

/**
 * Trace the method entry point or exit point
 * @param tp  tracking point (either ENTRY or Exit)
 * @param method method name
 */
  void trace(trackingPoint tp, string method);

    void     setLaTeXName(string latex);
    void     setReference(string ref);
    void     printReference(int ntab=0);

/**
 * Set the legend (for graphs, etc...
 * @param legend the legend
 */
    void     setLegend(string legend="");

/** 
 * Set the t-value legend
 * @param t t-value
 */
  void setTValueLegend(double t);

/** 
 * Set the mass legend
 * @param mass mass
 */
  void setMassLegend(double mass);

/** 
 * Set the band legend
 * @param band band number
 */
  void setBandLegend(int band);

/**
 * set the object name
 * This method can be overriden
 */
    virtual  void     setName(string n);
    virtual  string   getFrameName();

  /**
   * print the XeObject properties
   * This is a workaround for invoking "bool printIt(int)" in interactive ROOT
   * @return true if no error was found
   * @param level  print level
 */
    bool  print(int level=1);


/**
* Check the object,  recompute what is needed, and fill the status. 
* This is a workaround for invoking "bool checkIt(bool)" in interactive ROOT
* @return flag telling that everything is ok
* @param recomputeIfNeeded Should we compute if needed?
 */
    bool     check(bool recomputeIfNeeded=true);

 /**
     * @return  This object has a legend (for MultiGraphs)
 */
    bool     hasLegend();

 /**
    static XeTopObject* theInstance;
     * @return  This object has a reference
 */
    bool     isReferenced();

 /**
     * @return  This object has a name
 */
    bool     hasName();


 /**
     * get the Legend name 
 */
    string   getLegend();

 /**
     * get the a point to the name 
 */
    const char*   getNameChar();

 /**
     * get the name 
 */
    string   getName();

 /**
     * get the name and add a space
 */
    string   getNameSpace();


/**
  * get the father in the tree
  * @return a pointer to the father
*/
    XeObject* getUp();

/** 
 * get one object down in the tree
 * @param g sequence number of the object
 * @return pointer to the object
*/
    XeObject* getDown(int g);

/** 
 * get  the vector of down objects
 * @return pointer to the object
*/
    vector<XeObject*>* getDowns();


/**
  * get the number of objects below the current one
*/
    int       getNDowns();
 
 /**
     * get the name in LateX format
 */
    string   getLaTeXName();

 /**
     * get the name of the reference
 */
    string   getReference();

 /**
 * get create a graphical canvas with the object name
 * @param nx divide it in nx columns of sub canvases
 * @param ny divide it in ny rows of subcanvas
 * @param square make it square
 */
     void   openCanvas(int nx=1, int ny=1,bool square=false);
 
 /**
     * Simple wrapper to Write to an already opened file
     * @param name of the object to be written. If none, use the XeGraph name;
 */
    void     write(string name="");

 /**
     * Simple wrapper to read from an already opened file
     * @param name of the object to be read. If none, use the XeGraph name;
 */
    void     read(string name="");

 /**
     * Write an object to a newly recreared file containing this object only
     * @param fName name of the file to be written. If none, use the object name;
     * @param oName name of the object to be written. If none, use the object name;
 */
    void     writeToSimpleFile(string fName="",string oName="");

 /**
     * Read an object from a newly recreared file containing this object only
     * @param fName name of the file to be written. If none, use the object name;
     * @param nName name of the object to be read;
     * @param name of the object to be read. If none, use the file name;
     * 
 */
    void     readFromSimpleFile(string fName="", string oName="");

/**
 * Print the tree belowe current object
 * @param depth how deep we go
 * @param offset tabulation offset
 */
    void     printTree(int depth=12345,int offset=0);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    static void         setResultsDirectory(string dir="");
    static string       getResultsDirectory();
    static bool         withReference();
    static void         showTheReferences(bool show=true);

 /**
     * return a pointer to the top object
 */
    static XeTopObject *getTop();
    void                setFrameName(string frame);
    void                setXenon100Reference(string xenonRef);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/

    void     addName(string n);
    void     setName(const char* n);

/** 
 * Return the object status
 * @return UNCHECKED, CHECKED_OK or IN_ERROR
 */
    int  getStatus();

/** 
 * return the object status name
 */
    string  getStatusName();

/** 
 * return the object status name
 */
    string  getMiniStatusName();

/** 
 * Return the object status name, given its name
 * @param s status
 */
   static string  getStatusName(int s);

/** 
 * Return the object mini status name, given its name
 * @param s status
 */
   static string  getMiniStatusName(int s);

/** 
 * Is the object checked OK?
 */
    bool isCheckedOK();

/** 
 * Is the object either OK or unchecked (i.e. not in error) ?
 */
    bool isOK();

/** 
 * Is the object  in error) ?
 */
    bool isError();

/** 
 * mark the object in error
 */
    void markError();

/** 
 * mark the object OK 
 */
    void markOK();

/** 
 * mark the object unchecked
 */
    void markUnchecked();

/**
 *  mark the result of a check as CKECKED_OK, IN_ERROR
 *  @param ok the result went ok
 */
    void markFromCheck(bool ok);

/**
 * return cumulated status 
 */
    bool checkAll();

 /**
     * track changes in objects 
     * @return did the object changed?
 */
    bool     isChanged();

 /**
     * mark the object as changed
     * @param propagateUp: the upper objects are also marked as changed
 */
    void     markAsChanged(bool propagateUp=true);

 /**
     * Do nothing (introduced for legibility reasons)
 */
    void     doNotMarkAsChanged();
 
 /**
     * mark the object as recomputed
 */
    void     markAsRecomputed();

/**
 * attach a daughter object
 * @param daughter pointer to daughter object
 */
    void     attach(XeObject* daughter);

/**
 * detach a daughter object and reattach it to the top object 
 * @param daughter pointer to daughter object
 */
    void     detach(XeObject* daughter);

 /**
 * detach an object from its father and reattach it to the top object 
 * @param up pointer to new father XeObject
 */
    void     detachItself();

 /**
 * Update pointer to father. 
 * Warning! internal method
 * @param up pointer to father XeObject
 */
    void     linkUp(XeObject* up);

/**
 * Virtual method to print the object
 * @return true if no error was found
 * @param level  print level
 */
    virtual  bool printIt(int level=1);
                     
/** 
 *  * set experiment number
 *   * @param exp either ALL or given run number
 *    */
    void     setExperiment(int exp=NONE);

/** 
 *  * get experiment number
 *   * @return either ALL or given run number
 *    */
    int      getExperiment();

/**
 *  * make sure that experiment was set
 *   */
    bool     isExperimentSet();


  protected:



    static string resultsDirectory;
    static bool   showReference;
    string        Name;
    string        LaTeXName;
    string        FrameName;
    string        Reference;
    string        legend;
    bool          changed;
    int           status;
    int           experiment;  // for which experiment/number number it applies 
    string            upName;  // needed to restore the tree
    XeObject*         up;  //!    do not save
    vector<XeObject*> downs;   

/**
 * Internal method for initialization
 */ 
    void              initializeObject(string n=UNDEFINED_STRING);

/**
 * Virtual method to trace its internal flags
 */
    virtual  void traceTheFlags();



};

/**
   * XeTopObject : the top object to which everything is attached
*/

static XeTopObject* top=NULL;

class XeTopObject : public XeObject {

  public :


 /**
     * return a pointer to the top object
 */
    static void         newInstance();

    XeTopObject();
   ~XeTopObject();

  protected: 

};

/**
   * An extension to map<int,double>
*/
class XeTable : public map<int,double> , public XeObject {

  public :

             XeTable();
             XeTable(string name);
    virtual ~XeTable();

 /**
     * add an entry
     * @param i index
     * @param value value corresponding to index 
 */
    void     add(int i,double value);

 /**
     * @return first index value 
 */
    int      getFirst();

 /**
     * @return last index value 
 */
    int      getLast();
/**
 * method to print the object
 * @return true if no error was found
 * @param level  print level
 */
    bool     printIt(int level=1);

  protected :


    int      first;
    int      last;

};

typedef map<int,double>::iterator XeTableIterator;



/**
   *  S1/S2 flatterner used when this display mode is needed
*/

class Flattener : virtual public XeObject {

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  public :

     virtual ~Flattener();
             Flattener();
             Flattener(string name);

 /**
     * @return flattened log s2/s1
     * @param   S1 original S1
     * @param   S2 original S2
 */
    virtual double flatten(double S1, double S2)=0;

 /**
     * @return S2 from flattened log s2/s1
     * @param   S1 original S1
     * @param   flat flattened quantity
 */
    virtual double unflatten(double S1, double flat)=0;

 /**
     * Unflatten a graph of S1/S2
     * @param origin The original graph in S1, flattened log(S2/S1)
     * @param added   Output [cumulated] graph
 */
    void     unflatten(TGraph* origin,TGraph* added);

  
};

/**
   * Handling graphics styles
*/

class XeStyle : virtual public XeGraphics, virtual public XeObject  {

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   public :
 
 /**
     * Constructor
     * @param  params (see name of parameter)
 */

           XeStyle( int lineColor=DEFAULT_LINE_COLOR
                  , int lineWidth=DEFAULT_LINE_WIDTH
                  , int lineStyle=DEFAULT_LINE_STYLE  
                  , int markerColor=DEFAULT_MARKER_COLOR
                  , double markerSize=DEFAULT_MARKER_SIZE
                  , int markerStyle=DEFAULT_MARKER_STYLE
                  , int fillColor=DEFAULT_FILL_COLOR
                  , int fillStyle=DEFAULT_FILL_STYLE);
    void   setLineColor(int color=DEFAULT_LINE_COLOR);

 /**
     * Set Line color according to predefined scheme
     * @param 
 */
    void   setLineStyleColor(int index=0);
    void   setLineStyle(int style=DEFAULT_LINE_STYLE);
    void   setLineWidth(int width=DEFAULT_LINE_WIDTH);
    void   setLineProperties(int color=DEFAULT_LINE_COLOR
                            ,int width=DEFAULT_LINE_WIDTH
                            ,int style=DEFAULT_LINE_STYLE);
    int    getLineColor();
    int    getLineStyle();
    int    getLineWidth();

    void   setMarkerColor(int color=DEFAULT_MARKER_COLOR);
    void   setMarkerStyleColor(int index=0);
    void   setMarkerStyle(int style=DEFAULT_MARKER_STYLE);
    void   setMarkerSize(double size=DEFAULT_MARKER_SIZE);
    void   setMarkerProperties(int color=DEFAULT_MARKER_COLOR
                            ,double size=DEFAULT_MARKER_SIZE
                            ,int style=DEFAULT_MARKER_STYLE);
    int    getMarkerColor();
    int    getMarkerStyle();
    double getMarkerSize();

    void   setFillColor(int color=DEFAULT_FILL_COLOR);
    void   setFillStyle(int style=DEFAULT_FILL_STYLE);
    void   setFillProperties(int color=DEFAULT_FILL_COLOR
                            ,int style=DEFAULT_FILL_STYLE);
    int    getFillColor();
    int    getFillStyle();



 /**
     * get a color according to a "reasonnable" predefined scheme   
     * @return Color corresponding to index x
     * @param  index Color index
 */
    static int getStyleColor(int index);

 /**
     * get a color according the rainbow
     * @return Color corresponding to index, number of points, and symtri
     * @param  index index (0 is at purple, last one is red)
     * @param  nC number of indices  (if SYMMETRIC, nC=2*colors+1)
 */
    static int getRainbowColor(int index, int nC, int mode=ASYMMETRIC);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    

 /**
     * Get the 'normalized style' while plotting t-values graphs   
     * @return  pointer to the XeStyle
     * @param  t requested t value
 */
    static XeStyle* forTValue(double t);    

    static XeStyle* forEdge(int edge);    

          ~XeStyle();

/**
 * method to print the object
 * @return true if no error was found
 * @param level  print level
 */
    bool     printIt(int level=1);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
 protected :

    int    lineColor;
    int    lineWidth;
    int    lineStyle;
    int    markerColor;
    double markerSize;
    int    markerStyle;
    int    fillColor;
    int    fillStyle;
    
};

/**
   * An object which can be drawn with style (actually XeStyle)
*/

class XeStylized  :  virtual public XeGraphics , virtual public XeObject {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * Virtual method applying the XeStyle to the object
 */
    virtual void  applyStyle();

 /**
     * Pure virtual method applying the XeStyle to the object
     * @param  options drawing options
     * @param  flag multipurpose flag
 */
    virtual void   draw(string options="L",int flag=0)=0; 

 /**
 * Pure virtual method returning X maximum value
 */
    virtual double getMaxX()=0;

 /**
 * Pure virtual method returning X minimum value
 */
    virtual double getMinX()=0;

 /**
 * Pure virtual method returning Y maximum value
 */
    virtual double getMaxY()=0;

 /**
 * Pure virtual method returning Y minimum value
 */
    virtual double getMinY()=0;

 /**
 * Pure virtual method returning Y minimum value for those different from 0
 */
    virtual double getMinYNotZero()=0;

 /**
     * Draw the object with its frame
     * @param  options drawing options
     * @param  flag multipurpose flag
 */
    void  drawWithFrame(string options="C",int flag=0); 

 /**
     * Draw the object with its canvas and frame
     * @param  options drawing options
     * @param  flag multipurpose flag
 */
    void  drawWithCanvasAndFrame(string options="C",int flag=0); 


 /**
     * Draw a Frame to hold the object
     * @param  xmi    MinX  (or AUTOMATIC)
     * @param  xma    MaxX  (or AUTOMATIC)
     * @param  ymi    MinY  (or AUTOMATIC)
     * @param  yma    MaxY  (or AUTOMATIC)
     * @param  xScale can be LINEAR, LOG, AUTO
     * @param  yscale can be LINEAR, LOG, AUTO
 */
   void     drawFrame( double xmi=AUTOMATIC, double xma=AUTOMATIC
                     , double ymi=AUTOMATIC, double yma=AUTOMATIC
                     , int xScale=AUTO     , int yscale=AUTO
                     ) ;
 /**
     * automatic generation of X label
 */
    virtual string   generatedXLabel();

 /**
     * automatic generation of Y label
 */
    virtual string   generatedYLabel();

    XeStylized();
 /**
     * @return  constructor
     * @param   Object name
     * @param   xL    X-label when object drawing
     * @param   yL    Y-label when object drawing
     * @param   zL    Z-label when object drawing
 */

    XeStylized(string name,string xL="", string yL="",string zL="");

 /**
     * Set the labels
     * @param   xL    X-label when object drawing
     * @param   yL    Y-label when object drawing
     * @param   zL    Z-label when object drawing
 */

    void     setXYLabels(string xL="", string yL="", string zL="");


 /**
     * set default x scale 
     * @param scale LINEAR or LOG
 */
    void setDefaultXScale(int scale);

 /**
     * set default y scale 
     * @param scale LINEAR or LOG
 */
    void setDefaultYScale(int scale);


 /**
     * get default X scale 
     * @return LINEAR or LOG
 */
    int getDefaultXScale();

 /**
     * get default Y scale 
     * @return LINEAR or LOG
 */
    int getDefaultYScale();
 
/**
     * Set the object style
     * @param  style XeStyle to be applied to the object
 */
    void     setStyle(XeStyle* style=NULL);

    void     setLineProperties(int color=DEFAULT_LINE_COLOR
                              ,int width=DEFAULT_LINE_WIDTH
                              ,int style=DEFAULT_LINE_STYLE);
    void     setLineColor(int color=DEFAULT_LINE_COLOR);

 /**
     * set the color according the rainbow
     * @param  index index (0 is at purple, last one is red)
     * @param  nC number of indices  (if SYMMETRIC, nC=2*colors+1)
     * @param  mode ASYMMETRIC or SYMMETRIC, whether we have one or two rainbows
 */
    void     setRainbowColor(int index, int nC, int mode=ASYMMETRIC);


 /**
     * Set line color according to predefined "reasonable" style
     * @param  index Index in the color scheme
 */
    void     setLineStyleColor(int index=0);
    void     setLineStyle(int style=DEFAULT_LINE_STYLE);
    void     setLineWidth(int width=DEFAULT_LINE_WIDTH);
    void     setMarkerProperties(int color=DEFAULT_MARKER_COLOR
                                ,double size=DEFAULT_MARKER_SIZE
                                ,int style=DEFAULT_MARKER_STYLE);
    void     setMarkerColor(int color=DEFAULT_MARKER_COLOR);

 /**
     * Set marker color according to predefined "reasonable" style
     * @param  index Index in the color scheme
 */
    void     setMarkerStyleColor(int index=0);
    void     setMarkerStyle(int style=DEFAULT_MARKER_STYLE);
    void     setMarkerSize(double size=DEFAULT_MARKER_SIZE);
    void     setFillProperties(int color=DEFAULT_FILL_COLOR
                              ,int style=DEFAULT_FILL_STYLE);
    void     setFillColor(int color=DEFAULT_FILL_COLOR);
    void     setFillStyle(int width=DEFAULT_FILL_STYLE);

    int      getLineColor();
    int      getLineStyle();
    int      getLineWidth();
    int      getMarkerColor();
    int      getMarkerStyle();
    int      getFillColor();
    int      getFillStyle();
    double   getMarkerSize();
    string   getXLabel();
    string   getYLabel();
    string   getZLabel();
    XeStyle* getStyle();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
/**
 * Optional draw of the stylized object
 * @param plot either NONE, AUTO, LINEAR, or LOG
 * @param  options drawing options
 * @param  flag multipurpose flag
 */
    void  drawIt(int plot, string options="L",int flag=0); 


    virtual ~XeStylized();
    void     printStyle();
    void     setXLabel(string xL);
    void     setYLabel(string yL);
    void     setZLabel(string zL);

  protected :
    
    XeStyle  style;
    string   xLabel;
    string   yLabel;
    string   zLabel;
    int      defaultXScale;
    int      defaultYScale;

}; 


/**
   * Extension to TGraph 
*/

class XeGraph : public XeStylized {
   
  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     XeGraph(); 

 /**
     * Constructor   with minimum information
     * @param name name of the graph
     * @param xLab   X-Label
     * @param yLab   Y-Label
 */
     XeGraph( string name, string xLab="", string yLab="");

 /**
     * constructor reading directly from a text file
     * @param fn     name of the file to be read from
     * @aram  name   name of the XeObject
     * @param Xlab   X-Label
     * @param Ylab   Y-Label
     * @param st     XeStyle of the graph (useful in MultiGraphs)
     * @param legend Legend to be drawn (useful in MultiGraphs)
 */
             XeGraph( string fn, string name ,int dummy, string Xlab=""
                    , string Ylab="" , XeStyle* st=NULL , string legend="");

 /**
     * constructor  with all details  
     * @param name   name of the object
     * @param n      number of points in the graph
     * @param x      pointer to array of X's
     * @param y      pointer to array of Y's
     * @param Xlab   X-Label
     * @param Ylab   Y-Label
     * @param st     XeStyle of the graph (useful in MultiGraphs)
     * @param legend Legend to be drawn (useful in MultiGraphs)
     * @param norm   Normalisation of overall integral
 */
             XeGraph( string name, int n,double *x, double *y, string Xlab=""
                    ,string Ylab="" , XeStyle* st=NULL ,string legend=""
                    , double norm=1.);
 /**
     * constructor  with all details  but the actual values
     * @param name   name of the object
     * @param n      number of points in the graph
     * @param Xlab   X-Label
     * @param Ylab   Y-Label
     * @param st     XeStyle of the graph (useful in MultiGraphs)
     * @param legend Legend to be drawn (useful in MultiGraphs)
 */
             XeGraph( string name, int n, string Xlab="" ,string Ylab="" 
                     , XeStyle* st=NULL ,string legend="");


 /**
     * constructor from a TGraph 
     * @param name   name of the object
     * @param  gr    original TGraph
     * @param Xlab   X-Label
     * @param Ylab   Y-Label
     * @param st     XeStyle of the graph (useful in MultiGraphs)
     * @param legend Legend to be drawn (useful in MultiGraphs)
 */
             XeGraph(string name , TGraph& gr  ,  string Xlab="", string Ylab=""
                    , XeStyle* st=NULL, string legend="");
 /**
     * void Draw a graph in S1S2 plane, according to current display mode
     * @param  opt extra Drawing options
 */
    void     drawS1S2(string opt="");

    void     set(int i);
    void     setPoint(int i, double x, double y);
    int      getN();
    double*  getX();
    double*  getY();
    double   getX(int i);
    double   getY(int i);
    double   getMinX();
    double   getMaxX();
    double   getMinY();
    double   getMaxY();

 /**
     * Get the minimum value of all values>0 
 */
    double   getMinYNotZero();

 /**
     * Get the sum of all X values
 */
    double   getSumX();

 /**
     * Get the sum of all Y values
 */
    double   getSumY();

    XeGraph* divide(XeGraph* reference);

/**
 * Get pointer to TGraph
 * @reurn the pointer
*/
  TGraph *getTGraph();

 /**
     * Draw the graph
     * @param  options drawing options
     * @param  flag = WITH_LEGEND or NO_LEGEND
 */
    void     draw(string options="L",int flag=WITH_LEGEND); 
    virtual ~XeGraph();

/**
 * Print the graph
 * @return true for compatibility reasons
 * @param n number of points (0 means all)
 * @param stat also print min,max,sum of X and Y
 * @param xNorm  all x values are mutilplied by xNorm
 * @param yNorm  all y values are mutilplied by xNorm
 */
    bool  printGraph(int n=0,bool stat=false,double xNorm=1.,double yNorm=1.);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    void     addGraph(TGraph &gr);
    void     drawInsideS1S2(string opt="");

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    void     setTheNames(string n,string leg="");
    void     applyStyle();
    bool     isCompatible(XeGraph* reference);
    void     printHeader();

/**
 * method to print the object
 * @return true if no error was found
 * @param level  print level
 */
    bool     printIt(int level=1);
  protected :
    

    TGraph   graph;
 
} ;

/**
   *  A graph (actually XeGraph) whose all Y values are positive
*/

class XeGraphPositive : public XeGraph {

   public :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     ~XeGraphPositive();

 /**
     * @param  n number of points
     * @param  *x  values of abcissae (omitted when correspoind y is <=0)
     * @param  *y  values of ordinates (omitted when <=0)
     * @param xLab x-labels
     * @param yLab y-labels
     * @param style XeStyle used for drawing the object
     * @param legend Legend used in XeMultiGraph
 */
      XeGraphPositive( string name, int n, double *x,double *y,string xLab=""
                     ,string yLab="" , XeStyle* style=NULL ,string legend="");
 /**
     *  Constructs the object from an existing TGraph
     * @param gr Original graph
     * @param xLab x-labels
     * @param yLab y-labels
     * @param style XeStyle used for drawing the object
     * @param legend Legend used in XeMultiGraph
 */
      XeGraphPositive(string name, TGraph* gr,  string xLab="", string yLab=""
                     , XeStyle* style=NULL, string legend="");
      XeGraphPositive(string name, XeGraph* gr);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    

   protected:

     void build(int n, double *x, double *y, XeStyle* style,string legend);

} ;


/**
   * Xephyr style TMultiGraph.
   * It is made of XeGraphs with all the same x vectors (known as Vx vector)
   * Each graph differs from the other one by a z value (konwn as Vz vector)
*/

class XeMultiGraph : public XeStylized , public vector<XeGraph*> {


   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     ~XeMultiGraph();

      XeMultiGraph();
 /**
     * Constructor of an XeMultiGraph from scratch
     * @param name  name of the object
     * @param xA  x-Label
     * @param yA  y-Label
     * @param zA  z-Label  (name of parameter varying from graph to graph)
 */
      XeMultiGraph(string name,string xA="",string yA="", string zA="");

 /**
     * Constructor of an XeMultiGraph from a list of existing graphs 
     * @param name  name of the object
     * @param ng Number of pre-existing Xefraphs
     * @param g  pointer to pre-existing XeGraphs
     * @param z  z-values; if NULL assume 0,1,2,3 ... 
     * @param xA  x-Label
     * @param yA  y-Label
     * @param zA  z-Label  (name of parameter varying from graph to graph)
 */
      XeMultiGraph(string name, int ng, XeGraph** g, double* z=NULL
                  ,string xA="",string yA="", string zA="");

 /**
     * Constructor of an XeMultiGraph from a list of existing graphs 
     * @param name  name of the object
     * @param v  vector of pre-existing XeGraphs
     * @param vz  z-values; if NULL assume 0,1,2,3 ... 
     * @param xA  x-Label
     * @param yA  y-Label
     * @param zA  z-Label  (name of parameter varying from graph to graph)
*/
      XeMultiGraph(string name, vector<XeGraph*>& v, vector<double>* vz=NULL
                  ,string xA="",string yA="", string zA="");

 /**
     * Normalize to an  XeGraph, ie divide each of this graphs by this one
     * @return pointer to a newly created XeMultigraph
     * @param reference Pointer XeGraph in the original one
 */
      XeMultiGraph* normalize(XeGraph* reference);

 /**
     * Normarlize to one of its graphs, XeMultiGraph,  ie divide each of this graphs by this one
     * @return pointer to a newly created XeMultigraph
     * @param  i sequence number (if AUTO), use the reference histogram
 */
      XeMultiGraph* normalize(int i=AUTO);

 /**
     * Draw the Multigraph
     * @param  options drawing options
     * @param  withLegend  with the legend box
 */
      void     draw(string options="L",int withLegend=0); 

 /**
     * Add an existing XeGraph
     * @param g  XeGraph to be added
     * @param z  The "z" value in the graph
     * (i.e. the quantity changing from graph to graph).
     * If missing, it will be the graph number;
 */
      void     add(XeGraph* g,double z=UNDEFINED);

 /**
     * Add an existing XeGraph, setting its color and legend
     * @param g  XeGraph to be added
     * @param color the color
     * @param legend the legend
 */
      void     add(XeGraph* g,int color, string legend="");

/**
 * Get minimum overall X
 */      
      double   getMinX();

/**
 * Get maximum overall X
 */      
      double   getMaxX();

/**
 * Get minimum overall Y
 */      
      double   getMinY();

/**
 * Get maximum overall Y
 */      
      double   getMaxY();

/**
 * Get minimum overall Y, ignoring zeros
 */      
      double   getMinYNotZero();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

  /**
     * Set the colors according the rainbow
     * @param  mode ASYMMETRIC or SYMMETRIC, whether we have one or two rainbows
 */
    void     setRainbowColors(int mode=ASYMMETRIC);

/**
   * Set the colors according the symmetric rainbow
 */
    void     setSymmetricRainbowColors();
                   
/*
   * Set the colors according the asymmetric rainbow
 */
    void     setAsymmetricRainbowColors();


   void     modifyLineColor(int color=DEFAULT_LINE_COLOR);
   void     modifyLineStyle(int style=DEFAULT_LINE_STYLE);
   void     modifyLineWidth(int width=DEFAULT_LINE_WIDTH);

/**
     * Transpose the multigraph by changing its  rows and columns
     * @return  pointer to a newly created multigraph
*/
     XeMultiGraph* transpose();

 /**
     * Set the reference graph (for further normalization)
     * @param  g the graph
 */
      void     setReferenceGraph(XeGraph *g);

 /**
     * Set the reference graph from its own list
     * @param  i sequence number of reference graph
 */
      void     setReferenceGraph(int i);

 /**
     * Get a graph 
     * @param i sequence number
     * @return  a pointer to the requested XeGraph
 */
      XeGraph* getGraph(int i);

 /**
     * Get the reference graph 
     * @return  a pointer to the reference XeGraph
 */
      XeGraph* getReferenceGraph();

 /**
     * Forces the labels according to those of one of the graphs (usually the 1st one)
     * @param which index of the reference graph
 */
      void     forceAxes(int which=0);

 /**
     * Add a multiGraph to the current one
     * @param detailed in detail
 */
      void           add(XeMultiGraph* mg);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
      string         generatedXLabel();
      string         generatedYLabel();
      void           applyStyle();


 /**
     * Access to invidual z-value 
     * @return  z-value 
     * @param  g  graph index
 */
     double          getZ(int g);

 /**
     * Access to array of  z-values 
     * @return pointer to  z-values 
 */
     double*         getZ();

 /**
     * Access to vector of  z-values 
     * @return pointer to  the vector 
 */
     vector<double>* getVz();

 /**
     * Access to invidual x-value 
     * @return  x-value 
     * @param  g  graph index
 */
     double          getX(int g);

 /**
     * Access to array of  x-values 
     * @return pointer to  x-values 
 */
     double*         getX();

 /**
     * Access to vector of  x-values 
     * @return pointer to  the vector 
 */
     vector<double>* getVx();
 
/**
     * Dump the graphs
     * @param level   print level
 */
     bool     printIt(int level=1);

   protected :


/** 
 * Reset all legends
 */
     void            resetTheLegends();

/** 
 * Internal method : push one legent
 * @param leg Legend to be added
 */
     bool            addLegend(string leg);

     set<string>     theLegends;     /*!<  Legends of indiviudal graphs*/
     XeGraph*        referenceGraph;     /*!< Pointer to reference graph */
     vector<double>  Vx;     /*!< Vector of reference X value*/
     vector<double>  Vz;     /*!< Vector of "z" values*/
    
} ;

/**
   *  A purely static class to manage S1/S2 display modes 
*/


class S1S2Display : virtual public XeGraphics {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * Set the display mode
     * @param   mode S2_VS_S1, S2_OVER_S1_VS_S1, FLATTENED_S2_VS_S1 
     * or BAND_VS_SLICE
     *
 */
     static void       setMode(int mode); 

/**
 * Set the limits 
 * @param xmi  Lowest X
 * @param xma  Highest X
 * @param ymi  Lowest Y
 * @param yma  Highest Y
*/
     static void       setLimits(double xmi,double xma,double ymi,double yma);

/**
 * Set the default limits
 */ 
     static void       setTheDefaultLimits();

 /**
 * Open a canvas with default display mode name
 * @param nx divide it in nx columns of sub canvases
 * @param ny divide it in ny rows of subcanvas
 * @param square make it square
 */
     static void       openCanvas(int nx=1, int ny=1,bool square=false);

/**
 * Draw a frame (holder for objects to be drawn) with current limits
 * @param title frame title
 */
     static void       drawFrame(string title);

/**
 * Draw a frame (holder for 3D objects to be drawn) with current limits in x & y
 * @param title frame title
 * @param zmin lowest z scale
 * @param zmax highest z scale
 * @param zScale LINEAR or LOG
 */
     static void       drawFrameZ( string title,double zmin,double zmax 
                                 , int zScale=LINEAR);
   /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    

     static bool       areCutsHatched();
     static bool       checkFull();
     static bool       isFullyDefined();
     static bool       isModeDefined();
     static double     getXmax();
     static double     getXmin();
     static double     getY(double x,double y);
     static double     getYmax();
     static double     getYmin();
     static double     undoY(double x,double y);
     static int        getDefaultLogMode(int m=CURRENT_DISPLAY);
     static int        getMode(int m=CURRENT_DISPLAY);
     static string     getCanvasName();
     static string     getName();
     static string     getXLabel();
     static string     getYLabel();
     static void       hatchTheCuts(bool h=true);
     static void       setNBands(int nb);
     static void       setNSlices(int ns);
     static void       setXmax(double xma);
     static void       setXmin(double xmi);
     static void       setYmax(double yma);
     static void       setYmin(double ymi);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
      
/**
 * Draw a box in current frame (internal method)
 * @param x0 lowest x corner of box
 * @param x1 highest x corner of box
 * @param y0 lowest y corner of box
 * @param y1 highest y corner of box
 */
     static void       drawBox(double x0,double x1,double y0,double y1
                              ,int color, bool ratio=false);

               
     S1S2Display();
    ~S1S2Display();

     static bool       requiresComputation(int m=CURRENT_DISPLAY);
     static void       setCurrentFlattener(Flattener* flatten);
     static Flattener* getCurrentFlattener();

  protected :

     static int        displayMode;
     static int        nBands;
     static int        nSlices;
     static bool       cutsAreHatched;  
     static double     xmin;
     static double     xmax;
     static double     ymin;
     static double     ymax;
     static Flattener* currentFlattener;

};



/**
   * An Object which can be displayed in various S1,S2 representations
*/
class S1S2Object : virtual public XeGraphics, virtual public XeObject  {
    
  public : 
  
/**
 * Is such an object drawable in a given S1S2 Display mode ?
 * By default, return yes but for BAND_VS_SLICE.
 * @param mode S2_VS_S1, S2_OVER_S1_VS_S1, FLATTENED_S2_VS_S1 
 * or BAND_VS_SLICE
 */
  virtual bool isDrawable(int mode);

/**
 * Empty constructor for ROOT
 */
             S1S2Object();
  
/**
 * Regular constructor
 * @param name Object name
 */
             S1S2Object(string name);
    virtual ~S1S2Object();

/**
 * Method to be implemented: draw the object
 */
    virtual  void   draw()=0;

/**
 * Method to be implemented: return min X in current S1S2Display mode
 */
    virtual  double minX()=0;

/**
 * Method to be implemented: return min Y in current S1S2Display mode
 */
    virtual  double minY()=0;

/**
 * Method to be implemented: return max X in current S1S2Display mode
 */
    virtual  double maxX()=0;

/**
 * Method to be implemented: return max Y in current S1S2Display mode
 */
    virtual  double maxY()=0;

/**
 * Set drawing limits automatic (true) or manual (false)
 */
    virtual  void   setAutomaticLimits(bool automatic=true);

/**
 * Virtual method to override drawFrame
 */
  virtual void   drawTheFrame();

/**
 * Draw a frame to hold the object
 */
  void   drawFrame();
 
/**
 * Draw a frame to hold the object with z color code
 * @param zmin  minimum Z
 * @param zmax  maximum Z
 * @param zscale LINEAR or LOG
 */
  void   drawFrameWithZ(double zmin, double zmax, int zscale);

/**
 * Set the S1S2 mode
 * @param   mode S2_VS_S1, S2_OVER_S1_VS_S1, FLATTENED_S2_VS_S1 
 * or BAND_VS_SLICE
 */
  void setS1S2Mode(int mode);

/**
 * get the current S1S2 mode
 * @return  S2_VS_S1, S2_OVER_S1_VS_S1, FLATTENED_S2_VS_S1 or BAND_VS_SLICE
 */
  int getS1S2Mode();


/**
 * Draw the object, possibly superimposing on previous one
 * @param mode S2_VS_S1, S2_OVER_S1_VS_S1, FLATTENED_S2_VS_S1, CURRENT_DISPLAY
 * or BAND_VS_SLICE
 * @param same draw in previous frame
 */
  void drawS1S2(int mode=CURRENT_DISPLAY,bool same=false);

  protected :

     bool  automaticLimits;              /*!< Drawings are in automatic limits*/
   
};

/**
  * Helper Class to send consistently formatted messages. It allows different levels 
  of severity for the messages: Info, Warning and Error. Error actually throws an 
  error with the std::runtime_error() which will abort (no so) gracefully and quit
  Xephyr, whatever it is doing.
*/

class errorHandler {
   public:

    //	Constructor: @param name: Tag name, typically here you wanna put the
    //  name of the class for which you wanna handle errors.
	errorHandler(TString name);

    // Send a message with ERROR severity, Xephyr will quit.
    // @param functionName: Tag identifier, function or in general reference place 
    // in the code where this message has been produced. 
    // @param message: message to the user.
    void Error(TString functionName, TString message);
    // Send a message with ERROR severity, Xephyr will quit.
    // @param functionName: Tag identifier, function or in general reference place 
    // in the code where this message has been produced. 
    // @param message: message to the user.
    void Warning(TString functionName, TString message);
    // Send a message with ERROR severity, Xephyr will quit.
    // @param functionName: Tag identifier, function or in general reference place 
    // in the code where this message has been produced. 
    // @param message: message to the user.
	void Info(TString functionName, TString message);

	TString className;	
};

#endif

