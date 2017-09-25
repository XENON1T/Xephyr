#ifndef XePhys_h
#define XePhys_h
#include "XeStat.h"
#include "XeMath.h"
/* 
   Physics analysis for Xenon
   (c) Daniel.Lellouch@weizmann.ac.il
   Originally based on code by Kaixuan Ni 

*/

   enum  targetType {
         XE100_MIXTURE =  -3
       , DEPLETED      =  -2  
       , NATURAL       =  -1
       , NEUTRON       =   0
       , PROTON       
       , N_NUCLEONS   
       , XENON         =  54
       } ;
    
   const int N_XENON_MODES=3;
   const int XENON_MODES[N_XENON_MODES]={129,131,XE100_MIXTURE};

//  Units
    const double 
        unity         = 1.
      , GeV           = unity                  // unit for masses
      , TeV           = 1.e3 *GeV
      , MeV           = 1.e-3*GeV
      , KeV           = 1.e-6*GeV           // unit for recoil energy
      , eV            = 1.e-9*GeV
      , cm            = unity
      , fm            = 1.e-13*cm                  
      , cm2           = cm*cm               // unit for cross section
      , barn          = 1.e-24*cm2              
      , mbarn         = 1.e-3*barn       
      , mubarn        = 1.e-6*barn       
      , nbarn         = 1.e-9*barn    
      , pbarn         = 1.e-12*barn   
      , fbarn         = 1.e-15*barn   
      , kg            = unity               // unit for detector mass
      , p_e           = unity               // photo-electons
      , day           = unity               // unit for exposure
      ; 

//  Constants
    const double  
        PI            = 3.141592653589
      , C             = 299792.                    // km/s
      , GF            = 1.16637e-5                 // GeV^-2
      , HBAR_C_FM     = .197327                    // Gev*Fm
      , HBAR_C_SQ_CM2 = HBAR_C_FM*HBAR_C_FM*fm*fm  // Gev^2cm^2
      , ATOMIC_MASS   = 0.93149                    // GeV/c^2
      , PROTON_MASS   = 0.93827                    // GeV/c^2
      , NEUTRON_MASS  = 0.93956                    // GeV/c^2
      ;

//  Galactic models       
    const double 

        DEFAULT_V_0   = 220.                   // Maxwell DM velocity
      , DEFAULT_V_ESC = 544.                   // Escape velocity at earth
      , DEFAULT_V_E   = 220.                   // Earth velocity
      , DEFAULT_RHO   = 0.3                    // In Gev/cm^3

      , V_0_SCALE     = 230.                   // In rate formulae
      , RHO_SCALE     = 0.4                    //  ditto
      , SIGMA_0_SCALE = 1.*pbarn               //  ditto
      , RATE_0        = 503.                   //  ditto
      ;

//  Just for constructor
    const double  
        DEFAULT_WIMP_MASS     = 50.            // Default wimp mass
      , DEFAULT_SIGMA_NUCLEON = 1E-45          // Default cross section
      ;

// For Er Spectrum, needs to be very precise!
     const int    N_ERS_POINTS          =   1000 ;
     const double ERS_STEP              =    0.1 
                , ERS_MAX               =    ERS_STEP * N_ERS_POINTS
                ;

// For LEff and acceptance tabulations, ERecoil for 0 to 100 KeV
     const int    N_ER_POINTS      =   200 ;
     const double ER_STEP          =   0.5 
                , ER_MAX           =   ER_STEP * N_ER_POINTS
                ;

// Spectrum in PE
     const int    N_PE_POINTS             =  200  ;
     const double PE_STEP                 =  0.25  
                , PE_MAX                  =   PE_STEP*N_PE_POINTS
                ;

     const int    UNSMEARED_PE_MAX        =   60  
                , UNSMEARED_PE_MIN        =    1  
                , N_PE_NR_BACKGROUND =   75   // For n bkg estimation
                ;
;
// Default values of flags
     const bool     DEFAULT_SMEARING_PMT = true;

// For integrations
     enum           integrationMode{RATE};
 

// For Spin Factors and structure functions

   enum SPIN_CONTENT_MODEL   { RESSELL_BONN_SPIN
                             , RESSELL_NIJMEGEN_SPIN
                             , SCHWENK_SPIN
                             , SUHONEN_SPIN
                             , N_SPIN_CONTENT
                             } ;

    enum XENON_NUCLEUS       { XENON_129  
                             , XENON_131 
                             , N_S_NUCLEI 
                             } ;

    enum SI_INTERACTION_MODE { ENGEL_SI
                             , FRICK_SI
                             , SUHONEN_SI
                             , EXPONENTIAL_SI
                             , N_SI_FF
                             , DEFAULT_SI_FF = FRICK_SI
                             , N_COMPARE_SI = 3
                             } ;


   enum SD_INTERACTION_MODE { FF_SCHWENK_1BC        
                            , FF_SCHWENK_2BC       
                            , RESSELL_BONN        
                            , RESSELL_NIJMEGEN   
                            , SUHONEN_SD        
                            , FIRST_SD_SF
                            , SF_SCHWENK_1BC        =  FIRST_SD_SF
                            , SF_SCHWENK_2BC  
                            , N_SD_INTERACTIONS     =  SF_SCHWENK_2BC + 1
                            , N_SD_POLFF            =  RESSELL_NIJMEGEN + 1
                            , N_SD_FF               =  SUHONEN_SD + 1
                            , N_SD_SF               =  N_SD_INTERACTIONS-N_SD_FF
                            , N_SF_BANDS_COEFS      =  5 
                            , N_SCHWENK_Y           = 54
                            , N_SCHWENK_DELTA_A1    = 13 
                            , DEFAULT_SD_FF         =  FF_SCHWENK_1BC
                            , DEFAULT_SD_SF         =  SF_SCHWENK_2BC
                            , DEFAULT_SD_INTERACTION=  DEFAULT_SD_SF
                            } ;

      STATIC_CONST  int DEGREES_RESSELL        =  8 + 1 // 8th order + 1/(1+y)
                     , DEGREES_SCHWENK        =  9
                     , N_POLSD_COEFS          = 10
                     ;
 
      STATIC_CONST  int S00_SFACT              = 0
                     , S01_SFACT              = 1
                     , S11_SFACT              = 2
                     , N_SFACT                = 3
                     , S01_UPP                = 0
                     , S01_LOW                = 1
                     , S11_UPP                = 2
                     , S11_LOW                = 3
                     ;
/**
 * Physics class (for documentation purpose only)
 */

class XePhysics : virtual public XeMath, public XeObject {
  public:
   XePhysics();
   XePhysics(string name);
   virtual ~XePhysics();
};

/**
   * Model of a Galaxy.
   * Parameterized by density, Maxwell-Boltzman velocity,
   * escape and Earth velocities
*/


class GalaxyModel : public XePhysics {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     GalaxyModel( double rho=DEFAULT_RHO     , double v_0=DEFAULT_V_0
                , double v_esc=DEFAULT_V_ESC , double v_e=DEFAULT_V_E ) ;

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     double getRho();
     double getV_E();
     double getV_0();
     double getV_ESC();
     double getVMax();
     double getRoverR0();
     double getK0overK1();
     void   setRho(double rho);
     void   setV_0(double v_0);
     void   setV_E(double v_e);
     void   setV_ESC(double v_esc);
    ~GalaxyModel();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     double rateScale(double vmin);

   protected :

     double V_0;
     double V_ESC;
     double V_E;
     double RHO;

     double K0overK1;
     double RoverR0;
     void   updateKinematics();


} ;



/**
   * Description of a WIMP.
   * Currently parameterized by its mass only.
   * In the future other types of particles could be added
*/

class Wimp : public XePhysics {
  
  public  :
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    Wimp(double m=DEFAULT_WIMP_MASS);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    double getMass();
    void   setMass(double m);
   ~Wimp();

  protected :
  
    double mass;

};


/**
   * Description of a nucleus.
   * Defined by N and A.
   * Also characterized by mass and spin
*/

class Nucleus : public XePhysics {
  
  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @param  n  Atomic number, currently NEUTRON, PROTON, or XENON
     * @param  a  Nuclear number (needed for XENON)
 */

    Nucleus(int n, int a=0);

    int    getA();
    int    getN();
    bool   printIt(int level=1);
    double getJ();
    double getMass();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    double getOscillatorSizeSQ();
    double QtoEr(double Q);
    double QtoY(double Q);
    double ErtoQ(double Er);
    double ErtoY(double Er);
    static string nucleusName(int n,int a);
    static string spinName(double J);
    static double oscillatorParameterSQ(int A);
   ~Nucleus();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    static int    hashCode(int n, int a);
    static bool   isValid(int n, int a,bool quiet=false);
    static int    nFromHashCode(int hash);
    static int    aFromHashCode(int hash);

  protected:
    
    STATIC_CONST  double B0=  41.467;
    STATIC_CONST  double B1=  45.000;
    STATIC_CONST  double B2=  25.000; 

    int    A                ;
    int    N                ;
    double J                ;
    double mass             ; // in GeV
    double oscillatorSizeSQ ; // in fm

};


/**
   * Internal table of nuclei properties
*/


class LedererTable : public XePhysics {

   
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  public :

    static LedererTable* getInstance();

    LedererTable();
   ~LedererTable();
    Nucleus* getNucleus(int n, int a);
    Nucleus* getNucleus(int hash);
    bool     printIt(int level=1);
  
  protected :

    static LedererTable* instance;
    map<int,Nucleus*> nuclei;

};


// -------------------  Target: collection of nuclei -------------------------

class Interaction;

/**
   * Description of a target (sensitive gas/liquid).
   * Characterized by the element (currently Xenon only) and its abundance
   * pattern (either NATURAL or DEPLETED)
*/

class Target :public XePhysics, public integrable {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @param  n   Atomic number either XENON or other defined Nucleus
     * @param  a   abundance pattern, in case of Xenon: 
     *             either NATURAL or DEPLETED
 */

                      Target(int n, int a=NATURAL);
 /**
     * Set WIMP mass
     * @param mass The mass of the Wimp
 */
     void             setWimpMass(double mass);

 /**
     * Set the interaction
     * @param inter Interction
 */
     void             setInteraction(Interaction* inter);

 /**
     * Compute typical Maximum Energy recoil 
     * @return  maximum recoil energy 
 */
     double           ErMax();
 /**
     * Compute  differential ERecoil rate 
     * @return  differential rate per Kg per day per Kev 
 */
     double           dRate(double Er);
 /**
     * Compute  Integrated ERecoil rate 
     * @return  integrated rate per Kg per day between Ermin and Ermax
 */
     double           rate(double Ermin=0.,double Ermax=0.);
 /**
     * Compute a graph of Energu Recoil 
     * @return  a graph (actually XeGraph) of the rate per kg per day for
     * a given energy range (ErRange).
 */

     XeGraph*         newGraphOfRate(ErRange* er=NULL);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
                      Target();
 /**
     * set the Galaxy Model
     * @param  gal GalaxyModel
 */
     void             setGalaxyModel(GalaxyModel* gal);
     void             setWimp(Wimp *w);
     bool             printIt(int level=1);
     void             printInteraction();

     void             add(Target* target,double fraction);

 /**
     * Add a nucleus to the target compound
     * @param n Atomic number
     * @param a Mass number
     * @param f relative fraction 
 */
     void             add(int n, int a=0, double fraction=1.);
     double           getFraction(int n, int a);
     double           getWimpMass();
     map<int,double>& getFractions();
     Interaction*     getInteraction();
     GalaxyModel*     getGalaxyModel();
     Wimp*            getWimp();

     virtual         ~Target();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
/**
     * Check compatibility
     * @return flag telling that everything is ok
     * @param recomputeIfNeeded Recompute when thing aren't initialized
 */
    bool   checkIt(bool recomputeIfNeeded=true);

     double           getValue(int mode, double x);

   protected : 

     map<int,double>  fractions;
     Interaction     *interaction;
     GalaxyModel     *galaxyModel;
     Wimp            *wimp;
     void             init();
 
}  ;


/**
   * Spin content of protons and neutrons for various nuclei.
   * Used for Spin dependent calculations
   * This is a singleton
*/

 
class SpinContent :public XePhysics  {

  public:

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    static  string       getTheName(int m);
    static  string       getNucleusName(int n);
    static  SpinContent* getInstance(int model);
    static  int          getNucleus(int A);
    double               getProtonSpin(int A);
    double               getNeutronSpin(int A);

 /**
     * @return  spin content
     * @param  nucleon   either PROTON or NEUTRON
     * @param  A         nuclear number (current 129 or 131 only)
 */

    double               getSpin(int nucleon, int A);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  protected :

    static int           NUCLEI[N_S_NUCLEI];
    static double        SPIN_C[N_NUCLEONS][N_SPIN_CONTENT][N_S_NUCLEI];
           double        S_PROTON[N_S_NUCLEI];
           double        S_NEUTRON[N_S_NUCLEI];
           int           model;

    virtual             ~SpinContent();  
                         SpinContent(int model);  
    static  SpinContent* instances[N_SPIN_CONTENT];  

 
};


/**
   * Structure function for interactions (more general than form factors). 
   * This is a virtual class
*/


class StructureFunction :  public XePhysics {

  public :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    virtual         ~StructureFunction();
                     StructureFunction();
 /**
     * @param  w can be either SF_SCHWENK_1BC or SF_SCHWENK_2BC
 */
                     StructureFunction(int w);
    virtual void     printCoefficients();
 /**
     * @param   edge can be either NO_BAND , LOWER_EDGE , BAND_CENTER , UPPER_EDGE
 */

    virtual void     setBandEdge(int edge);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    string           getNameWithEdge();

    virtual XeStyle* getStyle();                             

  protected:

    int              edge;
    int              what;

 
} ;

/**
   * Spin Independent Structure function.
   * A virtual class
*/

class SIStructureFunction : public StructureFunction {

   public: 

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
                    SIStructureFunction(int w);
     virtual       ~SIStructureFunction();


 
};

/**
   *  Spin Dependent Structure function.
   * A virtual class
*/

class SDStructureFunction : public StructureFunction {


   public: 

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     XeGraph*       newGraphOfSA(Nucleus *nuc, int n, XeRange* yr=NULL);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     static  double Coefs[N_SD_SF][N_S_NUCLEI][N_NUCLEONS][N_POLSD_COEFS];
     static  SDStructureFunction* getStructureFunction(int p=UNDEFINED_INT);

                    SDStructureFunction(int p);
     virtual       ~SDStructureFunction();
     virtual void   resetCoefficients();
     virtual double SA(Nucleus *nuc, double y, int nucleon, bool print=false)=0;
   
     double         getYMax();

   protected :

     double         yMax;

};

/**
   * Spin Dependent Structure function a la Schwenk et al.
   * Defined by parameter SF_SCHWENK_1BC, SF_SCHWENK_2BC
*/

class SchwenkSDStructureFunction: public SDStructureFunction{

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
 public:

   static double    Bands[N_S_NUCLEI][N_NUCLEONS][2][N_POLSD_COEFS];

  ~SchwenkSDStructureFunction();
 /**
     * @param  p either SF_SCHWENK_1BC resp SF_SCHWENK_2BC (for 1 resp 2 currents) 
 */

   SchwenkSDStructureFunction(int p);

   void             printCoefficients();
 /**
     * @param   edge can be either NO_BAND , LOWER_EDGE , BAND_CENTER , UPPER_EDGE
 */

   void             setBandEdge(int edge);
   void             resetCoefficients();
   double*          getCoefficients(int A, int nucleon);
   double           SA(Nucleus *nuc, double y, int nucleon, bool print);
   XeStyle*         getStyle();

 protected:
  
   double           coefs[N_S_NUCLEI][N_NUCLEONS][N_POLSD_COEFS];

};



/**
   *  Form factor: an alternative to structure functions
*/

class FormFactor :  public XePhysics {


   public: 

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @param   edge can be either NO_BAND , LOWER_EDGE , BAND_CENTER , UPPER_EDGE
 */

     virtual void     setBandEdge(int edge);
     virtual void     printCoefficients();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     FormFactor();
     FormFactor(int w);
     virtual XeStyle* getStyle();                             
     virtual         ~FormFactor();
     string           getNameWithEdge();

   protected:

     int              edge;
     int              which;

 
};

/**
   * Internal table to keep Form Factors a la Suhonen
*/

class SuhonenTable : virtual public XeMath {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   public: 

     STATIC_CONST  int N_POINTS=2001; 
     static double    U_TABLE[N_S_NUCLEI][N_POINTS][2+N_SFACT]; 
     static double    OMEGA_0[N_S_NUCLEI];
     static double    OMEGA_1[N_S_NUCLEI];

     virtual         ~SuhonenTable();
                      SuhonenTable();
     int              getIndex(double u);
     double           interpolate(double u, int n, int index, int coef) ;
     
};


/**
   * Spin Independent Form Factor.
   * A virtual class
*/

class SIFormFactor :public FormFactor {

  /* virtual class */

   public: 

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     static SIFormFactor *getFormFactor(int w=UNDEFINED_INT); 
     virtual             ~SIFormFactor();
                          SIFormFactor(int w);
     virtual double       computeFF(double q,Nucleus *nuc)=0;

 
};

/**
   * Exponential Form Factor.
*/

class ExponentialFormFactor : public SIFormFactor{

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    protected: 
    
     STATIC_CONST  double Default_R0 = .3        ;
     STATIC_CONST  double Default_R1 = 0.91      ;
     STATIC_CONST  double Default_Coherence= 1.5 ;

     double R0       ;
     double R1       ;
     double Coherence;

   public: 

    ~ExponentialFormFactor();
     ExponentialFormFactor(double r0=Default_R0,double r1=Default_R1
                          ,double coherence=Default_Coherence); 
     void   setTheName();
     void   setR0(double r);
     void   setR1(double r);
     void   setCoherence(double c);
     double computeFF(double q, Nucleus* nuc); 

 
};

/**
   * Woods Saxon Form Factor.
*/

class WoodsSaxonFormFactor : public SIFormFactor {

  /* virtual class */

   public :
 
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
      virtual ~WoodsSaxonFormFactor();
               WoodsSaxonFormFactor(int w);
      void     setThickness(int t);
      double   computeFFWell(double q, double r2);

   protected:
     
      STATIC_CONST  double defaultThickness=1.0 ;
      virtual void        setTheName()=0;
      double              thickness;

 
};

/**
   * Spin Independent Form Factor a la Engel.
*/

class EngelSIFormFactor : public WoodsSaxonFormFactor {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   protected :
     STATIC_CONST  double defaultR0=1.2;

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     EngelSIFormFactor(double r0=defaultR0,double s=defaultThickness);
                    
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     void      setR0(double r0);
    ~EngelSIFormFactor();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     XeStyle*  getStyle();
     void      setTheName();
     double    computeFF(double q, Nucleus* nuc) ;


   protected :
     double R0;


};

/**
   * Spin Independent Form Factor a la Frick
*/

class FrickSIFormFactor : public WoodsSaxonFormFactor {

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   protected :

     STATIC_CONST  double defaultThicknessF=.9 ;
     STATIC_CONST  double defaultThicknessA=.52;
     STATIC_CONST  double defaultR0=1.23;
     STATIC_CONST  double defaultR1=-0.6;

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     FrickSIFormFactor(double r0=defaultR0, double r1=defaultR1
                    ,double s=defaultThicknessF, double a=defaultThicknessA);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     void      setR0(double r0);
     void      setR1(double r1);
     void      setThicknessA(double a);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     double    computeFF(double q, Nucleus* nuc);
     void      setTheName();
     XeStyle*  getStyle();
     virtual ~FrickSIFormFactor();

   protected :

     double R0;
     double R1;  
     double thicknessA;  


};

/**
   * Spin Independent Form Factor a la Suhonen
*/

class SuhonenSIFormFactor : virtual public SIFormFactor
                          , virtual public SuhonenTable{

   public:
     
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     SuhonenSIFormFactor();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/

    ~SuhonenSIFormFactor();
                    
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     double computeFF(double y,Nucleus *nuc);
     XeStyle*  getStyle();


};

/**
   * Spin Dependent Form Factor.
   * A virtual class
*/

                                          
class SDFormFactor : public FormFactor {

   public: 

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     static SDFormFactor* getFormFactor(int p=DEFAULT_SD_FF);
                          SDFormFactor(int p);
     XeGraph*             newGraphOfSA(Nucleus *nuc, int n, XeRange* yr=NULL);
     XeGraph*             newGraphOfFF(Nucleus *nuc, int n, XeRange* yr=NULL);
     XeGraph*             newGraphOfStructureFactor( Nucleus *nuc, int fc
                                                   , XeRange* yr=NULL);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
     static void          printStructureFactor(double *s);
     virtual             ~SDFormFactor();
   
     int                  preferredSpinModel(); 
     double               getYMax();
     void                 printDetailedStructureFactor( Nucleus *nuc
                                                      , XeRange* yr=NULL);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     virtual void         resetCoefficients();
     virtual void         computeStructureFactor(Nucleus *nuc, double y 
                                                , double *s)=0;
     double               FF(Nucleus *nuc, double y, double aN, double aP
                            , bool print=false);
     double               FF(Nucleus *nuc, double y, int nucleon
                            , bool print=false);
     double               SA(double* s ,double aN, double aP);
     double               SA(Nucleus *nuc, double y ,double aN, double aP
                            , bool print=false);
     double               SA(Nucleus *nuc, double y, int nucleon);

/**
     * Check compatibility
     * @return flag telling that everything is ok
     * @param recomputeIfNeeded Recompute when thing aren't initialized
 */
    bool   checkIt(bool recomputeIfNeeded=true);

     bool                 checkIt(Nucleus* nuc);
     static double        getAn(int mode);
     static double        getAp(int mode);
     static string        getStructureFactorName(int c);
     static string        getStructureFactorLaTeXName(int c);

   protected :

     double               yMax;


};
 
/**
   * Spin Dependent Form Factor parameterized by a polynomial.
   * Defined by RESSELL_BONN, RESSELL_NIJMEGEN, FF_SCHWENK_1BC, FF_SCHWENK_2BC
*/

                                          

class PolynomialSDFormFactor : public SDFormFactor {

   public:

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    static double Coefs[N_SD_POLFF][N_S_NUCLEI][N_SFACT][N_POLSD_COEFS];

    virtual ~PolynomialSDFormFactor();
             PolynomialSDFormFactor(int p);

    void     resetCoefficients();
    void     computeStructureFactor(Nucleus *,double y,double *s);
    void     printCoefficients();
    double*  getCoefficients(int A, int which);

   protected :

     int    nDegrees;
     void   setParameters(string nam, int ncoefs, double yMax);
     double coefs[N_S_NUCLEI][N_SFACT][N_POLSD_COEFS];


};
                    

/**
   * Spin Dependent Form Factor a la Ressel Bonn
*/

class RessellBonnSDFormFactor: public PolynomialSDFormFactor{
  public:

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    RessellBonnSDFormFactor();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
   ~RessellBonnSDFormFactor();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    XeStyle* getStyle();


};

/**
   * Spin Dependent Form Factor a la Ressel Nijmegen
*/

class RessellNijmegenSDFormFactor: public PolynomialSDFormFactor{

  public:
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    RessellNijmegenSDFormFactor();
    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
   ~RessellNijmegenSDFormFactor();
    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    XeStyle* getStyle();

};

/**
   * Spin Dependent Form Factor a la Schwenk, one body current
*/

class Schwenk1bcFormFactor: public PolynomialSDFormFactor{

 public:
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   Schwenk1bcFormFactor();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  ~Schwenk1bcFormFactor();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   XeStyle* getStyle();


};

/**
   * Spin Dependent Form Factor a la Schwenk, two bodies current
*/

class Schwenk2bcFormFactor: public PolynomialSDFormFactor{

 public:
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   Schwenk2bcFormFactor();


    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  ~Schwenk2bcFormFactor();
 /**
     * @param   edge can be either NO_BAND , LOWER_EDGE , BAND_CENTER , UPPER_EDGE
 */

   void     setBandEdge(int edge);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   XeStyle* getStyle();
   static double Bands[N_S_NUCLEI][N_SF_BANDS_COEFS][N_POLSD_COEFS];


};


class SchwenkDeltaA1  : public  XePhysics {

 public :
  
    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
   SchwenkDeltaA1();
   static XeGraph* newGraphOfSA( Nucleus*nuc, int nucleon, int delta_a1
                        , int normalizeTo=UNDEFINED_INT);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  ~SchwenkDeltaA1();
   static void     getSA(double *sa, Nucleus* nuc, int nucleon, int delta_a1
                        , int normalizeTo=UNDEFINED_INT);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
   static double yValues[N_SCHWENK_Y];
   static double deltaA1Values[N_SCHWENK_DELTA_A1];
   static double SA_NP[N_S_NUCLEI][N_NUCLEONS][N_SCHWENK_Y][N_SCHWENK_DELTA_A1];

   static XeStyle* getStyle();

   

};

class SchwenkYRange : public GeneralRange {
   public :

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     SchwenkYRange(int nb=N_SCHWENK_Y);
    ~SchwenkYRange();


};

/**
   * Spin dependent form factor a la Suhonen.
*/

class SuhonenSDFormFactor : virtual public SDFormFactor
                         , virtual public SuhonenTable {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
     SuhonenSDFormFactor();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    ~SuhonenSDFormFactor();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
     XeStyle* getStyle();
     void computeStructureFactor(Nucleus *nuc,double y,double *s);

};


/**
   * virutal class describing an interaction.
   * Chacracterized by a  Wimp, a nucleus, either a form factor or a structure ffunction
*/

class Interaction :  public XePhysics  {

  /* virtual class */
   
  public :


    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    virtual           ~Interaction();
                       Interaction(double s);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    void               setSigmaNucleon(double s);
    virtual double     SigmaErNucleus(Wimp* wimp, Nucleus* nuc,double Er)=0;
    virtual            SpinContent* getSpinContent();
    bool               isSpinDependent();
    double             getSigmaNucleon();
    double             SigmaEr(Wimp* wimp, Target* target, double Er) ;
    double             QMax(Wimp* wimp,double V, Target* target);
    double             ErMax(Wimp* wimp, double V, Target* target);
    double             YMax(Wimp* wimp, double V, Target* target);
    FormFactor        *getFormFactor();
    void               setFormFactor(FormFactor* f);
    StructureFunction *getStructureFunction();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    bool               withStructureFunction();
    bool               withFormFactor();
    void               updateKinematics(Wimp* wimp, Nucleus* nuc);
    string             getNameWithFF();
    string             getLegendFromFF();    
    XeStyle           *getStyleFromFF();

   protected :
  
    double             sigmaNucleon;
    double             m_D;
    int                A;
    double             m_A;
    double             mu_A;
    double             mu_p; 
    double             mu_n; 
    double             kin_p; 
    double             kin_n; 
    bool               spinDependent;
    FormFactor        *form;
    StructureFunction *structure;
 

} ;

/**
   * Spin independent interaction.
   * Characterized by nucleon cross section and Form Factor
*/

class SIInteraction : public Interaction {

  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @param  s sigma-nucleon
     * @param  f form factor
 */

    SIInteraction(double s, SIFormFactor*f);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
   ~SIInteraction();
 /**
     * @param  s sigma-nucleon
     * @param  ff ENGEL_SI , FRICK_SI , SUHONEN_SI , EXPONENTIAL_SI 
     * or DEFAULT_SI_FF = FRICK_SI
 */


    SIInteraction(double s=DEFAULT_SIGMA_NUCLEON, int ff=DEFAULT_SI_FF);
 /**
     * @return Differential cross section
     * @param  Wimp Incoming particle
     * @param  nuc  Nucleus
     * @param  Er   Energy recoil 
 */

    double               SigmaErNucleus(Wimp* wimp, Nucleus* nuc,double Er);
    SIStructureFunction *getSIStructureFunction();
    SIFormFactor        *getSIFormFactor();

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
  protected :

    void setParameters(SIFormFactor *f);


};

/**
   * Spin dependent interaction.
   * This is a virtual class
   * Defined by Neutron and Proton coupling
*/

class SDInteraction : public Interaction {
  
  public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
 /**
     * @param s nucleus cross section
     * @param an neutron coupling
     * @param ap proton coupling
     * @param what either FF_SCHWENK_1BC    , FF_SCHWENK_2BC    , RESSELL_BONN    , RESSELL_NIJMEGEN  , SUHONEN_SD    , SF_SCHWENK_1BC    , SF_SCHWENK_2BC 
 */

    SDInteraction(double s,double an,double ap,int what=UNDEFINED_INT);
 /**
     * @param s nucleus cross section
     * @param an neutron coupling
     * @param ap proton coupling
     * @param f Spin dependent form factor
 */
    SDInteraction(double s,double an,double ap,SDFormFactor* f);
 /**
     * @param s nucleus cross section
     * @param an neutron coupling
     * @param ap proton coupling
     * @param sf Spin dependent structure function
 */
    SDInteraction(double s,double an,double ap,SDStructureFunction* sf);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    virtual             ~SDInteraction();

 /**
     * @param   edge can be either NO_BAND , LOWER_EDGE , BAND_CENTER , UPPER_EDGE
 */

    void                 setBandEdge(int edge);
    void                 setSpinContentModel(int m);
    void                 setFormFactor(SDFormFactor *f);
    void                 setStructureFunction(SDStructureFunction *sf);
    double               SigmaErNucleus(Wimp* wimp,Nucleus* nuc,double Er);
    SpinContent         *getSpinContent();
    SDFormFactor        *getSDFormFactor();
    SDStructureFunction *getSDStructureFunction();

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
  protected :

    SpinContent         *spinContent;
    double               aN;
    double               aP;
    int                  nucleon;
    void                 setParameters(double an, double ap,SDFormFactor* f
                                      ,SDStructureFunction* sf);


};
   

/**
   * Spin dependent cross section, either pure neutron or pure proton
*/

class SDPureInteraction: public SDInteraction {

  public  :

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
    static   string getModeName(int m);

    virtual ~SDPureInteraction();
             SDPureInteraction(double s, int nucleon, SDFormFactor*f);
             SDPureInteraction(double s, int nucleon, SDStructureFunction* sf);
             SDPureInteraction(double s, int nucleon, int what);

    /* -------------------------------------------------------------
     *                Internal methods (not for user)
     * ------------------------------------------------------------*/
                    
    string   getModeName();

  protected :

    void     setNucleon(int m);


}; 

/**
   * Spin dependent cross section,  pure proton
*/

class SDPureProtonInteraction: public SDPureInteraction {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    SDPureProtonInteraction(double s,SDFormFactor *f=NULL);
    SDPureProtonInteraction(double s,SDStructureFunction* sf);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
   ~SDPureProtonInteraction();


};

/**
   * Spin dependent cross section, pure neutron 
*/

class SDPureNeutronInteraction: public SDPureInteraction {

   public :

    /* -------------------------------------------------------------
     *                     Basic methods 
     * ------------------------------------------------------------*/
                    
    SDPureNeutronInteraction(double s, SDFormFactor *f=NULL);
    SDPureNeutronInteraction(double s, SDStructureFunction* sf);

    /* -------------------------------------------------------------
     *                     Advanced methods 
     * ------------------------------------------------------------*/
                    
   ~SDPureNeutronInteraction();


};
 

#endif
