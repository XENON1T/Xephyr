/* 
   Physics analysis for Xenon
   (c) Daniel.Lellouch@weizmann.ac.il
   Originally based on code by Kaixuan Ni 

*/
#include "XePhys.h"
 
XePhysics::XePhysics() : XeObject() {}
XePhysics::XePhysics(string name) : XeObject(name) {}
XePhysics::~XePhysics(){}

// ----------------------- Galactic model ---------------------------------

GalaxyModel::~GalaxyModel(){}
GalaxyModel::GalaxyModel(double rho, double v_0, double v_esc, double v_e) :
 XePhysics() {
   V_0   = v_0;
   V_ESC = v_esc;
   V_E   = v_e;
   RHO   = rho;
   updateKinematics();
}

void GalaxyModel::updateKinematics(){
   markAsChanged();
   if(V_0<0.) {
     cout<<"Invalid V_0"<<endl;
     K0overK1=0.;
   }
   else {
     double x=V_ESC/V_0;
     K0overK1=1./(erf(x)-2/sqrt(PI)*x*exp(-x*x));       
     if(V_E<0.) {
       cout<<"V_E"<<endl;
       RoverR0=0.;
     }
     else {
     double y=V_E/V_0;
     RoverR0=K0overK1*(.5*(sqrt(PI)*(y+.5/y)*erf(y)+exp(-y*y))
                       -(x*x+y*y/3.+1.)*exp(-x*x));
     }
   }
   char t[100];
   sprintf(t,"Galaxy model : Rho=%.1f Gev/cm^3 V_0=%.1f V_esc=%.1f km/s"
         ,RHO,V_0,V_ESC);
   setName(t);
 }

 double GalaxyModel::rateScale(double Vmin){
   if(Vmin<V_ESC-V_E){
     return K0overK1*(
            sqrt(PI)/4*V_0/V_E*(erf((Vmin+V_E)/V_0)-erf((Vmin-V_E)/V_0))
           -exp(-V_ESC*V_ESC/V_0/V_0)
            );
   }
   else if(Vmin<V_E+V_ESC){
     return K0overK1/2./V_E*(
            sqrt(PI)/2.*V_0*(erf(V_ESC/V_0)-erf((Vmin-V_E)/V_0))
           -(V_ESC+V_E-Vmin)*exp(-V_ESC*V_ESC/V_0/V_0)
            );
   }
   return 0;
 }

 double GalaxyModel::getRho()           {return RHO;}
 double GalaxyModel::getV_E()           {return V_E;}
 double GalaxyModel::getV_0()           {return V_0;}
 double GalaxyModel::getV_ESC()         {return V_ESC;}
 double GalaxyModel::getVMax()          {return V_ESC+V_E;}
 double GalaxyModel::getK0overK1()      {return K0overK1;}
 double GalaxyModel::getRoverR0()       {return RoverR0;}
 void   GalaxyModel::setRho(double r)   {
   if(RHO==r) return;
   RHO=r;
   updateKinematics();
 }
 void   GalaxyModel::setV_E(double v)   {
   if(V_E==v) return; 
   V_E=v; 
   updateKinematics();
 }
 void   GalaxyModel::setV_0(double v)   {
   if(V_0==v) return; 
   V_0=v; 
   updateKinematics();
 }
 void   GalaxyModel::setV_ESC(double v) {
   if(V_ESC==v) return; 
   V_ESC=v; 
   updateKinematics();
 }
  
// -------------------------- WIMP --------------------------------------

Wimp::Wimp(double m) : XePhysics()  {setMass(m);}
Wimp::~Wimp(){}
   
double Wimp::getMass()         {return mass;}
void   Wimp::setMass(double m) { 
  if(mass==m) return;
  mass=m;
  char t[80];
  sprintf(t,"Wimp of mass %.1f GeV/c2",m);
  setName(t);
  Reference="Lewin & Smith, Astroparticle Physics 6 (1996) 87";
  markAsChanged();
}

// ----------------------------------  Nucleus -------------------------
 
bool Nucleus::isValid(int n, int a,bool quiet){
 switch(n) {
    case PROTON  :
    case NEUTRON : if(a==1) return true;
                   break;              
    case XENON   :
         switch(a){
           case 124: 
           case 126:
           case 128:
           case 129: 
           case 130: 
           case 131: 
           case 132: 
           case 134: 
           case 136: return true; 
         }
         break;
  }
  if(!quiet) cout<<"Can't handle nucleus with N="<<n<<", A="<<a<<endl;
  return false;
}

double Nucleus::oscillatorParameterSQ(int A){
  double a3=pow((double)A,-1./3.);
  return B0/a3/(B1-B2*a3);
}

int    Nucleus::getA()                {return A;}
int    Nucleus::getN()                {return N;}
double Nucleus::getJ()                {return J;}
double Nucleus::getMass()             {return mass;}
double Nucleus::getOscillatorSizeSQ() {return oscillatorSizeSQ;}

double Nucleus::QtoEr(double Q)  {return Q*Q/2./mass/KeV;}
double Nucleus::QtoY(double Q)   {
  double q=Q/HBAR_C_FM;
  return q*q*getOscillatorSizeSQ()/4.;
}
double Nucleus::ErtoQ(double Er) {return sqrt(2.*mass*Er*KeV);}
double Nucleus::ErtoY(double Er) {return QtoY(ErtoQ(Er));}

string Nucleus::nucleusName(int n,int a){
  switch(n){
    case NEUTRON  : return "n" ; 
    case PROTON   : return "p"  ; 
    case XENON    :
         switch(a) { 
           case DEPLETED      : return "Depleted Xe"; 
           case NATURAL       : return "Natural Xe";  
           case XE100_MIXTURE : return "Xenon100 mixture"; 
         }
         return "Xe "+formatI(a)   ;
  }
  return UNDEFINED_STRING;
}
  
string Nucleus::spinName(double J){
  int j=(int)(J*10+.1);
  switch(j){
    case  0: return "  0";
    case  5: return "1/2";
    case 10: return "  1";
    case 15: return "3/2";
  }
  return UNDEFINED_STRING;
}

bool Nucleus::printIt(int){
  cout<<"  Spin "<<spinName(J)<<", mass "<<setw(5)<<setprecision(1)<<fixed
       <<setw(0)<<mass<<" GeV/c2"<<endl;
  return true;
}

int  Nucleus::hashCode(int n, int a)  {return (n<<16)+a;}
int  Nucleus::nFromHashCode(int hash) {return hash>>16;}
int  Nucleus::aFromHashCode(int hash) {return hash&0xFFFF;}

Nucleus::~Nucleus() {}
Nucleus::Nucleus(int n, int a) : XePhysics(nucleusName(n,a)) {
  if(!isValid(n,a)) {
    markError();
    char ti[100];
    sprintf(ti,"Invalid nucleus with N=%d, A=%d",n,a);
    setName(ti);
    return;
  }
  switch(n) {

    case NEUTRON : J=1./2. ; 
         mass=NEUTRON_MASS; 
         setFrameName("Neutron");
         setLaTeXName("n");
         break;

    case PROTON  : J=1./2. ; 
         mass=PROTON_MASS;  
         setFrameName("Proton");
         setLaTeXName("p");
         break;

    case XENON   :
         switch(a){
           case 124: mass=123.906*ATOMIC_MASS; J=0   ; break;
           case 126: mass=125.904*ATOMIC_MASS; J=0   ; break;
           case 128: mass=127.904*ATOMIC_MASS; J=0   ; break;
           case 129: mass=128.905*ATOMIC_MASS; J=0.5 ; break;
           case 130: mass=129.903*ATOMIC_MASS; J=0   ; break;
           case 131: mass=130.905*ATOMIC_MASS; J=1.5 ; break;
           case 132: mass=131.904*ATOMIC_MASS; J=0   ; break;
           case 134: mass=133.905*ATOMIC_MASS; J=0   ; break;
           case 136: mass=135.907*ATOMIC_MASS; J=0   ; break;
                     return;
         }
         break; 
  }
  N=n;
  A=a;
  oscillatorSizeSQ=oscillatorParameterSQ(A);
  if(A==1) Reference="PDG";
  else {
    Reference="www.webelements.com";
  }
  LedererTable::getInstance()->attach(this);
}

// ------------------ LedererTable of isotopes ---------------------

LedererTable* LedererTable::instance=NULL;

LedererTable* LedererTable::getInstance() {
  if(instance==NULL) instance=new LedererTable();
  return instance;
}

LedererTable::~LedererTable(){
  for(map<int,Nucleus*>::iterator it=nuclei.begin();it!=nuclei.end();it++){
    Nucleus *nuc=it->second;
    delete nuc;
  }
  nuclei.clear();
}

LedererTable::LedererTable(): XePhysics("Private Lederer table"){}


Nucleus* LedererTable::getNucleus(int n, int a){
  return getNucleus(Nucleus::hashCode(n,a));
}

Nucleus* LedererTable::getNucleus(int hash){
  int n=Nucleus::nFromHashCode(hash);
  int a=Nucleus::aFromHashCode(hash);
  Nucleus *nuc=nuclei[hash];
  if(nuc==NULL) {
    nuclei[hash]=new Nucleus(n,a);
    nuc=nuclei[hash];
  }
  if(!nuc->isOK()) {
    cout<<"Can't find nucleus n="<<n<<", a="<<a<<endl;
    return NULL;
  }
  return nuc;
}

bool LedererTable::printIt(int level) {
  if(level<1) return true;
  int n=nuclei.size();
  cout<<"---"<<n<<" nuclei defined in "<<getName()<<endl;
  if(level<2) return true;
  for(map<int,Nucleus*>::iterator it=nuclei.begin();it!=nuclei.end();it++){
    Nucleus *nuc=it->second;
    if(nuc->isOK()){
      cout<<nuc->getName()<<" ";
      nuc->printIt();
    }
  }
  cout<<endl;
  return true;
}

// -------------------  Target: collection of------------

void Target::add(int n, int a, double f) {
  if(f<=0.) return;
  if(LedererTable::getInstance()->getNucleus(n,a)!=NULL) {
    int hash=Nucleus::hashCode(n,a);
    fractions[hash]+=f;
  }
}
   
Target::Target() : XePhysics("Target"), integrable()  {init();}

void Target::init(){
  interaction=NULL;
  galaxyModel=NULL;
  wimp=NULL;
}

bool Target::checkIt(bool recomputeIfNeeded){
  bool ok=true;
  if(wimp==NULL) {
    if(recomputeIfNeeded){
      cout<<"... no WIMP defined, take default, mass="
          <<DEFAULT_WIMP_MASS<<" GeV/c2"<<endl;
      setWimp(new Wimp());
    }
    else {
      cout<<" ERROR: no wimp defined"<<endl;
      ok=false;
    }
  }
  if(interaction==NULL) {
    if(recomputeIfNeeded){
      cout<<"... no interaction defined, take default SI, sigma-N="
         <<DEFAULT_SIGMA_NUCLEON<<")"<<endl;
      setInteraction(new SIInteraction());
    }
    else {
      cout<<" ERROR: no interaction defined"<<endl;
      ok=false;
    }
  }
  if(galaxyModel==NULL) {
    if(recomputeIfNeeded){
      cout<<"... no galaxy model defined, take default"<<endl;
      setGalaxyModel(new GalaxyModel());
    }
    else {
      cout<<" ERROR: no galaxy model defined"<<endl;
      ok=false;
    }
  }
  if(ok) markOK();
  else markError();
  return ok;
}

void Target::setInteraction(Interaction* inter){
  if(interaction==inter) return;
  detach(interaction); 
  interaction=inter;
  attach(interaction); 
  markAsChanged();
}

void Target::setGalaxyModel(GalaxyModel* gal)  {
  if(galaxyModel==gal) return;
  detach(galaxyModel);
  galaxyModel=gal;
  attach(galaxyModel);
}

void Target::setWimp(Wimp *w) {
  if(wimp==w) return;
  detach(wimp);
  wimp=w;
  attach(wimp);
}

Interaction* Target::getInteraction()                  {return interaction;}
GalaxyModel* Target::getGalaxyModel()                  {return galaxyModel;}
Wimp*        Target::getWimp()                         {return wimp;}

double Target::getWimpMass() {
  if(wimp==NULL) {
    cout<<"No wimp is defined!"<<endl;
    return 0.;
  }
  return wimp->getMass();
}

void  Target::setWimpMass(double mass) {
  if(wimp==NULL) wimp=new Wimp();
  wimp->setMass(mass);
}

double Target::ErMax(){
  if(!checkIt()) {
    cout<<"Target "<<getName()<<" isn't well defined!"<<endl;
    return 0.;
  }
  double vMax=galaxyModel->getVMax();
  return interaction->ErMax(wimp,vMax,this);
}

Target::~Target(){}

void Target::printInteraction(){
cout << "debug 1" << endl;
  if(!checkIt()) return;
cout << "debug 2" << endl;
  if(wimp==NULL) cout<<"No wimp defined"<<endl;
  else cout<<"Wimp            : "<<wimp->getName()<<endl;
cout << "debug 3" << endl;

  if(interaction==NULL) cout<<"No interaction defined"<<endl;
  else {
cout << "debug 4" << endl;
    cout<<"Interaction     : "<<interaction->getName()<<endl;
//    cout<<"Form factor     : "<<interaction->getFormFactor()->getName()<<endl;
cout << "debug 5" << endl;
    if(interaction->isSpinDependent()){
cout << "debug 6" << endl;
//    cout<<"Spin content    : "<<interaction->getSpinContent()->getName()<<endl;
cout << "debug 7" << endl;
    }
  }

  if(galaxyModel==NULL) cout<<"No galactic model defined"<<endl;
  else cout<<"Galactic Model  : "<<galaxyModel->getName()<<endl;
cout << "debug 8" << endl;
}

bool Target::printIt(int){
  double sumA=0.;
  double sumMass=0.;
  for(map<int,double>::iterator it=fractions.begin();it!=fractions.end();it++){
    Nucleus *nuc=LedererTable::getInstance()->getNucleus(it->first);
    double f=it->second;
    if(f>0.){
      cout<<setw(32)<<fixed<<setprecision(2)<<100.*f<<setw(0)
          <<" % "<<nuc->getName();
      nuc->printIt();
      sumA    += f*nuc->getA();
      sumMass += f*nuc->getMass();
    }
  }
  cout<<"                           Total"<<setw(9)<<fixed<<setprecision(2)
      <<sumA<<setw(23)<<fixed<<setprecision(1)<<sumMass<<setw(0)<<endl;
  return true;
}


XeGraph* Target::newGraphOfRate(ErRange* range){
  if(!checkIt()) return NULL;
  if(range==NULL) range=XeRange::getDefaultErRange();
  int n=range->getNBins();
  Interaction *inter=getInteraction();
  XeGraph* gr=new XeGraph("Rate for "+getName(),n,LABEL_ER,LABEL_EVT_KG_DAY_KEV
                         ,inter->getStyleFromFF(),inter->getLegendFromFF());
  gr->setName("Differential rate");
  for(int i=0;i<n;i++){
    pair<double,double> e=range->getInterval(i);
    double eMin=e.first;
    double eMax=e.second;
    double r=rate(eMin,eMax)/(eMax-eMin);
    if(getDebugLevel()>1) {
      cout<<"Emin="<<eMin<<", Emax="<<eMax<<", r="<<r<<endl;
    }
    gr->setPoint(i,.5*(eMin+eMax),r);
  }
  if(getDebugLevel()>0) gr->printIt();
  gr->setDefaultYScale(LOG);
  return gr;
}


double Target::getFraction(int n, int a){
  int h=Nucleus::hashCode(n,a);
  if(fractions.find(h)==fractions.end()) return 0.;
  return fractions[h];
}

map<int,double>& Target::getFractions(){return fractions;}

void Target::add(Target* target,double f){
  map<int,double> frac=target->getFractions();
  for(map<int,double>::iterator it=frac.begin();it!=frac.end();it++){
    int hash=it->first;
    double fr=it->second;
    add(Nucleus::nFromHashCode(hash),Nucleus::aFromHashCode(hash),f*fr);
  }
} 

Target::Target(int n, int a) : XePhysics("Target: "+Nucleus::nucleusName(n,a)) {
  init();
  if(a==XE100_MIXTURE) return;
  else if(n==XENON && (a==DEPLETED || a==NATURAL)) {
    if(a==DEPLETED) {
      add(XENON,124,0.0001); 
      add(XENON,126,0.0001); 
      add(XENON,128,0.0063); 
      add(XENON,129,0.2494); 
      add(XENON,130,0.0475); 
      add(XENON,131,0.2575); 
      add(XENON,132,0.3338); 
      add(XENON,134,0.1048); 
      add(XENON,136,0.0007); 
    }
    else {
      add(XENON,124,0.0009); 
      add(XENON,126,0.0009); 
      add(XENON,128,0.0192); 
      add(XENON,129,0.2644); 
      add(XENON,130,0.0408); 
      add(XENON,131,0.2118); 
      add(XENON,132,0.2689); 
      add(XENON,134,0.1044); 
      add(XENON,136,0.0887); 
    }
    Reference="Private communication by G. Plante";
    return;
  }
  else if(Nucleus::isValid(n,a)){
    add(n,a,1.);
    Nucleus *nuc=LedererTable::getInstance()->getNucleus(n,a);
    setName(nuc->getName());
    setFrameName(nuc->getFrameName());
    setLaTeXName(nuc->getLaTeXName());
    return;
  }
  cout<<"Don't know how to make a target with n="<<n<<", a="<<a<<endl;
}

double Target::rate(double Ermin,double Ermax){ // in Events/min/Kg
  if(Ermax==0.) Ermax=ErMax();
  if(checkIt()) return integrate(RATE,Ermin,Ermax);
  return 0.;
}

double Target::getValue(int mode, double x) {
  switch(mode){
    case RATE: double r=dRate(x);
               return r;
  }
  cout<<"Invalid mode :"<<mode<<endl;
  return 0.;
}

double Target::dRate(double Er){  // in Events/day/Kg/KeV
  

  double m_D  = wimp->getMass();
  double v0=galaxyModel->getV_0();
  double rho=galaxyModel->getRho();
  double E_0=.5/KeV*m_D*v0*v0/C/C;    // Wimp kinetic energy in keV

  double rat=0;
  for(map<int,double>::iterator it=fractions.begin();it!=fractions.end();it++){
    int hash=it->first;
    double f=it->second;
    if(f<=0.) continue;
    int n=Nucleus::nFromHashCode(hash);
    int a=Nucleus::aFromHashCode(hash);
    Nucleus *nuc=LedererTable::getInstance()->getNucleus(n,a);
    if(nuc==NULL) continue;
    double sigma0=interaction->SigmaErNucleus(wimp,nuc,Er)/SIGMA_0_SCALE;
    if(sigma0==0.) continue;
    double m_T=nuc->getMass();
    double E0r= 4.*m_D*m_T/(m_D+m_T)/(m_D+m_T)*E_0;
    double vMin=sqrt(Er/E0r)*v0;
    rat+= sigma0*f*(rho/RHO_SCALE) *(v0/V_0_SCALE)*RATE_0/m_D/m_T/E0r
          * galaxyModel->rateScale(vMin);
  }
  return rat;
}

//
//---- Direct parametrisation of structure function ------------------------

StructureFunction::StructureFunction()       : XePhysics() {}       
StructureFunction::StructureFunction(int w)  : XePhysics() {what=w;}
StructureFunction::~StructureFunction()                   {}

XeStyle* StructureFunction::getStyle()          {return new XeStyle();}
void     StructureFunction::setBandEdge(int )   {edge=NO_BAND;}
void     StructureFunction::printCoefficients() {
  cout<<"Print coefficients of "<<getName()<<" not yet implemented"<<endl;
}
string   StructureFunction::getNameWithEdge()   {
  return getNameSpace()+getEdgeName(edge);
}

SIStructureFunction::SIStructureFunction(int w)  : StructureFunction(w) {}
SIStructureFunction::~SIStructureFunction()                             {}

SDStructureFunction* SDStructureFunction::getStructureFunction(int p){
  switch(p){
    case SF_SCHWENK_1BC : 
    case SF_SCHWENK_2BC : return new SchwenkSDStructureFunction(p);
  }
  cout<<"Unknown type of SDStructureFunction: "<<p<<endl;
  return NULL;
}

SDStructureFunction::SDStructureFunction(int p)  : StructureFunction(p) {
  edge=NO_BAND;
  if(p<0){
    p=DEFAULT_SD_SF;
    cout<<"(taking default SD structure function)"<<endl;
  }
  switch(what){
    case SF_SCHWENK_1BC : setName("Schwenk SF (1 bc)")  ; break;
    case SF_SCHWENK_2BC : setName("Schwenk SF (1+2 bc)"); break;
    default          : 
         cout<<"Unknown type of SchwenkSDStructureFunction: "<<what<<endl;
         return;
  }
  setLegend();
  resetCoefficients();
}
SDStructureFunction::~SDStructureFunction()     {}

void   SDStructureFunction::resetCoefficients() {}
double SDStructureFunction::getYMax()           {return yMax;}

XeGraph* SDStructureFunction::newGraphOfSA(Nucleus *nuc,int np,XeRange* yr){
  if(yr==NULL) yr=XeRange::getDefaultYRange();
  int n=yr->getNPoints();
  XeGraph* gr=new XeGraph(
              "SA for "+getNameSpace()+ SDPureInteraction::getModeName(np)
             +" on "+nuc->getNameSpace()+getEdgeName(edge)
             , n,LABEL_Y_UOVER2,LABEL_SA ,getStyle(),getName());
  for(int i=0;i<n;i++){
    double y=yr->getValue(i);
    double sa=SA(nuc,y,np);
    gr->setPoint(i,y,sa);
  }
  gr->setDefaultYScale(LOG);
  return gr;
}

SchwenkSDStructureFunction::~SchwenkSDStructureFunction(){}
SchwenkSDStructureFunction::SchwenkSDStructureFunction(int p)
                           :SDStructureFunction(p){
  resetCoefficients();
}

void SchwenkSDStructureFunction::resetCoefficients(){
  int p=what-FIRST_SD_SF;
  for(int ns=0;ns<N_S_NUCLEI;ns++) {
    for(int np=0;np<N_NUCLEONS;np++) {
      for(int c=0;c<N_POLSD_COEFS;c++) coefs[ns][np][c]=Coefs[p][ns][np][c];
    }
  }
  edge=what==SF_SCHWENK_1BC? NO_BAND:BAND_CENTER;
}

void SchwenkSDStructureFunction::setBandEdge(int e){
  if(what==SF_SCHWENK_1BC) {
    edge=NO_BAND;
    return;
  } 
  edge=e;
  if(edge==BAND_CENTER) resetCoefficients();
  else {
    int p=edge==LOWER_EDGE? 0:1;
    for(int ns=0;ns<N_S_NUCLEI;ns++) {
      for(int np=0;np<N_NUCLEONS;np++) {
        for(int c=0;c<N_POLSD_COEFS;c++) coefs[ns][np][c]=Bands[ns][np][p][c];
      }
    }
  }
}

double SchwenkSDStructureFunction::SA( Nucleus *nuc, double y, int nucleon
                                     , bool prnt) {
  double *c=getCoefficients(nuc->getA(),nucleon);
  if(c==NULL) return 0.;
  double u=2.*y;
  double x=exp(-u);
  double sa=computePolynomial(u,DEGREES_SCHWENK,c)*x;
  if(prnt) {
    cout<<"Computing SA for "<<getName()<<" A="<<nuc->getA()<<", y="<<y
        <<" nucleon="<<nucleon<<" : "<<sa<<endl;   
  }
  return sa;
}

void SchwenkSDStructureFunction::printCoefficients(){
  cout<<endl<<"Coefficients of "<<getName()<<endl;
  int e1=what==SF_SCHWENK_1BC? BAND_CENTER:LOWER_EDGE;
  int e2=what==SF_SCHWENK_1BC? BAND_CENTER:UPPER_EDGE;
  for(int e=e1;e<=e2;e++) {
    setBandEdge(e);
    cout<<getEdgeName(e)
        <<"       Xe 129               Xe 131     " <<endl
        <<"        neutron  proton        neutron  proton "
                                   <<endl;
    for(int p=0;p<=DEGREES_SCHWENK;p++){
      for(int ns=0;ns<N_S_NUCLEI;ns++){ 
        for(int np=0;np<N_NUCLEONS;np++)printf(" %12.4e",coefs[ns][np][p]);
        printf("  ");
      }
      cout<<endl;
    }   
    cout<<endl;
  }
  cout<<endl;
}

double* SchwenkSDStructureFunction::getCoefficients(int A, int nucleon){
  int index=SpinContent::getNucleus(A);
  if(index<0) return NULL;
  return (double*) &coefs[index][nucleon][0];
}

XeStyle* SchwenkSDStructureFunction::getStyle() {
  int col=what==SF_SCHWENK_1BC? kBlack:kBlue;
  return new XeStyle(col);
}

// -------------------  Form Factors ------------------------------------

FormFactor::~FormFactor(){}
FormFactor::FormFactor(): XePhysics(){}
FormFactor::FormFactor(int w): XePhysics(){
  edge=NO_BAND;
  which=w;
}
string   FormFactor::getNameWithEdge()   {
   return getNameSpace()+getEdgeName(edge);
}
void     FormFactor::printCoefficients() {
  cout<<"Print coefficients of "<<getName()<<" not yet implemented"<<endl;
}
XeStyle* FormFactor::getStyle()          {return new XeStyle();}
void     FormFactor::setBandEdge(int ) {edge=NO_BAND;}

SuhonenTable::~SuhonenTable(){}
SuhonenTable::SuhonenTable(){} 

int SuhonenTable::getIndex(double u){
  int i=0;
  if(u<0. || u>=20.) {
    cout<<"Error in SuhonenSDFormFactor::compute: y="<<u/2.<<endl;
    return UNDEFINED_INT;
  }
  if(u<5.) i=(int)(200.*u);
  else     i=(int)(1000+(u-5.)*200./3.);
  if(i==0) i=1;
  else if(i==2000) i-=1; 
  return i;
}

double SuhonenTable::interpolate(double u, int n, int index, int coef) {
  double x[3];
  double s[3];
  x[0]=U_TABLE[n][index-1][0];
  x[1]=U_TABLE[n][index  ][0]; 
  x[2]=U_TABLE[n][index+1][0];
  s[0]=U_TABLE[n][index-1][coef];
  s[1]=U_TABLE[n][index  ][coef];
  s[2]=U_TABLE[n][index+1][coef];
  return XeInterpolation::parabolic(x,s,u);
}
     
// -------------------------- Spin independent form factor ----------------

SIFormFactor::~SIFormFactor(){}
SIFormFactor::SIFormFactor(int w): FormFactor(w){}

SIFormFactor* SIFormFactor::getFormFactor(int w){
  if(w<0){
    cout<<"(taking default SI form factor)"<<endl; 
    w=DEFAULT_SI_FF;
  }
  switch(w) {
    case ENGEL_SI       : return new EngelSIFormFactor();
    case FRICK_SI       : return new FrickSIFormFactor();
    case SUHONEN_SI     : return new SuhonenSIFormFactor();
    case EXPONENTIAL_SI : return new ExponentialFormFactor();
  }
  cout<<"Unknown SI Form factor :"<<w<<endl;
  return NULL;
}

void ExponentialFormFactor::setTheName(){
  char t[100];
  sprintf(t,"Exponential Form Factor, R0=%.2f *R1=%.2f Coherence=%.1f"
         ,R0,R1,Coherence);
  setName(t);
} 
  
ExponentialFormFactor::~ExponentialFormFactor()  {}

ExponentialFormFactor::ExponentialFormFactor(double r0,double r1
                      ,double coherence): SIFormFactor(EXPONENTIAL_SI) {
  R0=r0;
  R1=r1;
  Coherence=coherence;
  setTheName();
  Reference="Jungman et al., Phys. Reports 267 (1996) 195 (eq. 7.31)";
}

double ExponentialFormFactor::computeFF(double q, Nucleus* nuc) {
  int A=nuc->getA();
  double mA= A*ATOMIC_MASS;
  double R= R0 + R1*pow((double)A,1./3.);
  return exp(-q*mA*R*R/2./Coherence*MeV);
}
void ExponentialFormFactor::setR0(double r)        {R0=r; setTheName();}
void ExponentialFormFactor::setR1(double r)        {R1=r; setTheName();}
void ExponentialFormFactor::setCoherence(double c) {Coherence=c; setTheName();}

 
WoodsSaxonFormFactor::~WoodsSaxonFormFactor(){}
WoodsSaxonFormFactor::WoodsSaxonFormFactor(int w): SIFormFactor(w){}
void WoodsSaxonFormFactor::setThickness(int t){
  thickness=t;
  setTheName();
}

double WoodsSaxonFormFactor::computeFFWell(double q, double r2){
  if(q<0.) return 0.;
  else if(q==0.) return 1.;
  double t=q*sqrt(r2);
  return 3*(sin(t)-t*cos(t))/t/t/t;
}


EngelSIFormFactor::~EngelSIFormFactor(){}
EngelSIFormFactor::EngelSIFormFactor(double r0,double s)
                 : WoodsSaxonFormFactor(ENGEL_SI) {
  R0=r0;
  thickness=s;
  setTheName();
  Reference="Engel, Phys Letters B264(1991) 114";
  setLegend("Engel");
}
XeStyle* EngelSIFormFactor::getStyle()     {return new XeStyle(kGreen+1) ;}
void   EngelSIFormFactor::setR0(double r0) {R0=r0; setTheName();}
double EngelSIFormFactor::computeFF(double q, Nucleus* nuc) {
  int A=nuc->getA();
  double r=R0*pow((double)A,1./3.);
  return computeFFWell(q,r*r-5.*thickness*thickness);
}
void EngelSIFormFactor::setTheName(){
  char t[80];
  sprintf(t,"Engel SI form factor, R0=%.2f Fm,  thickness=%.1f"
         ,R0,thickness); 
  setName(t);
}

FrickSIFormFactor::~FrickSIFormFactor(){}
FrickSIFormFactor::FrickSIFormFactor(double r0, double r1 ,double s, double a)
               : WoodsSaxonFormFactor(FRICK_SI) {
  R0=r0;
  R1=r1;
  thickness=s;
  thicknessA=a;
  setTheName();
  setLegend("Frick");
  Reference="Lewin et al., Astro. Physics 6(1996)";
}
XeStyle* FrickSIFormFactor::getStyle()   {return new XeStyle(kBlue);}
void FrickSIFormFactor::setR0(double r0) {R0=r0; setTheName();}
void FrickSIFormFactor::setR1(double r1) {R1=r1; setTheName();}
void FrickSIFormFactor::setThicknessA(double a) {thicknessA=a; setTheName();}

double FrickSIFormFactor::computeFF(double q, Nucleus* nuc) {
  int A=nuc->getA();
  double r=R0*pow((double)A,1./3.)+R1;
  double rn=r*r+7./3.*PI*PI*thicknessA*thicknessA-5.*thickness*thickness;
  return computeFFWell(q,rn);
}
void FrickSIFormFactor::setTheName(){
  char t[120];
  sprintf(t
    ,"Frick FF, R0=%.2f Fm, R1= %.2f Fm  t=%.2f and a=%.2f"
    ,R0,R1,thickness,thicknessA); 
  setName(t);
}

     
SuhonenSIFormFactor::~SuhonenSIFormFactor(){}    
SuhonenSIFormFactor::SuhonenSIFormFactor()
  : SIFormFactor(SUHONEN_SI), SuhonenTable(){
  setName("Suhonen SI");
  setLegend("Suhonen SI");
}    

XeStyle*    SuhonenSIFormFactor::getStyle(){return new XeStyle(32);}
double SuhonenSIFormFactor::computeFF(double y,Nucleus *nuc){
  int A=nuc->getA();
  if(A%2==0) return 0.;
  int n=SpinContent::getNucleus(A);
  if(n<0) return 0.;
  double u=2.*y;
  int i=getIndex(u);
  if(i<0) return 0.;
  return interpolate(u,n,i,2);
}


// -------------------------- Spin dependent form factor ----------------

SDFormFactor::~SDFormFactor(){}
SDFormFactor::SDFormFactor(int w): FormFactor(w) {}

double SDFormFactor::getYMax(){return yMax;}

int SDFormFactor::preferredSpinModel(){
  switch(which) {
    case RESSELL_BONN       : return RESSELL_BONN_SPIN;
    case RESSELL_NIJMEGEN   : return RESSELL_NIJMEGEN_SPIN;
    case FF_SCHWENK_1BC     : 
    case FF_SCHWENK_2BC     : return SCHWENK_SPIN;
    case SUHONEN_SD         : return SUHONEN_SPIN;
  }
  cout<<"Invalid parametrization in "<<getName()<<endl;
  return -1;
}

SDFormFactor* SDFormFactor::getFormFactor(int p){
  if(p<0){
    p=DEFAULT_SD_FF;
    cout<<"(taking default SD form factor)"<<endl;
  }
  switch(p){
    case RESSELL_BONN       : return new RessellBonnSDFormFactor();
    case RESSELL_NIJMEGEN   : return new RessellNijmegenSDFormFactor();
    case FF_SCHWENK_1BC      : return new Schwenk1bcFormFactor();
    case FF_SCHWENK_2BC      : return new Schwenk2bcFormFactor();
    case SUHONEN_SD         : return new SuhonenSDFormFactor();
  }
  cout<<"Can't create SDFormFactor of type "<<p<<endl;
  return NULL;
}

string SDFormFactor::getStructureFactorName(int c){
  switch(c) {
    case S00_SFACT : return "S00" ;
    case S01_SFACT : return "S01" ;
    case S11_SFACT : return "S11" ;
  }
  return UNDEFINED_STRING;
}

string SDFormFactor::getStructureFactorLaTeXName(int c){
  switch(c) {
    case S00_SFACT : return "S_{00}" ;
    case S01_SFACT : return "S_{01}" ;
    case S11_SFACT : return "S_{11}" ;
  }
  return UNDEFINED_STRING;
}

void SDFormFactor::printStructureFactor(double *s){
  for(int c=0;c<N_SFACT;c++) cout<<getStructureFactorName(c)<<":"<<s[c]<<" ";
}

bool SDFormFactor::checkIt(bool){
  Target xenon(XENON,NATURAL);
  bool ok=true;
  for(int xe=129;xe<132;xe+=2) {
    Nucleus*nuc=LedererTable::getInstance()->getNucleus(XENON,xe);
    ok=ok && checkIt(nuc);
  }
  return ok;
}

bool SDFormFactor::checkIt(Nucleus* nuc) {
  saveWarnings();
  bool ok=true;
  for(int p=NEUTRON;p<=PROTON;p++) {
    suppressWarnings();
    for(int i=0;i<=100;i++) {
      double y=i*0.01*yMax;
      if(FF(nuc,y,p)<=0.) {
        showWarnings();
        FF(nuc,y,p,true);       
        ok=false;
        break;
      }
    }
  }
  restoreWarnings();
  return ok;
}

XeGraph* SDFormFactor::newGraphOfFF(Nucleus *nuc, int nucleon, XeRange* yr){
  if(yr==NULL) yr=XeRange::getDefaultYRange();
  int n=yr->getNPoints();
  XeGraph* gr=new XeGraph("FF for "+getName(),n,LABEL_Y_UOVER2,LABEL_FF
            ,getStyle(),getName());
  for(int i=0;i<n;i++){
    double y=yr->getValue(i);
    double f2=FF(nuc,y,nucleon);
    gr->setPoint(i,y,f2);
  }
  gr->setDefaultYScale(LOG);
  return gr;
}

XeGraph* SDFormFactor::newGraphOfSA(Nucleus *nuc, int nucleon, XeRange* yr){
  if(yr==NULL) yr=XeRange::getDefaultYRange();
  int n=yr->getNPoints();
  XeGraph* gr=new XeGraph(
              "SA for "+getNameSpace()+ SDPureInteraction::getModeName(nucleon)
             +" on "+nuc->getNameSpace()+getEdgeName(edge)
             ,n,LABEL_Y_UOVER2,LABEL_SA,getStyle(),getName());

  for(int i=0;i<n;i++){
    double y=yr->getValue(i);
    double sa=SA(nuc,y,nucleon);
    gr->setPoint(i,y,sa);
  }
  gr->setDefaultYScale(LOG);
  return gr;
}

XeGraph* SDFormFactor::newGraphOfStructureFactor(Nucleus *nuc, int fc
                                                , XeRange* yr){
  if(yr==NULL) yr=XeRange::getDefaultYRange();
  int n=yr->getNPoints();
  XeGraph* gr=new XeGraph(
              "SF for "+getNameSpace()+"on "
              +nuc->getNameSpace()+getEdgeName(edge)
              ,n,LABEL_Y_UOVER2,getStructureFactorLaTeXName(fc)
              ,getStyle(),getName());
  for(int i=0;i<n;i++){
    double y=yr->getValue(i);
    double s[N_SFACT];
    computeStructureFactor(nuc,y,s);
    gr->setPoint(i,y,s[fc]);
  }
  gr->setDefaultYScale(LOG);
  return gr;
}

void SDFormFactor::printDetailedStructureFactor(Nucleus *nuc,XeRange* yr){
  cout<<"Detailed calculation of "<<getName()<<" for "<<nuc->getName()<<endl;
  if(yr==NULL) yr=XeRange::getDefaultYRange();
  int n=yr->getNPoints();
  for(int i=0;i<n;i++){
    double y=yr->getValue(i);
    double s[N_SFACT];
    computeStructureFactor(nuc,y,s);
    cout<<y<<" ";
    printStructureFactor(s);
    double sqn=SA(s,1.,0.);
    double sqp=SA(s,0.,1.);
    cout<<" n:"<<sqn<<" p:"<<sqp<<endl;
  }
}

void SDFormFactor::resetCoefficients(){}

double SDFormFactor::FF(Nucleus *nuc, double y, int nucleon,bool prnt){
  return FF(nuc,y,getAn(nucleon),getAp(nucleon),prnt);
}

double SDFormFactor::getAn(int mode){
  switch(mode) {
    case PROTON : return 0.;
    case NEUTRON: return 1.;
  }
  cout<<"Invalid mode "<<mode<<endl;
  return 0.;
}

double SDFormFactor::getAp(int mode){
  switch(mode) {
    case PROTON : return 1.;
    case NEUTRON: return 0.;
  }
  cout<<"Invalid mode "<<mode<<endl;
  return 0.;
}

double SDFormFactor::FF(Nucleus *nuc,double y, double aN, double aP,bool prnt){
  double numerator=SA(nuc,y,aN,aP,prnt); 
  double denominator=SA(nuc,0.,aN,aP,false); 
  if(denominator==0. && doWarn()){
    cout<<"Invalid S(0) for "<<getName()<<endl; 
    return 0.;
  }
  double f2=numerator/denominator;
  if(f2<0. && doWarn()){
    cout<<getName()<<" returns a non-physical result!"<<endl
        <<"FF(y="<<y<<", aN="<<aN<<", aP="<<aP<<") = "<<f2<<endl;
  }
  return f2;
}

double SDFormFactor::SA(double *s, double aN,double aP){
  return (aP+aN)*(aP+aN)*s[S00_SFACT] 
       + (aP-aN)*(aP-aN)*s[S11_SFACT]
       + (aP+aN)*(aP-aN)*s[S01_SFACT]
       ;
}

double SDFormFactor::SA(Nucleus *nuc, double y,  int nucleon){
  double aN=0.;
  double aP=0.;
  if(nucleon==NEUTRON)     aN=1.;
  else if(nucleon==PROTON) aP=1.;
  return SA(nuc,y,aN,aP,false);
}

double SDFormFactor::SA(Nucleus *nuc, double y, double aN,double aP,bool prnt){
  double s[N_SFACT];
  computeStructureFactor(nuc,y,s);
  double sq=SA(s,aN,aP);
  if(prnt){
    cout<<getName()<<" on "<<nuc->getName()<<" y="<<y
        <<" an="<<aN<<" ap="<<aP<<endl;
    printStructureFactor(s);
    cout<<" sa="<<sq<<endl;
  }
  return sq;
}

PolynomialSDFormFactor::~PolynomialSDFormFactor(){}
PolynomialSDFormFactor::PolynomialSDFormFactor(int p) : SDFormFactor(p){
  switch(which) {
    case RESSELL_BONN :
         setParameters("Ressell/Bonn" ,DEGREES_RESSELL,10.); 
         break;
    case RESSELL_NIJMEGEN :
         setParameters("Resell/Nijmegen",DEGREES_RESSELL,10.);
         break;
    case FF_SCHWENK_1BC :
         setParameters("Schwenk (1 bc)" ,DEGREES_SCHWENK,10.);
         break;
    case FF_SCHWENK_2BC :
         setParameters("Schwenk (1+2 bc)" ,DEGREES_SCHWENK,10.);
         break;
  }
  resetCoefficients();
}

void PolynomialSDFormFactor::resetCoefficients(){
  for(int n=0;n<N_S_NUCLEI;n++){
    for(int f=0;f<N_SFACT;f++){
      for(int c=0;c<N_POLSD_COEFS;c++) coefs[n][f][c]=Coefs[which][n][f][c];
    }
  }
  edge=which==FF_SCHWENK_2BC? BAND_CENTER:NO_BAND ;
}

void PolynomialSDFormFactor::setParameters(string nam, int n, double y){
  yMax=y;
  nDegrees=n;
  setName(nam);
  setLegend(nam);
  if(which==FF_SCHWENK_1BC || which==FF_SCHWENK_2BC)  {
    Reference="Achim Schwenk et al, http://xxx.lanl.gov/abs/1208.1094";
  }
  else {
    Reference="Ressel & Dean, Phys Rev C, 56-1 (1997) 535";
  }
}

void PolynomialSDFormFactor::printCoefficients(){
  cout<<endl<<"Coefficients of "<<getName()<<endl;
  for(int ns=0;ns<N_S_NUCLEI;ns++) {
    cout<<"                      "<<SpinContent::getNucleusName(ns)<<"      ";
  }
  cout<<endl<<"Power ";
  for(int n=0;n<N_S_NUCLEI;n++) {
    cout<<"   S_00        S_01        S_11      ";
  }
  cout<<endl;
  for(int p=0;p<=nDegrees;p++){
    printf("%3d",p);
    for(int n=0;n<N_S_NUCLEI;n++) {
      printf(" ");
      for(int s=0;s<N_SFACT;s++) printf("%12.4e",coefs[n][s][p]); 
    }
    cout<<endl;
  }
}

double* PolynomialSDFormFactor::getCoefficients(int A, int w){
  int index=SpinContent::getNucleus(A);
  if(index<0) return NULL;
  return (double*) &coefs[index][w][0];
}

void PolynomialSDFormFactor::computeStructureFactor(Nucleus *nuc,double y
                                                   ,double* s){

  int A=nuc->getA();
  if(A<0) {
    cout<<getName()<<" can't handle nucleus "<<nuc->getName()<<endl; 
    for(int i=0;i<N_SFACT;i++) s[i]=0.;
    return;
  }

  double *c00=getCoefficients(A,S00_SFACT);
  double *c01=getCoefficients(A,S01_SFACT);
  double *c11=getCoefficients(A,S11_SFACT);
  double u=2.*y;
  double x=exp(-u);
  switch(which){
    case FF_SCHWENK_1BC: 
    case FF_SCHWENK_2BC: 
          s[S00_SFACT]=computePolynomial(u,nDegrees,c00)*x;
          s[S01_SFACT]=computePolynomial(u,nDegrees,c01)*x;
          s[S11_SFACT]=computePolynomial(u,nDegrees,c11)*x;
          break;
         
    case RESSELL_BONN    :
    case RESSELL_NIJMEGEN:
         {
          double r=1./(1.+y);
          s[S00_SFACT]=(computePolynomial(y,nDegrees-1,c00)+c00[nDegrees]*r)*x;
          s[S01_SFACT]=(computePolynomial(y,nDegrees-1,c01)+c01[nDegrees]*r)*x;
          s[S11_SFACT]=(computePolynomial(y,nDegrees-1,c11)+c11[nDegrees]*r)*x;
          break;
         }
   }
}


RessellBonnSDFormFactor::~RessellBonnSDFormFactor(){}
RessellBonnSDFormFactor::RessellBonnSDFormFactor()
 : PolynomialSDFormFactor(RESSELL_BONN) {}
XeStyle* RessellBonnSDFormFactor::getStyle(){return new XeStyle(kGreen+1);}

XeStyle* RessellNijmegenSDFormFactor::getStyle(){return new XeStyle(kRed);}
RessellNijmegenSDFormFactor::~RessellNijmegenSDFormFactor(){}
RessellNijmegenSDFormFactor::RessellNijmegenSDFormFactor()
 : PolynomialSDFormFactor(RESSELL_NIJMEGEN) {}

Schwenk1bcFormFactor::~Schwenk1bcFormFactor(){}
Schwenk1bcFormFactor::Schwenk1bcFormFactor() 
  : PolynomialSDFormFactor(FF_SCHWENK_1BC) {}
XeStyle* Schwenk1bcFormFactor::getStyle(){return new XeStyle(kBlack,2,2);}


Schwenk2bcFormFactor::~Schwenk2bcFormFactor(){}
Schwenk2bcFormFactor::Schwenk2bcFormFactor() 
  : PolynomialSDFormFactor(FF_SCHWENK_2BC) {}
XeStyle* Schwenk2bcFormFactor::getStyle(){return new XeStyle(kBlue,2,2);}

void Schwenk2bcFormFactor::setBandEdge(int e){
  resetCoefficients();  
  edge=e;
  if(edge==BAND_CENTER) return;
  for(int b=0;b<2;b++){  
    int q=0;
    int w=0;
    if(b==0) { 
      w=S01_SFACT; 
      q=edge==LOWER_EDGE? S01_LOW:S01_UPP;
    }
    else { 
      w=S11_SFACT; 
      q=edge==LOWER_EDGE? S11_LOW:S11_UPP;
    }
    for(int n=0;n<N_S_NUCLEI;n++){
      for(int c=0;c<N_POLSD_COEFS;c++) coefs[n][w][c]=Bands[n][q][c];  
    }
  }
 
}

SchwenkYRange::SchwenkYRange(int nb)
  : GeneralRange("Y",nb,SchwenkDeltaA1::yValues){}

SchwenkYRange::~SchwenkYRange(){}

SchwenkDeltaA1::~SchwenkDeltaA1(){}
SchwenkDeltaA1::SchwenkDeltaA1() :XePhysics() {}

void SchwenkDeltaA1::getSA(double *sa, Nucleus* nuc, int nucleon, int delta_a1
                          ,int norm) {
  SDStructureFunction *sf = (norm==UNDEFINED_INT) ? NULL
                          : SDStructureFunction::getStructureFunction(norm);
  SchwenkYRange syr;
  double* y=syr.getValues();
  for(int i=0;i<N_SCHWENK_Y;i++)  sa[i]=UNDEFINED;
  int n=SpinContent::getNucleus(nuc->getA());
  if(n<0) return;
  if(nucleon!=PROTON && nucleon!=NEUTRON) return;
  if(delta_a1<0 || delta_a1 >=N_SCHWENK_DELTA_A1) return;
  for(int i=0;i<N_SCHWENK_Y;i++) {
    sa[i]=SA_NP[n][nucleon][i][delta_a1];
    if(sa[i]!=UNDEFINED && sf!=NULL) sa[i]/=sf->SA(nuc,y[i],nucleon);    
  }
}

XeGraph* SchwenkDeltaA1::newGraphOfSA(Nucleus*nuc,int nucleon,int delta_a1
                                     ,int normalizedTo) {
  int n=SpinContent::getNucleus(nuc->getA());
  if(n<0) return NULL;
  if(nucleon!=PROTON && nucleon!=NEUTRON) return NULL;
  if(delta_a1<0 || delta_a1 >=N_SCHWENK_DELTA_A1) return NULL;
  double sa[N_SCHWENK_Y];
  getSA(sa,nuc,nucleon,delta_a1,normalizedTo);
  SchwenkYRange syr;
  double* y=syr.getValues();
  XeGraph* gr=new XeGraphPositive("SA for SchendDeltaA1",N_SCHWENK_Y,y,sa
                    ,LABEL_Y_UOVER2,LABEL_SA,getStyle(),"Schwenk data");
  gr->setDefaultYScale(LOG);
  return gr;
}
  
XeStyle* SchwenkDeltaA1::getStyle(){return new XeStyle(32,1,1);}

// --------------- Suhonen ------------------------------
SuhonenSDFormFactor::~SuhonenSDFormFactor() {}
SuhonenSDFormFactor::SuhonenSDFormFactor() : 
  SDFormFactor(SUHONEN_SD),SuhonenTable() {
  setName("Suhonen");
  setLegend("Suhonen");
  Reference="Holmlund et al, Phys. Letters B 584 (2004) 31";
  yMax=10.;
}
XeStyle* SuhonenSDFormFactor::getStyle(){return new XeStyle(32);}


void SuhonenSDFormFactor::computeStructureFactor(Nucleus *nuc,double y
                                                , double *s){ 
  for(int i=0;i<N_SFACT;i++) s[i]=0.;
  int n=SpinContent::getNucleus(nuc->getA());
  if(n<0) return;
  double factor=(1.+2.*nuc->getJ())/8./PI;

  if(y==0.){
    s[S00_SFACT]=factor/2.*OMEGA_0[n]*OMEGA_0[n];
    s[S01_SFACT]=factor*OMEGA_0[n]*OMEGA_1[n];
    s[S11_SFACT]=factor/2.*OMEGA_1[n]*OMEGA_1[n];
    return; 
  }       

  // find the proper line
  double u=2.*y;
  int i=getIndex(u);
  if(i<0) return;
  s[S00_SFACT]=interpolate(u,n,i,2+S00_SFACT)*factor/2.*OMEGA_0[n]*OMEGA_0[n];
  s[S01_SFACT]=interpolate(u,n,i,2+S01_SFACT)*factor*OMEGA_0[n]*OMEGA_1[n];
  s[S11_SFACT]=interpolate(u,n,i,2+S11_SFACT)*factor/2.*OMEGA_1[n]*OMEGA_1[n];
}
    
#include "XePhys.data"

// -----------------------Spin Content --------------------------------------

SpinContent* SpinContent::instances[N_SPIN_CONTENT]={NULL,NULL,NULL,NULL};
int          SpinContent::NUCLEI[N_S_NUCLEI]={129,131};

double SpinContent::getProtonSpin(int A)  {return S_PROTON[getNucleus(A)];}
double SpinContent::getNeutronSpin(int A) {return S_NEUTRON[getNucleus(A)];}

double SpinContent::getSpin(int nucleon, int A){
  switch(nucleon) {
    case PROTON  : return getProtonSpin(A);
    case NEUTRON : return getNeutronSpin(A);
  }
  cout<<"Invalid nucleon type in SpinContent::getSpin("<<nucleon
       <<","<<A<<")"<<endl; 
  return 0;
}

SpinContent::~SpinContent(){}
SpinContent::SpinContent(int m) : XePhysics(getTheName(m)){
  model=m;
  for(int i=0;i<N_S_NUCLEI;i++) {
    S_PROTON[i]=SPIN_C[PROTON][model][i];
    S_NEUTRON[i]=SPIN_C[NEUTRON][model][i];
  }
}

string SpinContent::getNucleusName(int n){
  switch(n){ 
    case XENON_129 : return "Xe129";
    case XENON_131 : return "Xe131";
  }
  return UNDEFINED_STRING;
}

string SpinContent::getTheName(int mode){
  switch(mode){
    case RESSELL_BONN_SPIN     : return "Ressel-Bonn";
    case RESSELL_NIJMEGEN_SPIN : return "Ressel-Nijmegen";
    case SCHWENK_SPIN          : return "Schwenk";
    case SUHONEN_SPIN          : return "Suhonen";
  }
  return UNDEFINED_STRING;
}

SpinContent* SpinContent::getInstance(int m){
  if(m<0 || m>=N_SPIN_CONTENT) {
    cout<<"Wrong model for spin content:"<<m<<endl;
    return NULL;
  }
  if(instances[m]==NULL) instances[m]=new SpinContent(m);
  return instances[m];
}

int SpinContent::getNucleus(int A){
  for(int i=0;i<N_S_NUCLEI;i++) if(NUCLEI[i]==A) return i;
  cout<<"Can't handle A="<<A<<endl;
  return UNDEFINED_INT;
}

// ----------------------- Interaction --------------------------------------


Interaction::~Interaction(){}
Interaction::Interaction(double s) : XePhysics("Interaction") {
  sigmaNucleon=s;
}

void   Interaction::setSigmaNucleon(double s) {
  if(sigmaNucleon==s) return;
  sigmaNucleon=s;
  markAsChanged();
}

bool               Interaction::isSpinDependent()      {return spinDependent;}
bool               Interaction::withStructureFunction(){return structure!=NULL;}
bool               Interaction::withFormFactor()       {return form!=NULL;}
double             Interaction::getSigmaNucleon()      {return sigmaNucleon;}
SpinContent*       Interaction::getSpinContent()       {return NULL;}
StructureFunction* Interaction::getStructureFunction() {return structure;}
FormFactor*        Interaction::getFormFactor()        {return form;}
void               Interaction::setFormFactor(FormFactor* f) {form=f;}

double Interaction::SigmaEr(Wimp* wimp, Target* target, double Er) {
  double s=0;
  map<int,double> frac=target->getFractions();
  for(map<int,double>::iterator it=frac.begin();it!=frac.end();it++){
    int hash=it->first;
    double f=it->second;
    if(f<=0.) continue;
    int n=Nucleus::nFromHashCode(hash);
    int a=Nucleus::aFromHashCode(hash);
    Nucleus *nuc=LedererTable::getInstance()->getNucleus(n,a);
    s += SigmaErNucleus(wimp,nuc,Er)*f;
  }
  return s;
}

void Interaction::updateKinematics(Wimp* wimp, Nucleus* nuc){
  A     = nuc->getA();
  m_A   = nuc->getMass();
  m_D   = wimp->getMass();
  mu_A  = m_D*m_A/(m_D+m_A);
  mu_p  = m_D*PROTON_MASS/(PROTON_MASS+m_D); 
  mu_n  = m_D*NEUTRON_MASS/(NEUTRON_MASS+m_D); 
  kin_p = (mu_A/mu_p)*(mu_A/mu_p);
  kin_n = (mu_A/mu_n)*(mu_A/mu_n);
}


// returns max Q in GeV/c
double Interaction::QMax(Wimp* wimp,double V, Target* target){
  double v=V/C;
  double s=0.;
  map<int,double> frac=target->getFractions();
  for(map<int,double>::iterator it=frac.begin();it!=frac.end();it++){
    int hash=it->first;
    double f=it->second;
    if(f<=0.) continue;
    Nucleus *nuc=LedererTable::getInstance()->getNucleus(hash);
    updateKinematics(wimp, nuc);
    s=max(s,2.*mu_A*v);
  }
  return s;
}

// returns max recoil energy in KeV/c^2
double Interaction::ErMax(Wimp* wimp, double V, Target* target){
  double v=V/C;
  double s=0.;
  map<int,double> frac=target->getFractions();
  for(map<int,double>::iterator it=frac.begin();it!=frac.end();it++){
    int hash=it->first;
    double f=it->second;
    if(f<=0.) continue;
    Nucleus* nuc=LedererTable::getInstance()->getNucleus(hash);
    updateKinematics(wimp,nuc);
    double er=nuc->QtoEr(2.*mu_A*v);
    s=max(s,er);
  }
  if(s<=0.) {
    cout<<"Problem in computing Ermax, target was: "<<endl;
    target->printIt();
  }
  return s;
}

// returns max y
double Interaction::YMax(Wimp* wimp, double V, Target* target){
  double v=V/C;
  double s=0.;
  map<int,double> frac=target->getFractions();
  for(map<int,double>::iterator it=frac.begin();it!=frac.end();it++){
    int hash=it->first;
    double f=it->second;
    if(f<=0.) continue;
    Nucleus* nuc=LedererTable::getInstance()->getNucleus(hash);
    updateKinematics(wimp,nuc);
    s=max(s,nuc->QtoY(2.*mu_A*v));
  }
  return s;
}

XeStyle* Interaction::getStyleFromFF() {
  if(form!=NULL)      return form->getStyle(); 
  if(structure!=NULL) return structure->getStyle(); 
  return NULL;
}

string Interaction::getLegendFromFF() {
  if(form!=NULL)           return form->getLegend(); 
  else if(structure!=NULL) return structure->getLegend(); 
  return "";
}

string Interaction::getNameWithFF() {
  string res=getNameSpace();
  if(form!=NULL)           res+=", "+form->getNameWithEdge();
  else if(structure!=NULL) res+=", "+structure->getNameWithEdge(); 
  return res;
}

SIInteraction::~SIInteraction(){}
SIInteraction::SIInteraction(double s, SIFormFactor*f) : Interaction(s) {
  if(f==NULL){
    f=SIFormFactor::getFormFactor(DEFAULT_SI_FF);
    cout<<"(taking default form factor)"<<endl;
  }
  setParameters(f);
}

void SIInteraction::setParameters(SIFormFactor *f){
  spinDependent=false;
  structure=NULL;
  form=f;
  char n[80];
  sprintf(n,"Interaction: Spin independent, sigma_0=%.2e cm^2"
         ,sigmaNucleon);
  setName(n);
  setLegend(getLegendFromFF());
}

SIInteraction::SIInteraction(double s, int w) : Interaction(s) {
  if(w<0){
    cout<<"(taking default form factor for SI interaction)"<<endl;
    w=DEFAULT_SI_FF;
  }
  setParameters(SIFormFactor::getFormFactor(w));
}

SIStructureFunction* SIInteraction::getSIStructureFunction() {
  return (SIStructureFunction*)structure;
}
SIFormFactor* SIInteraction::getSIFormFactor() {
  return (SIFormFactor*) form;
}

double SIInteraction::SigmaErNucleus(Wimp* wimp, Nucleus* nuc,double Er){
  updateKinematics(wimp,nuc);
  double q = nuc->ErtoQ(Er)/HBAR_C_FM;
  double FF_SI=((SIFormFactor*)form)->computeFF(q,nuc);
  return sigmaNucleon*kin_n*(A*FF_SI)*(A*FF_SI);
}

SDInteraction::~SDInteraction(){}
SDInteraction::SDInteraction(double s,double an, double ap, SDFormFactor* f)
  : Interaction(s){setParameters(an,ap,f,NULL);}

SDInteraction::SDInteraction(double s,double an,double ap
  ,SDStructureFunction* sf) : Interaction(s) {setParameters(an,ap,NULL,sf);}


void SDInteraction::setParameters(double an, double ap,SDFormFactor* f
                                      ,SDStructureFunction* sf){
  aN=an;
  aP=ap;
  nucleon=UNDEFINED_INT; 
  spinDependent=true;
  if(sf!=NULL) {
    structure=sf;
    form=NULL;
    spinContent=NULL;
  }
  if(f!=NULL){
    structure=NULL;
    form=f;
    setSpinContentModel(((SDFormFactor*)form)->preferredSpinModel());
  }
  setLegend(getLegendFromFF());
}


SDInteraction::SDInteraction(double s,double an,double ap,int what)
  : Interaction(s){
  if(what<0){
    what=DEFAULT_SD_INTERACTION;
    cout<<"(taking default SD form factor / structure function)"<<endl;
  }
  if(what<N_SD_FF) {
    setParameters(an,ap,SDFormFactor::getFormFactor(what),NULL);
  }
  else {
    setParameters(an,ap,NULL,SDStructureFunction::getStructureFunction(what));
  }
}

SDStructureFunction* SDInteraction::getSDStructureFunction() {
  return (SDStructureFunction*)structure;
}
SDFormFactor* SDInteraction::getSDFormFactor() {return (SDFormFactor*) form;}
SpinContent*  SDInteraction::getSpinContent()  {return spinContent;}

void SDInteraction::setSpinContentModel(int model){
  spinContent=SpinContent::getInstance(model);
}

void SDInteraction::setFormFactor(SDFormFactor *f){
  form=f;
  structure=NULL;
}

void SDInteraction::setStructureFunction(SDStructureFunction *s) {
  structure=s;
  form=NULL;
}

void SDInteraction::setBandEdge(int edge){
  if(withFormFactor())             form->setBandEdge(edge);
  else if(withStructureFunction()) structure->setBandEdge(edge);
  else cout<<"Neither structure function nor form factor in"<<getName()<<endl; 
}

double SDInteraction::SigmaErNucleus(Wimp* wimp,Nucleus* nuc,double Er){
  if(nucleon!=PROTON && nucleon!=NEUTRON){
    cout<<getName()<<" isn't pure neutron or proton"<<endl;
    return UNDEFINED;
  }
  updateKinematics(wimp,nuc);
  if(A%2==0) return 0.;
  double y = nuc->ErtoY(Er);
  double J = nuc->getJ();
  double spinFactor=1.;
  double kin=nucleon==NEUTRON? kin_n:kin_p;
  if(withFormFactor()){
    double spin= nucleon==NEUTRON ? spinContent->getNeutronSpin(A)
                                  : spinContent->getProtonSpin(A);
    spinFactor=((SDFormFactor*)form)->FF(nuc,y,aN,aP)*spin*spin*(J+1)/J;
  }
  else if(withStructureFunction()) { 
    spinFactor=((SDStructureFunction*)structure)->SA(nuc,y,nucleon)*PI/(2*J+1);
  }
  else {
    cout<<"Neither structure function nor form factor in"<<getName()<<endl; 
    return UNDEFINED;
  }
  return 4./3.*sigmaNucleon*kin*spinFactor;
}

SDPureInteraction::~SDPureInteraction(){}

SDPureInteraction::SDPureInteraction(double s, int m , SDFormFactor*f) 
  : SDInteraction(s,SDFormFactor::getAn(m),SDFormFactor::getAp(m),f) {
  setNucleon(m);
}

SDPureInteraction::SDPureInteraction(double s, int m , SDStructureFunction *f) 
  : SDInteraction(s,SDFormFactor::getAn(m),SDFormFactor::getAp(m),f) {
  setNucleon(m);
}

SDPureInteraction::SDPureInteraction(double s, int m , int what) 
  : SDInteraction(s,SDFormFactor::getAn(m),SDFormFactor::getAp(m),what) {
  setNucleon(m);
}

void SDPureInteraction::setNucleon(int m){
  nucleon=m;
  string nn=getModeName(nucleon);
  setName("SD "+getModeName(nucleon)+", sigma_0="
               +formatG(sigmaNucleon,2,6));
} 

string SDPureInteraction::getModeName(int m){
  switch(m) {
    case NEUTRON: return "pure neutron";
    case PROTON : return "pure proton"; 
  }
  cout<<"Invalid nucleon for SD pure interaction"<<endl;
  return UNDEFINED_STRING;
}

string SDPureInteraction::getModeName(){return getModeName(nucleon);}



SDPureProtonInteraction::~SDPureProtonInteraction(){}
SDPureProtonInteraction::SDPureProtonInteraction(double s,SDFormFactor *f)
   :SDPureInteraction(s,PROTON,f) {}
SDPureProtonInteraction::SDPureProtonInteraction(double s
    ,SDStructureFunction *f) :SDPureInteraction(s,PROTON,f) {}

SDPureNeutronInteraction::~SDPureNeutronInteraction(){}
SDPureNeutronInteraction::SDPureNeutronInteraction(double s,SDFormFactor *f)
   :SDPureInteraction(s,NEUTRON,f) {}
SDPureNeutronInteraction::SDPureNeutronInteraction(double s
   ,SDStructureFunction *f) :SDPureInteraction(s,NEUTRON,f) {}


