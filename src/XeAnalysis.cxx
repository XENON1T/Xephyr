#include "XeAnalysis.h"
#include "XeMath.h"
#include "XeStat.h"


//-------------- Stuff linked to a run ------------------------------

RunComponent::RunComponent() : XeObject(){run=NULL;}
RunComponent::~RunComponent(){}
RunComponent::RunComponent(string n) : XeObject(){run=NULL;setName(n);}

XeRun* RunComponent::getRun()             {return run;}
int    RunComponent::getRunNumber()       {return runNumber;}
bool   RunComponent::isRunCompatible(int) {return true;}

void   RunComponent::traceTheFlags()      {
  if(run==NULL) cout<<"no run defined yet"<<endl;
  run->traceTheFlags();
}

 
void RunComponent::setRun(XeRun* r){
  run=r;
  runNumber=run->getNumber();
  setExperiment(runNumber);
  run->attach(this);
  checkRunCompatibility();
}

bool RunComponent::checkRunCompatibility(bool warn){
  bool ok=isRunCompatible(runNumber);
  if(!ok && warn){
    cout<<"Warning, "<<getName()<<" isn't compatible with run "<<runNumber
        <<endl;
  }
  return ok;
}
 
//------------- Stuff linked to a detector -------------------------

DetectorComponent::DetectorComponent() : XeObject(){}
DetectorComponent::DetectorComponent(string n) : XeObject(){setName(n);}
DetectorComponent::~DetectorComponent(){}

//--------------------- Qy --------------------------------------

Qy::Qy() :DetectorComponent()  {
  tValue=0.;
  maxTval = 1.;  //Harmonize this for the two signal model... LEFF_TVALUE_MAX is a const!!  
  TvalStepSize = 0.2;

}

 
// -------------------- LEFF --------------------------------------
int              LEff::tValueMode=LINEAR;
LinearRange*     LEff::tRange=new LinearRange("LEff t-value",N_LEFF_TABULATED
                                             ,-LEFF_TVALUE_MAX,LEFF_TVALUE_MAX);
LinearRange*     LEff::getTRange()          {return tRange;}
void             LEff::setTValueMode(int m) {
  switch(m) {
   case LINEAR :
   case EXPONENTIAL :
        tValueMode=m;
        return;
  }
  cout<<"Invalid LEff t-value handling mode"<<endl;
}
pair<int,double> LEff::getBinAndFraction(double t){
  return tRange->getIndexAndFraction(t);
}
LEff::LEff() :DetectorComponent()  {
  tValue=0.;
  ErMin=DEFAULT_LEFF_ER_MIN;
  maxTval = LEFF_TVALUE_MAX;  //Harmonize this for the two signal model... LEFF_TVALUE_MAX is a const!!  
}

LEff::~LEff() {}

double LEff::getTValue()         {return tValue;}
void   LEff::setTValue(double t) {
  int level=getDebugLevel();
  if(level>1) {
    cout<<getName()<<": t-value set to "<<t<<" ;  ";
  }
  if(tValue==t) {
    if(level>1) cout<<"no need to mark it as changed"<<endl;
    return;
  }
  tValue=t;
  if(level>1) cout<<"mark it as changed"<<endl;
  markAsChanged(false);
}

double LEff::getErMin()         {return ErMin;}
void   LEff::setErMin(double e) {
  if(getDebugLevel()>1) {
    cout<<getName()<<": ErMin set to "<<e<<endl;
  }
  if(ErMin==e) return;
  ErMin=e;
  markAsChanged();
}


double LEff::getLEff(double Er)  {return getLEff(Er,tValue);}

double LEff::getLEff(double Er,double t){
 if(Er<ErMin) return 0.;
 double l=0;
  if(t>0.)      l=XeInterpolation::interpolate(tValueMode,t
                           ,getLEffFromTable(Er,CENTRAL)
                           ,getLEffFromTable(Er,ONE_SIGMA_ABOVE));
  else if(t<0.) l=XeInterpolation::interpolate(tValueMode,-t
                           ,getLEffFromTable(Er,CENTRAL)
                           ,getLEffFromTable(Er,ONE_SIGMA_BELOW));
  else          l=getLEffFromTable(Er,CENTRAL);
  
  return max(0.,l);
}

XeGraph*  LEff::newGraphOfCurrentLEff(ErRange *er,int plot){
  return newGraphOfLEff(er,tValue,plot);
}

XeGraph*  LEff::newGraphOfLEff(ErRange *er,double t,int plot){
  if(er==NULL) er=XeRange::getDefaultErRange();
  int n=er->getNPoints();
  XeGraph *g=new XeGraph("LEff for "+getName(),n,LABEL_ER,LABEL_LEFF);
  for(int i=0;i<n;i++){
    double e=er->getValue(i);
    double leff=getLEff(e,t);
    g->setPoint(i,e,leff);
  }
  g->setDefaultYScale(LINEAR);
  g->drawIt(plot);
  return g;
}

double LEff::tabulated(int T_INDEX){
  return (T_INDEX-SIDE_LEFF_TABULATED)*LEFF_TVALUE_STEP;
}

string LEff::tabulatedName(int t) {return tValueName(tabulated(t));}

XeMultiGraph* LEff::newMultiGraphOfLEff(ErRange *er,int step,int plot){
  XeMultiGraph *mg=new XeMultiGraph("LEff values",LABEL_ER,LABEL_LEFF
                                   ,LABEL_LEFF_TVALUE);
  for(int i=0;i<N_LEFF_TABULATED;i+=step){
    double le=tabulated(i);
    XeGraph *g=newGraphOfLEff(er,le);
    g->setTValueLegend(le);
    g->setLineStyleColor(i-SIDE_LEFF_TABULATED);
    mg->add(g,le);
  }
  mg->setDefaultYScale(LINEAR);
  mg->drawIt(plot);
  return mg;
}

LEffRuns8And10::LEffRuns8And10()  : LEff() {
  setName("LEff for runs 8 and 10");
  setXenon100Reference("run10ubp:leff");
}
LEffRuns8And10::~LEffRuns8And10() {}

double LEffRuns8And10::getLEffFromTable(double Er,int mode)  {
  return XeInterpolation::linear(_er,&(_leff[mode][0]),N_LEFF_RUNS_8_AND_10,Er);
}
 
LEff2011::LEff2011()  : LEff() {
  setName("LEff as of end 2011");
  setReference("arXiv:1104.2587");
}
LEff2011::~LEff2011() {}


double LEff2011::getLEffFromTable(double Er,int mode)  {
  double x[3];
  double y[3];
  int index=(int)(Er/ER_STEP);
  if(index<1) return 0.;
  else if(index>=N_ER_POINTS-1) return _leff[N_ER_POINTS-1][mode];
  x[1]=index*ER_STEP;
  x[0]=x[1]+ER_STEP;
  x[2]=x[1]-ER_STEP;
  y[0]=_leff[index-1][mode];
  y[1]=_leff[index  ][mode];
  y[2]=_leff[index+1][mode];
  double e=XeInterpolation::parabolic(x,y,Er);
  return max(e,0.);
}

LEff* LEff::newDefault(int runNumber){
  switch (runNumber){
    case RUN_08:
    case RUN_10:  return new LEffRuns8And10(); 
    case RUN_12:  return new LEffRuns8And10(); 
  }
  cout<<"No default LEff implemented for run "<<runNumber<<endl;
  return NULL;
}

//--------------- General cut , both HW and Anaylsis ---------------------

XeCut::~XeCut(){}
XeCut::XeCut() : S1S2Object(), RunComponent(){}
XeCut::XeCut(XeRun* r,string name) : S1S2Object(), RunComponent()  {
  setName(name);
  cutMode=CUT_UNKNOWN;
  nature=USER_GENERAL_CUT;
  enabled=true;
  cutFunction=NULL;
  overallAcceptance=1.;
  setRun(r);
  markAsChanged();
}

void     XeCut::enable(bool ena)            {enabled=ena;markAsChanged();}
int      XeCut::getNature()                 {return nature;}
bool     XeCut::passIt(double,double)       {return true;}   // /!\  this never cuts so --> S1S2Data[_cut_] = S1S2Data[]
bool     XeCut::passes(double S1,double S2) {return enabled?passIt(S1,S2):true;}
bool     XeCut::isEnabled()                 {return enabled;}
double   XeCut::getOverallAcceptance()      {return overallAcceptance;}
double   XeCut::maxX()                      {return run->getAnalysisS1Min();}
double   XeCut::minX()                      {return run->getAnalysisS1Max();}
XeGraph* XeCut::newGraph(int )              {return NULL;}
TF1*     XeCut::getCutFunction()            {return cutFunction;}

double   XeCut::maxY()  {
  return S1S2Display::getY(maxX(),maxX()*run->getHighestS2OverS1());
}

double   XeCut::minY()  {
  return S1S2Display::getY(minX(),minX()*run->getLowestS2OverS1());
}

void XeCut::draw() {
  if(S1S2Display::getMode()<0) return;
  XeGraph* g=newGraph();
  if(g==NULL) return;
  int n=g->getN();
  if(debugLevel>1){
    cout<<"Drawing the cut "<<getName()<<", a graph of "<<n<<" points"<<endl;
    g->printIt(10);
  }
  g->setLineColor(CUT_LINE_COLOR);
  if(S1S2Display::areCutsHatched()){
    if(cutMode!=CUT_UNKNOWN){
      int w=cutMode==CUT_ABOVE? -1:1;
      g->setLineWidth(w*CUT_FILL_WIDTH);  
      g->setFillProperties(CUT_FILL_COLOR,CUT_FILL_STYLE);  
    }
  }
  if(n>1) g->drawInsideS1S2("L");
}

// ------------------ Set of Xe cuts --------------------------

XeSetOfCuts::~XeSetOfCuts()  {reset();}
XeSetOfCuts::XeSetOfCuts() : S1S2Object(), vector<XeCut*>(), RunComponent() {
  setName("No cuts");
  init();
}

void   XeSetOfCuts::init()                 {run=NULL;runNumber=UNDEFINED_INT;}
int    XeSetOfCuts::getNCuts()             {return size();}
XeCut* XeSetOfCuts::getCutBySequence(int i){
  return checkIndex(i)?(*this)[i]:NULL;
}

bool XeSetOfCuts::isRunCompatible(int r) {
  int n=getNCuts();
  for(int i=0;i<n;i++) {
    XeCut* cut=getCutBySequence(i); 
    if(!cut->isRunCompatible(r)) return false;
  }
  return true;
}

XeCut* XeSetOfCuts::getCutByNature(int nat){
  int n=getNCuts();
  for(int i=0;i<n;i++) {
    XeCut* cut=getCutBySequence(i); 
    if(cut==NULL) continue;
    if(cut->getNature()==nat) return cut;
  }
  cout<<"Could not find cut of given nature"<<endl;
  return NULL;
}

bool XeSetOfCuts::checkIndex(int i){
  if(i<0 || i>=getNCuts()) {
    cout<<"Invalid Cut number "<<i<<endl;
    return false;
  }
  return true;
}

double XeSetOfCuts::minX(){
  double r=VERY_LARGE;
  int n=getNCuts();
  for(int i=0;i<n;i++) {
    XeCut* cut=getCutBySequence(i);
    if(cut==NULL) continue;
    r=min(r,cut->minX());
  }
  return r;
}

double XeSetOfCuts::minY(){
  double r=VERY_LARGE;
  int n=getNCuts();
  for(int i=0;i<n;i++) {
    XeCut* cut=getCutBySequence(i);
    if(cut==NULL) continue;
    r=min(r,cut->minY());
  }
  return r;
}

double XeSetOfCuts::maxX(){
  double r=VERY_SMALL;
  int n=getNCuts();
  for(int i=0;i<n;i++) {
    XeCut* cut=getCutBySequence(i);
    if(cut==NULL) continue;
    r=max(r,cut->maxX());
  }
  return r;
}

double XeSetOfCuts::maxY(){
  double r=VERY_SMALL;
  int n=getNCuts();
  for(int i=0;i<n;i++) {
    XeCut* cut=getCutBySequence(i);
    if(cut==NULL) continue;
    r=max(r,cut->maxY());
  }
  return r;
}


string XeSetOfCuts::getCutNameBySequence(int i){
  return checkIndex(i)? getCutBySequence(i)->getName():UNDEFINED_STRING;
}

void XeSetOfCuts::addCut(XeCut* cut) {
  int rn=cut->getRunNumber();
  if(run==NULL) {
     run=cut->getRun();
     runNumber=rn;
     run->attach(this);
  }
  else if(rn!=runNumber){
    cout<<"Trying to add cut of run "<<rn<<" to set of run "<<runNumber<<endl;
    return;
  }
  push_back(cut);
  cut->detachItself();
  attach(cut);  
  markAsChanged();
  setTheName();
}

void XeSetOfCuts::addCuts(XeSetOfCuts* cuts) {
  int n=cuts->getNCuts();
  for(int i=0;i<n;i++) addCut(cuts->getCutBySequence(i));
  cuts->detachItself();
}

void XeSetOfCuts::reset(){
  clear();
  setTheName();
}

void XeSetOfCuts::setTheName() {setName("Set of "+formatI(getNCuts())+" cuts");}

double XeSetOfCuts::getOverallAcceptance(){
  int n=getNCuts();
  double a=1.;
  for(int i=0;i<n;i++) a *= getCutBySequence(i)->getOverallAcceptance();
  return a;
}

bool XeSetOfCuts::printIt(int level){
  if(level<1) return true;
  cout<<getName()<<endl;
  int n=getNCuts();
  if(level>1) {
    for(int i=0;i<n;i++) {
      XeCut* cut=getCutBySequence(i);
      double overall=cut->getOverallAcceptance();
      cout<<"                  - "<<cut->getName();
      if(overall!=1.) cout<<" (accept. "<<overall<<")";
      cout<<endl;
    }
  }
  return true;
}

bool XeSetOfCuts::isEnabled(int i){
  return checkIndex(i)? getCutBySequence(i)->isEnabled():false;
}

void XeSetOfCuts::enable(int i,bool enabled){
  if(checkIndex(i)) getCutBySequence(i)->enable(enabled);
}
void XeSetOfCuts::disable(int i)          {enable(i,false);}
void XeSetOfCuts::disableAll()            {enableAll(false);}
void XeSetOfCuts::enableAll(bool enabled) {enableFirst(enabled?getNCuts():0); }
void XeSetOfCuts::enableFirst(int n) {
  int s=getNCuts();
  for(int i=0;i<s;i++) enable(i,i<n);
}


bool XeSetOfCuts::passes(double s1,double s2){
  int nc=getNCuts();
  for(int c=0;c<nc;c++) {
    if(!getCutBySequence(c)->passes(s1,s2)) {
      if(debugLevel>0){
        cout<<"s1="<<s1<<", s2="<<s2<<" fails cut "
            <<getCutBySequence(c)->getName()<<endl;
      }
      return false;
    }
  }
  return true;
}

XeMultiGraph* XeSetOfCuts::newMultiGraph(int plot){
  XeMultiGraph *mg=new XeMultiGraph("cuts of "+getName());
  for(int i=0;i<getNCuts();i++) mg->add(getCutBySequence(i)->newGraph());
  mg->drawIt(plot);
  return mg;
}


void XeSetOfCuts::draw() {
  if(S1S2Display::getMode()<0) return;
  if(debugLevel>1) cout<<"Drawing "<<getName()<<endl;
  for(int i=0;i<getNCuts();i++) {
    XeCut* cut=getCutBySequence(i);
    if(debugLevel>1) cout<<"   ...drawing "<<cut->getName()<<endl;
    cut->draw();
  }
}

XeSetOfCuts* XeSetOfCuts::newDefaultXeSetOfDarkMatterCuts(XeRun* r){
  XeSetOfCuts* sof=new XeSetOfCuts();
  if(XeStat::isCutsBased()) sof->addCut(new S1OverS2Cut(r,CENTRAL));
  return sof;
}

// ----------------------- SelectionCut ---------------------------------

int SelectionCut::tValueMode=LINEAR;

SelectionCut::SelectionCut(XeRun* r,int mode, bool all) : XeCut(r,"") {
  variable=all;
  run=r;
  tabulated=false;
  tValue=0.;
  smearMode=mode;
  nature=USER_SELECTION_CUT;
}

SelectionCut::SelectionCut() : XeCut(){}
SelectionCut::~SelectionCut(){}

void  SelectionCut::setTValueMode(int m) {
  switch(m) {
   case LINEAR :
   case EXPONENTIAL :
        tValueMode=m;
        return;
  }
  cout<<"Invalid SelectionCut t-value handling mode"<<endl;
}


void SelectionCut::setTValue(double t) {
  if(tValue==t) return;
  tValue=t;
  markAsChanged();
}

int       SelectionCut::getSmearMode()             {return smearMode;}
bool      SelectionCut::isVariable()               {return variable;}
double    SelectionCut::getTValue()                {return tValue;}

string SelectionCut::getSmearModeName(){
  switch(getSmearMode()){
    case SELECTION_CUT_UNKNOWN         : return "n/a";
    case SELECTION_CUT_ON_UNSMEARED_S1 : return "uS1";
    case SELECTION_CUT_ON_SMEARED_S1   : return " S1";
  }
  return UNDEFINED_STRING;
}

void SelectionCut::tabulate(){
  if(tabulated) return;
  int mode1=variable? ONE_SIGMA_BELOW:CENTRAL;
  int mode2=variable? ONE_SIGMA_ABOVE:CENTRAL;
  for(int i=0;i<N_PE_POINTS;i++){
    double x=PE_STEP*i;
    if(smearMode==SELECTION_CUT_ON_UNSMEARED_S1 && x<run->getSelectionUS1Min()){
      for(int m=mode1;m<=mode2;m++) _acc[i][m]=0.;
    }
    else {
      for(int m=mode1;m<=mode2;m++) _acc[i][m]=computeAcceptance(x,m);
    }
  }
  tabulated=true;
  markAsChanged();
}

bool SelectionCut::printIt(int level){
  if(level<1) return true;
  tabulate();
  cout<<endl<<getName()<<endl;
  for(int i=0;i<N_PE_POINTS;i++){
    double x=PE_STEP*i;
    cout<<setw(6)<<x<<" "<<_acc[i][0]<<" "<<_acc[i][1]<<" "<<_acc[i][2]<<endl; 
  }
  cout<<endl;
  return true;
}


double SelectionCut::getAcceptance(double S1){return getAcceptance(S1,tValue);}

double SelectionCut::getAcceptance(double S1,double t){
  if(!enabled) return 1.;
  double acceptance=getAcceptance(S1,CENTRAL);
  if(variable){
    if(t>0.){
      acceptance=XeInterpolation::interpolate(tValueMode,t,acceptance
                                 ,getAcceptance(S1,ONE_SIGMA_ABOVE));
    }
    else if(t<0.){
      acceptance=XeInterpolation::interpolate(tValueMode,-t,acceptance 
                                 ,getAcceptance(S1,ONE_SIGMA_BELOW));
    }
  }
  else if(t!=0.) {
    cout<<"Can't get acceptance of "<<getName()<<" for t!=0"<<endl;
    return 0.;
  }
  return min(1.,max(0.,acceptance));
}

double SelectionCut::getAcceptance(double S1,int mode){
  tabulate();
  double dex=S1/PE_STEP;
  int index=(int)dex;
  if(index<1) return 0.;
  else if(index>=N_PE_POINTS-1) return _acc[N_PE_POINTS-1][mode];
  double x=dex-index;
  return x*_acc[index+1][mode]+(1.-x)*_acc[index][mode];
}

XeGraph* SelectionCut::newGraphOfCurrentAcceptance(S1Range* pe,int plot){
  return newGraphOfAcceptance(pe,tValue,plot);
}

XeGraph* SelectionCut::newGraphOfAcceptance(S1Range* s1,double t,int plot){
  if(s1==NULL) s1=XeRange::getDefaultS1Range();
  int n=s1->getNPoints();
  XeGraph* g=new XeGraph("Acceptance of "+getName(),n,LABEL_S1);
  for(int i=0;i<n;i++){
    double p=s1->getValue(i);
    double a=getAcceptance(p,t);
    g->setPoint(i,p,a);
  }
  attach(g);
  g->drawIt(plot);
  return g;
}


SelectionCutS1::~SelectionCutS1(){}
SelectionCutS1::SelectionCutS1() : SelectionCut() {}
SelectionCutS1::SelectionCutS1(XeRun* r,bool var)
 : SelectionCut(r,SELECTION_CUT_ON_SMEARED_S1,var){
  nature=USER_SELECTION_CUT_ON_S1;
}

SelectionCutUnsmearedS1::~SelectionCutUnsmearedS1(){}
SelectionCutUnsmearedS1::SelectionCutUnsmearedS1(): SelectionCut(){}
SelectionCutUnsmearedS1::SelectionCutUnsmearedS1(XeRun* r,bool var)
 : SelectionCut(r,SELECTION_CUT_ON_UNSMEARED_S1,var){
  nature=USER_SELECTION_CUT_ON_UNSMEARED_S1;
}

bool s1sTotCut::isRunCompatible(int r){return r==8 ||r==10;}
double s1sTotCut::defaultS1Min(int run){
  switch(run){
    case RUN_08 : return DEFAULT_MIN_S1_RUN_08;
    case RUN_10 : return DEFAULT_MIN_S1_RUN_10;
    case RUN_12 : return DEFAULT_MIN_S1_RUN_10;
  }
  cout<<"Xs1sTot isn't defined for run "<<run<<endl;
  return UNDEFINED;
}

double s1sTotCut::defaultS1Max(int run){
  switch(run){
    case RUN_08 : return DEFAULT_MAX_S1_RUN_08;
    case RUN_10 : return DEFAULT_MAX_S1_RUN_10;
    case RUN_12 : return DEFAULT_MAX_S1_RUN_10;
  }
  cout<<"Xs1sTot isn't defined for run "<<run<<endl;
  return UNDEFINED;
}

s1sTotCut::s1sTotCut()         : SelectionCutS1(){}
s1sTotCut::s1sTotCut(XeRun* r) : SelectionCutS1(r,false){
  MinS1=defaultS1Min(runNumber);
  MaxS1=defaultS1Max(runNumber);
  cutMode=CUT_BELOW;
  nature= S1_STOT_SELECTION_CUT;
  setTheName();
  tabulate();
}

void s1sTotCut::setTheName(){
  char n[100];
  sprintf(n,"Xs1sTot   %.1f<S1<%.1f pe",MinS1,MaxS1);
  setName(n);
}

s1sTotCut::~s1sTotCut() {}

void   s1sTotCut::setMinS1(double minS1) {MinS1=minS1;tabulate();setTheName();}
void   s1sTotCut::setMaxS1(double maxS1) {MaxS1=maxS1;tabulate();setTheName();}
double s1sTotCut::getMinS1()             {return MinS1;}
double s1sTotCut::getMaxS1()             {return MaxS1;}

XeGraph* s1sTotCut::newGraph(int plot) {
  double x[4];
  double y[4];
  x[0]=MinS1      ; y[0]=TINY;
  x[1]=MinS1+.001 ; y[1]=VERY_LARGE;
  x[2]=MaxS1+.001 ; y[2]=VERY_LARGE;
  x[3]=MaxS1      ; y[3]=TINY;
  XeGraph* g=new XeGraph("Representation of "+getName(),4,x,y,LABEL_S1);
  attach(g);
  g->drawIt(plot);
  return g;
}

bool   s1sTotCut::passIt(double S1,double ){return S1>=MinS1 && S1<=MaxS1;}

double s1sTotCut::computeAcceptance(double S1, int mode){
  if(mode!=CENTRAL){
    cout<<"Can't compute "<<getName()<<" outside central value"<<endl;
    return 0.;
  }
  return passIt(S1,0.0)? 1.:0.;
}

void s2peaks0Cut::setTheName(){
  char title[132];
  sprintf(title,"Xs2peaks0     S2>%.1f pe",minS2);
  setName(title);
}

s2peaks0Cut::~s2peaks0Cut(){}
s2peaks0Cut::s2peaks0Cut()         : SelectionCutUnsmearedS1() {}
s2peaks0Cut::s2peaks0Cut(XeRun* r) : SelectionCutUnsmearedS1(r,false) {
  minS2=s2Min(runNumber);
  cutMode=CUT_ABOVE;
  nature=S2_PEAK_S0_SELECTION_CUT;
  setName("Xs2peaks0");
  switch(runNumber){
    case RUN_08 : setXenon100Reference(
                  "pl_spin_dependent_for_run_08_paper#acceptance");
                   break;
    case RUN_10 : setXenon100Reference("run10:s2gt150acc"); 
                   break;
    case RUN_12 : setXenon100Reference("run10:s2gt150acc"); 
                   break;
    default : cout<<"s2peaks0Cut can't handle run "<<runNumber<<endl;
  }
  setTheName();
  tabulate();
}
bool s2peaks0Cut::isRunCompatible(int r) {return r==8 || r==10;}
 
XeGraph* s2peaks0Cut::newGraph(int) { 
  /* do not mix up s2Peak with corrected S2!*/
  return NULL;
}

double s2peaks0Cut::computeAcceptance(double uS1, int mode) {
  if(mode!=CENTRAL){
    cout<<"Can't compute "<<getName()<<" outside central value"<<endl;
    return 0.;
  }
  switch(run->getNumber()) {
    case  RUN_08 :
       { int i=(int)(.001+uS1/PE_STEP);
         if(i<0 || i>=N_PE_POINTS) return 0.;
         return acc_run_08[i]; 
       }
    case RUN_10 :
       { double a= 1.00 - 0.08*exp(-0.33*(uS1-1.25));
         return max(a,0.);
       }
    case RUN_12 :
       { double a= 1.00 - 0.08*exp(-0.33*(uS1-1.25));
         return max(a,0.);
       }
  }
  cout<<getName()<<" can't handle run "<<run->getNumber()<<endl;
  return 0.;
}

double s2peaks0Cut::s2Min(int r){
  switch(r){
    case RUN_08 : return S2_MIN_RUN_8;
    case RUN_10 : return S2_MIN_RUN_10;
    case RUN_12 : return S2_MIN_RUN_10;
  }
  cout<<"Unknown run number for s2peaks0Cut :"<<r<<endl;
  return 0.;
}

S1coin2Cut::~S1coin2Cut(){}
S1coin2Cut::S1coin2Cut()         : SelectionCutS1()        {}
S1coin2Cut::S1coin2Cut(XeRun * r): SelectionCutS1(r,false) {
  nature=S1_COIN2_SELECTION_CUT;
  setName("Xs1coin2");
  setXenon100Reference(
    "cuts:xs1coin_veto#updateacceptance_for_final_dm_cuts_and_34kg_fv"
  );
  tabulate();
}
 
bool S1coin2Cut::isRunCompatible(int r) {return r==10;}
double S1coin2Cut::computeAcceptance(double uS1, int mode) {
  if(mode!=CENTRAL){
    cout<<"Can't compute "<<getName()<<" outside central value"<<endl;
    return 0.;
  }
  if(uS1<1.) return 0.;
  double a=1.- 0.427124*exp(-0.387355*(uS1-2.75385));
  return max(a,0.);
}


// --- this is for run 10, and exclude s1coin2
OtherSelectionCutsS1::~OtherSelectionCutsS1(){}
OtherSelectionCutsS1::OtherSelectionCutsS1()          : SelectionCutS1() {}
OtherSelectionCutsS1::OtherSelectionCutsS1(XeRun * r) : SelectionCutS1(r,true) {
  setName("OtherSelectionCutsS1");
  setXenon100Reference("run10ubp:acceptance");
  tabulate();
}
       
double OtherSelectionCutsS1::computeAcceptance(double uS1,int mode) {
  STATIC_CONST  double a0[N_SIGMA_MODES]={.771773,.777342,.782911};
  STATIC_CONST  double a1[N_SIGMA_MODES]={.00157241,.00183791,.00210341};

  double a=a0[mode]+a1[mode]*uS1;
  return min(1.,a);
}

//-- this is for run 8, and includes s1coin2

AllSelectionCutsS1::~AllSelectionCutsS1(){}
AllSelectionCutsS1::AllSelectionCutsS1() : SelectionCutS1(){}
AllSelectionCutsS1::AllSelectionCutsS1(XeRun * r): 
  SelectionCutS1(r,false) {
  setName("AllSelectionCutsS1");
  setXenon100Reference("pl_spin_dependent_for_run_08_paper#acceptance");
  tabulate();
}
bool   AllSelectionCutsS1::isRunCompatible(int r) {return r==8;}
double AllSelectionCutsS1::computeAcceptance(double uS1,int mode) {
  if(mode!=CENTRAL){
    cout<<"Can't compute "<<getName()<<" outside central value"<<endl;
    return 0.;
  }
  int i=(int)(.001+uS1/PE_STEP);
  if(i<0 || i>=N_PE_POINTS) return 0.;
  return acc_run_08[i];
}



// ------------------ Set of cuts -----------------------------------------

XeSetOfSelectionCuts::~XeSetOfSelectionCuts()  {reset();}
XeSetOfSelectionCuts::XeSetOfSelectionCuts() : XeSetOfCuts() {}

SelectionCut* XeSetOfSelectionCuts::getSelectionCutBySequence(int i){
  return (SelectionCut*)getCutBySequence(i);
}

SelectionCut* XeSetOfSelectionCuts::getSelectionCutByNature(int nature){
  return (SelectionCut*)getCutByNature(nature);
}

void XeSetOfSelectionCuts::setTValue(int i, double t){
  if(checkIndex(i)) (getSelectionCutBySequence(i))->setTValue(t);
}

void XeSetOfSelectionCuts::addSelectionCut(SelectionCut* cut) {
  addCut(cut);
}

void XeSetOfSelectionCuts::addSelectionCuts(XeSetOfSelectionCuts* cuts) {
  int n=getNCuts();
  for(int i=0;i<n;i++) addSelectionCut(cuts->getSelectionCutBySequence(i));
  setTheName();
}

void XeSetOfSelectionCuts::setTheName(){
  setName("Set of "+formatI(getNCuts())+" selection cuts");
}

double XeSetOfSelectionCuts::getAcceptance(double s1, bool smeared){
  int n=getNCuts();
  double a=1.;
  int want= smeared? SELECTION_CUT_ON_SMEARED_S1: SELECTION_CUT_ON_UNSMEARED_S1;  
  for(int i=0;i<n;i++) {
    SelectionCut *cut=getSelectionCutBySequence(i);
    if(cut->getSmearMode()==want) {
      a *= cut->getAcceptance(s1);
    }
  }
  return a;
}

double XeSetOfSelectionCuts::getMinS1(){
  SelectionCut* cut=getSelectionCutByNature(S1_STOT_SELECTION_CUT);
  if(cut!=NULL) return ((s1sTotCut*)cut)->getMinS1();
  cout<<"No Xs1sTot cut in current set of cuts!"<<endl;
  return UNDEFINED;
}

double XeSetOfSelectionCuts::getMaxS1(){
  SelectionCut* cut=getSelectionCutByNature(S1_STOT_SELECTION_CUT);
  if(cut!=NULL) return ((s1sTotCut*)cut)->getMaxS1();
  cout<<"No Xs1sTot cut in current set of cuts!"<<endl;
  return UNDEFINED;
}

bool XeSetOfSelectionCuts::printIt(int level) {return printIt(level,NULL);}
bool XeSetOfSelectionCuts::printIt(int level, vector<double>* remaining){
  if(level<1) return true;
  int n=getNCuts();
  cout<<getName()<<" composed "<<n<<" cuts"<<endl;
  printReference();
  
  if(remaining!=NULL) {
     cout<<"                          "
         <<leftJustify("Before all cuts",30)
         <<"    "<<setprecision(4)<<(*remaining)[0]<<RESET_PRECISION
         <<" events"<<endl;;
  }
  if(level<2) return true;
  for(int i=0;i<n;i++) {
     SelectionCut *c=getSelectionCutBySequence(i);
     cout<<"                  - ("<<c->getSmearModeName()
         <<") "<<leftJustify(c->getName(),30);
     if(remaining!=NULL) {
       cout<<"    "<<setprecision(4)<<(*remaining)[i+1]<<RESET_PRECISION
           <<" events";
     }
     cout<<endl;
     c->printReference(11);
     if(level>3) c->printIt();
  }
  return true;
}

XeSetOfSelectionCuts* XeSetOfSelectionCuts::newDefault(XeRun *r){
  XeSetOfSelectionCuts *ssc=new XeSetOfSelectionCuts();
  ssc->addCut(new s1sTotCut(r));
  ssc->addCut(new s2peaks0Cut(r));
  int rn=r->getNumber();
  switch(rn){
    case RUN_08 : ssc->addCut(new AllSelectionCutsS1(r));
                   break;
    case RUN_10 : ssc->addCut(new S1coin2Cut(r));
                  ssc->addCut(new OtherSelectionCutsS1(r));
                  break;
    case RUN_12 : ssc->addCut(new S1coin2Cut(r));
                  ssc->addCut(new OtherSelectionCutsS1(r));
                  break;
    default : cout<<" No default set of selection cuts for run "
                  <<r->getName()
                  <<endl;
              delete ssc;
              return NULL;
  }
  return ssc;
}


//------ cuts on S1/S2 ----------------------------------------------------

DarkMatterCut::~DarkMatterCut(){}
DarkMatterCut::DarkMatterCut() : XeCut() {} 
DarkMatterCut::DarkMatterCut(XeRun *r,string name,double eff) 
   : XeCut(r,getTheName(name,eff)){
  efficiency=eff;
  nature=DARK_MATTER_CUT;
}

DarkMatterCut::DarkMatterCut(XeRun* r, string name,int mode) : 
  XeCut(r,getTheName(name,mode)){
  efficiency=Rejection::getRejection(mode);
  nature=DARK_MATTER_CUT;
}

string DarkMatterCut::getTheName(string name, double eff){
  return name+" "+Rejection::getModeName(eff);
}

string DarkMatterCut::getTheName(string name, int mode){
  return name+" "+Rejection::getModeName(mode);
}

//---  Single cuts on S2/S1, basically S2<S2Max(s1) --------------------------

S1S2SingleCut::S1S2SingleCut() : DarkMatterCut() {} 
S1S2SingleCut::S1S2SingleCut(XeRun* r,string name,bool a ,double eff) 
  : DarkMatterCut(r,name,eff) {setMode(a);}
S1S2SingleCut::S1S2SingleCut(XeRun* r,string name,bool a, int mode)   
  : DarkMatterCut(r,name,mode) {
  setMode(a);
}

S1S2SingleCut::~S1S2SingleCut(){}

bool S1S2SingleCut::passIt(double S1,double S2){
  return above? S2>s2Cut(S1) : S2<s2Cut(S1);
}
void S1S2SingleCut::setMode(bool a){
  above=a;
  cutMode=above? CUT_ABOVE:CUT_BELOW;
  flattener=run->getFlattener();
}

XeGraph* S1S2SingleCut::newGraph(int plot){
  int p=S1S2Display::getMode();
  if(p<0) return NULL;
  S1Range pe(N_PE,PE_MIN,PE_MAX);
  XeGraph* g=new XeGraph("Representation of "+getName(),N_PE,LABEL_S1);
  for(int i=0;i<N_PE;i++){
    double s1=max(.0001,pe.getValue(i));
    double s2=s2Cut(s1);
    if(spaceMode==FLATTENED_S2_VS_S1) s2=flattener->unflatten(s1,s2);
    g->setPoint(i,s1,s2); 
  }
  attach(g);
  g->drawIt(plot);
  return g;
}

S1OverS2Cut::~S1OverS2Cut() { deleteWithPointer(cutFunction);}
S1OverS2Cut::S1OverS2Cut()  : S1S2SingleCut() {} 
S1OverS2Cut::S1OverS2Cut(XeRun* r,int mode) 
  : S1S2SingleCut(r,getOfficialName(r->getNumber()),false,mode) {
  setRunAndMode(r,mode);
}
  
S1OverS2Cut::S1OverS2Cut(XeRun* r,double e) 
  : S1S2SingleCut(r,getOfficialName(r->getNumber()),false,e){
  int mode=Rejection::getMode(e);
  setRunAndMode(r,mode);
}

void S1OverS2Cut::setRunAndMode(XeRun* r, int mode){
  run=r;
  setCutFunction(mode);
  overallAcceptance=getTheOverallAcceptance();
}

double S1OverS2Cut::getTheOverallAcceptance(){
   switch(runNumber){
    case RUN_10: 
      switch(cutMode) {
        case CENTRAL : return 10633./24527.; 
        default      : cout<<Name<<"  defined for CENTRAL mode only";
      }
      break;
    default: cout<<"Does not know the acceptance of "<<Name;
  }
  cout <<", assume 1."<<endl;
  return 1.;
}

string S1OverS2Cut::getOfficialName(int run){
  char title[80];
  sprintf(title," S1/S2 for run %d",run);
  return string(title);
}

void S1OverS2Cut::setCutFunction(int mode){
  
  STATIC_CONST  double run8Params[3][10]={
         {2.576/3.09*(-0.438131)+2.93418854 ,2.576/3.09*0.000597874-0.251661
          ,3.922036e-02 ,-3.777535e-03 ,2.2598315e-04 ,-8.5156134e-06
          ,2.01609e-07 ,-2.8995e-09 ,2.30693e-11 ,-7.76065e-14}
         ,{ 2.807/3.09*(-0.438131)+2.93418854,2.807/3.09*0.000597874-0.251661
          ,3.922036e-02,-3.777535e-03,2.2598315e-04,-8.5156134e-06
          ,2.01609e-07,-2.8995e-09,2.30693e-11,-7.76065e-14}
         ,{ -0.438131+2.93418854,0.000597874-0.251661,3.922036e-02
          ,-3.777535e-03,2.2598315e-04,-8.5156134e-06,2.01609e-07
          ,-2.8995e-09,2.30693e-11,-7.76065e-14}
         };

  cutFunction=NULL;
  rejectionMode=mode;
  int m=-1;
  switch(getRun()->getNumber()){
     case RUN_08: 
       switch(rejectionMode) {
          case REJECT995:  m=0; break;
          case REJECT9975: m=1; break;
          case REJECT999 : m=2; break;
          default:         cout<<"Can't set "<<Name<<" for rejection mode "
                               <<Rejection::getModeName(rejectionMode)<<endl;
                               return;     
       } 
       cutFunction = new TF1("cutFunction","pol9",2,50);
       cutFunction->SetParameters(run8Params[m][0],run8Params[m][1],run8Params[m][2]
                            ,run8Params[m][3],run8Params[m][4],run8Params[m][5]
                            ,run8Params[m][6],run8Params[m][7],run8Params[m][8]
                            ,run8Params[m][9]
                            );
       spaceMode=S2_VS_S1;
       break;
     case RUN_10:
       spaceMode=FLATTENED_S2_VS_S1;
       cutFunction= new TF1("cutFunction", "pol6", 3, 32);
       cutFunction->SetParameters(-0.252889, -0.0699762, 0.0138153, -0.00131545
                             , 6.57013e-05, -1.6588e-06, 1.66591e-08);
       break;
     default: 
       cout<<"Official cuts n.y.i for run "<<run<<endl;
       spaceMode=UNDEFINED_INT;       
       return;
  }
}

bool S1OverS2Cut::isRunCompatible(int r)       {return r==8 || r==10;}
bool S1OverS2Cut::passIt(double S1, double S2) {
  double s2=spaceMode==FLATTENED_S2_VS_S1 ? flattener->flatten(S1,S2) : S2;
  return s2<=s2Cut(S1);
}
 
double S1OverS2Cut::s2Cut(double S1)     {
  if(cutFunction!=NULL) {
    switch(spaceMode){
      case S2_VS_S1           : return S1*pow(10.,cutFunction->Eval(S1));
      case FLATTENED_S2_VS_S1 : return cutFunction->Eval(S1);
    }
    cout<<"Major inconsistecy, 'spaceMode' is invalid for "<<getName()<<endl;
  }
  else {
    cout<<"Major inconsistecy, 'cutFunction' isn't defined for "<<getName()<<endl;
  }
  return UNDEFINED;
}

//--------------  Run Flattener ------------------ 

RunFlattener::~RunFlattener(){ deleteWithPointer(eBandFlat);}
RunFlattener::RunFlattener(): Flattener(), RunComponent(){}

RunFlattener::RunFlattener(XeRun* r): Flattener(), RunComponent() {
  setRun(r);
  setName(getTheName(run));
  eBandFlat=NULL;
  eBandMin=.5;
  eBandMax=50.;
  switch(runNumber){
    case RUN_08 :
      eBandFlat = new TF1("ebandflat","pol9",eBandMin,eBandMax);
      eBandFlat->SetParameters( 2.66695+0.267047+0.00019154
                              ,-0.0960315-0.155731+0.000101306
                              , 0.00604496+0.0331754
                              ,-0.000200455-0.00357708
                              , 3.29415e-06+0.000222689,-2.10034e-08-8.49461e-06
                              , 2.01609e-07
                              ,-2.8995e-09
                              , 2.30693e-11
                              ,-7.76065e-14
                              ) ;
      break;

    case RUN_10 :
      eBandFlat = new TF1("ebandflat"
                ,"exp(-0.513894-0.274465*x)+2.25458+0.0131485+(-0.0148529+-0.0166898)*x+ (0.000198616+0.00460857)*x**2-0.000523096*x**3+2.85031e-05*x**4-7.41271e-07*x**5+7.38937e-09*x**6"
                ,eBandMin,eBandMax);
      break;

    case RUN_12 :
      eBandFlat = new TF1("ebandflat"
                ,"exp(-0.513894-0.274465*x)+2.25458+0.0131485+(-0.0148529+-0.0166898)*x+ (0.000198616+0.00460857)*x**2-0.000523096*x**3+2.85031e-05*x**4-7.41271e-07*x**5+7.38937e-09*x**6"
                ,eBandMin,eBandMax);
      break;

    default:
      cout<<"Don't know how to flatten bands for run "<<runNumber<<endl;
  }
  checkRunCompatibility();
}

bool RunFlattener::isRunCompatible(int r)   {return r==8 || r==10;}
string RunFlattener::getTheName(XeRun* run) {
  return "Flattener for "+run->getName();
}
 
double RunFlattener::flatten(double S1, double S2){
  if(eBandFlat==NULL){
    cout<<"Don't know how to flatten bands for this run "<<endl;
    return UNDEFINED;
  }
  else if(S1<eBandMin || S1 >eBandMax){
    cout<<"Don't know how to flatten bands for S1="<<S1<<", should be between "
        <<eBandMin<<" and "<<eBandMax<<endl;
    return UNDEFINED;
  }
  return S1<=0.? UNDEFINED : log10(S2/S1)-eBandFlat->Eval(S1);
}

double RunFlattener::unflatten(double S1, double flat){
  if(eBandFlat==NULL){
    cout<<"Don't know how to flatten bands for this run "<<endl;
    return UNDEFINED;
  }
  else if(S1<eBandMin || S1 >eBandMax){
    cout<<"Don't know how to flatten bands for S1="<<S1<<", should be between "
        <<eBandMin<<" and "<<eBandMax<<endl;
    return UNDEFINED;
  }
  return S1<=0.? UNDEFINED : S1* pow(10.,flat+eBandFlat->Eval(S1));
}


//-------------------------------------------------------- Set of S1 and S2


S1S2Data::~S1S2Data(){}
S1S2Data::S1S2Data() : S1S2Object() , DataSet(), RunComponent() {}

S1S2Data::S1S2Data(string name, int t, XeRun* r, TGraph* gr) : 
  S1S2Object() , DataSet(2), RunComponent(){
  setName(name);
  setRun(r);
  setParameters(t);
  if(gr!=NULL) add(gr);
}

S1S2Data::S1S2Data(S1S2Data* orig, DarkMatterCut* cut) : 
   S1S2Object() ,  DataSet(2) , RunComponent() {
  setName(getTheName(orig,cut));
  setRun(orig->getRun());
  setParameters(orig->getDataType());
  add(orig,cut);
  orig->setCutData(this);
}

S1S2Data::S1S2Data(S1S2Data* orig, XeSetOfCuts* cuts) :
  S1S2Object() ,  DataSet(2) , RunComponent() {
  setName(getTheName(orig,cuts));
  setRun(orig->getRun());
  setParameters(orig->getDataType());
  add(orig,cuts);
  orig->setCutData(this);
}

void S1S2Data::setParameters(int t){
  dataType=t;
  normalization=1.;
  cutData=NULL;
  s1Min=UNDEFINED;
  s1Max=UNDEFINED;
  s2Min=UNDEFINED;
  s2Max=UNDEFINED;
  s2overs1Min=UNDEFINED;
  s2overs1Max=UNDEFINED;
}

void S1S2Data::printSummaryHeader(){
  string line(80,'-');
  cout<<endl<<line<<endl
      <<"Dataset        N.events       S1         S2/S1    N.events.after.normalization"<<endl
      <<"              All   Cut    Min   Max   Min   Max      Factor   All      Cut"<<endl
      <<line<<endl;
}

void S1S2Data::printSummary(string header){
  trace(ENTRY,"S1S2Data::printSummary()");
  string ne=string(6,' ');
  string no=string(12,' ');
  string nn=string(8,' ');
  if(cutData!=NULL){
    ne=formatI(cutData->getNEvents(),6);
    no=formatG(cutData->getNormalization(),12,3);
    nn=formatF(cutData->getNEventsNormalized(),8,2);
  } 
  cout<<header
      <<setw(6)<<getNEvents()<<ne
      <<formatF(minS1(),7,2)  
      <<formatF(maxS1(),6,2)
      <<formatF(minS2OverS1(),6,1)
      <<formatF(maxS2OverS1(),7,1)
      <<no
      <<formatF(getNEventsNormalized(),9,2)
      <<nn
      <<RESET_PRECISION<<endl;
  trace(EXIT,"S1S2Data::printSummary()");
}


bool S1S2Data::printIt(int level) {
  int n=getNEvents();
  bool ok=n>0;
  if(level<1) return ok;
  cout<<getName()<<" has "<<n<<" entries"<<endl;
  if(level<=1) return ok;
  for(int e=0;e<n;e++) {
    double s1=getS1(e);
    double s2=getS2(e);
    double f=run->flatten(s1,s2);
    printf(" %3d S1=%4.1f  S2=%6.1f  Flattened:%6.3f\n",e,s1,s2,f);
  }
  return ok;
}

void S1S2Data::printComputedBands(int level,int bandMax){
  if(level<=0) return;
  if(analysisMode!=PL_ANALYSIS) { 
    printIt(level);
    return;
  }
  int n=getNEvents();
  if(bandMax==ALL) cout<<getName()<<" has "<<n<<" entries"<<endl;
  else cout<<"Events with band<="<<bandMax<<" in "<<getName()<<endl;;
  if(level<1) return ;
  for(int e=0;e<n;e++) {
    if(bandMax==ALL || bands[e]> bandMax) continue;
    double s1=getS1(e);
    double s2=getS2(e);
    double f=run->flatten(s1,s2);
    printf(" %3d S1=%4.1f  S2=%6.1f  Band=%2d   Slice=%3d  Flattened:%6.3f\n"
          ,e,s1,s2,bands[e],slices[e],f);
  }
}

S1S2Data* S1S2Data::getCutData() {return cutData;}
 
void S1S2Data::setCutData(S1S2Data* c) {cutData=c;}
void S1S2Data::reset()                 {clear();slices.clear();bands.clear();}
void S1S2Data::normalize(double ratio) {normalization=ratio;update();}
void S1S2Data::normalizeToEvents(double events){
  nEvents=getNEvents();
  if(nEvents==0){
    cout<<"Can't normalize empty "<<getName()<<endl; 
    return;
  }
  normalize(events/nEvents);
}

XeGraph* S1S2Data::newGraph(int plot){
  if(S1S2Display::getMode()<0) return NULL;
  int n=getNEvents();
  XeGraph *gr=new XeGraph(getName(),n,getS1(),getS2(),LABEL_S1,LABEL_S2);
  gr->setMarkerColor(defaultColor());  
  attach(gr);
  gr->drawIt(plot);
  return gr;
}

void S1S2Data::draw(){
  XeGraph *g=newGraph();
  int n=g->getN();
  if(debugLevel>1){
    cout<<"Drawing "<<getName()<<" consisting of "<<n<<" events"<<endl;
    g->printIt();
  }
  if(n>0) g->drawS1S2("P");
}

void S1S2Data::make(TGraph* gr) {clear();add(gr);}
void S1S2Data::add(TGraph* gr) {
  int n=gr->GetN();
  if(n<=0 ||n >MAX_NUMBER_OF_EVENTS){
    cout<<"Too many events in "<<getName()<<" : "<<n<<endl;
    n=0;
    markError();
    return;
  }
  double* s1=gr->GetX();
  double* l=gr->GetY();
  for(int i=0;i<n;i++) {
    double s2=s1[i]*pow(10.,l[i]);
    add(s1[i],s2);
  }
  printIt(printLevel-1);
  update();
}

string S1S2Data::getTheName(S1S2Data* original, XeSetOfCuts* cuts){
  return original->getName()+" after "+cuts->getName();
}

string S1S2Data::getTheName(S1S2Data* original, DarkMatterCut* cut){
  return original->getName()+" after "+cut->getName();
}

void S1S2Data::make(S1S2Data* orig, XeSetOfCuts* cuts)  {clear();add(orig,cuts);}
void S1S2Data::add(S1S2Data* original, XeSetOfCuts* cuts){
  int n=original->getNEvents();
  XeRun* r=original->getRun();
  if(n<=0 ||n >MAX_NUMBER_OF_EVENTS){
    cout<<"Too many events in "<<getName()<<" : "<<n<<endl;
    n=0;
    markError();
    return;
  }

  if(run==NULL) run=r;
  else if(run!=r){
    cout<<"Mismatch between runs : "<<r->getName()
                                 <<" vs "<<run->getName()<<endl;
    n=0;
    markError();
    return;
  }

  for(int i=0;i<n;i++){
    double s1=original->getS1(i);
    double s2=original->getS2(i);
    if(cuts->passes(s1,s2)) add(s1,s2);
    else {
      if(printLevel>2) {
        cout<<"Event failing the cut in S1S2Data  " << original->getName() 
             << " : s1="<<s1<<", s2="<<s2<<endl;
      }
    }
  }
  update();
  if(printLevel>0) {
    cout<<getName()<<" has "<<getNEvents()<<" events"<<endl;
  }
}

EquiContentBins* S1S2Data::newEquiContentS1Bins(int nS1){
  vector<double>*  s1=getColumn(0);
  EquiContentBins* eq=new EquiContentBins("S1 in "+getName(),nS1,&s1[0]);
  return eq;
}

EquiContentBins* S1S2Data::newEquiContentS1Bins(int nS1,double s1Mi
                                                          , double s1Ma){
  if(debugLevel>0){
    cout<<"Getting"<<nS1<<" equicontent bins of "<<getName()
        <<" which has a total of "<<getNEvents()<<" events"<<endl;
  }
  vector<double>*  s1=getColumn(0);
  EquiContentBins* eq=new EquiContentBins("S1 in "+getName()
                                        ,nS1,&s1[0],s1Mi,s1Ma);
  return eq;
}

void S1S2Data::make(S1S2Data* original, DarkMatterCut* cut) {
  clear();
  add(original,cut);
}

void S1S2Data::add(S1S2Data* original, DarkMatterCut* cut){
  int n=original->getNEvents();
  for(int i=0;i<n;i++){
   double s1=original->getS1(i);
   double s2=original->getS2(i);
   if(cut->passes(s1,s2)) add(s1,s2);
  }
}

bool S1S2Data::isAfterCuts(int t) {
  return t>=N_DATA_TYPES && t<N_CUT_DATA_TYPES;
}

int S1S2Data::typeAfterCut(int type){
  switch(type){
    case DM_DATA          : 
    case DM_CUT_DATA      : return DM_CUT_DATA;
    case AM_BE_DATA       : 
    case AM_BE_CUT_DATA   : return AM_BE_CUT_DATA;
    case E_GAMMA_DATA     : 
    case E_GAMMA_CUT_DATA : return E_GAMMA_CUT_DATA;
  }
  cout<<"Don't know the type after cut for type "<<type<<endl;
  return UNDEFINED_INT;
}

int S1S2Data::typeBeforeCut(int type){
  switch(type){
    case DM_DATA          : 
    case DM_CUT_DATA      : return DM_DATA;
    case AM_BE_DATA       : 
    case AM_BE_CUT_DATA   : return AM_BE_DATA;
    case E_GAMMA_DATA     : 
    case E_GAMMA_CUT_DATA : return E_GAMMA_DATA;
  }
  cout<<"Don't know the type before cut for type "<<type<<endl;
  return UNDEFINED_INT;
}

string S1S2Data::formatContent(int t, double c,double w, bool total){
  
  if(total) {
    switch(t) {  
      case DM_CUT_DATA_INTEGRATED           :  
      case DEFAULT_SIGNAL_INTEGRATED        :  
      case ALL_BACKGROUNDS_INTEGRATED  : return string(w,' ');         
    }
  }
  switch(t) {  
    case DM_DATA                :
    case DM_CUT_DATA            :
    case DM_CUT_DATA_INTEGRATED :
         return formatI((int)(c+.5),w);
  }
  return formatF(c,w,2);
}

bool   S1S2Data::isBackground(int type, bool warn){
  bool ok= type>=FIRST_BACKGROUND && type<=LAST_BACKGROUND;
  if(warn && !ok) cout<<"Warning ! type '"<<type<<"' isn't background"<<endl;
  return ok;
}

string S1S2Data::getShortTypeName(int m){
  switch(m){
    case UNSPECIFIED_DATA                : return "Unspecified";
    case DM_DATA                         : return "Dark matter";
    case AM_BE_DATA                      : return "Am-Be";
    case E_GAMMA_DATA                    : return "E Gamma";
    case AM_BE_CUT_DATA                  : return "Am-Be after cuts";
    case E_GAMMA_CUT_DATA                : return "E Gamma after cuts";
    case DM_CUT_DATA                     : return "Data";       
    case ER_BACKGROUND                   : return "Gauss.";
    case ER_LEAKAGE                      : return "Anomal.";
    case NR_BACKGROUND                   : return "Neutron";
    case ALL_BACKGROUNDS                 : return "All Bkg";   
    case DM_CUT_DATA_INTEGRATED          : return "Integr.";         
    case ALL_BACKGROUNDS_INTEGRATED      : return "Integr.";         
    case DEFAULT_SIGNAL                  : 
         return formatG(DEFAULT_SIGMA_NUCLEON,8,2);
    case DEFAULT_SIGNAL_INTEGRATED       : return "Integr.";         
  }
  return UNDEFINED_STRING;
}

string S1S2Data::getTypeName(int m){
  switch(m){
    case UNSPECIFIED_DATA                : return "Unspecified";
    case DM_DATA                         : return "Dark matter data";
    case AM_BE_DATA                      : return "Am-Be data";
    case E_GAMMA_DATA                    : return "E Gamma data";
    case AM_BE_CUT_DATA                  : return "Am-Be data after cuts";
    case E_GAMMA_CUT_DATA                : return "E Gamma data after cuts";
    case DM_CUT_DATA                     : return "DM data after cuts";
    case ER_BACKGROUND                   : return "ER gaussian";
    case ER_LEAKAGE                      : return "ER leakage";
    case NR_BACKGROUND                   : return "NR background";
    case ALL_BACKGROUNDS                 : return "All backgrounds";
    case DM_CUT_DATA_INTEGRATED          : return "(integrated)";         
    case ALL_BACKGROUNDS_INTEGRATED      : return "(integrated)";         
    case DEFAULT_SIGNAL                  : return "Expected signal";
    case DEFAULT_SIGNAL_INTEGRATED       : return "(integrated)";         
  }
  return UNDEFINED_STRING;
}

int S1S2Data::mapS1Shape(int s1s2Shape){
  switch (s1s2Shape) {
    case SHAPE_SIGNAL          : return DEFAULT_SIGNAL;
    case SHAPE_ER_BACKGROUND   : return ER_BACKGROUND;
    case SHAPE_ER_LEAKAGE      : return ER_LEAKAGE;
    case SHAPE_NR_BACKGROUND   : return NR_BACKGROUND;
    case SHAPE_DATA            : return DM_CUT_DATA;
  }
  return UNDEFINED_INT;
}

int S1S2Data::defaultColor(int m){
  switch(m){
    case AM_BE_DATA         : return kRed;
    case E_GAMMA_DATA       : return kBlue;
  }
  return kBlack;
}

int    S1S2Data::getDataType()  {return dataType;}
int    S1S2Data::defaultColor() {return defaultColor(dataType);}

bool   S1S2Data::update()       {

  nEvents=getNEvents();
  nEventsNormalized=nEvents*normalization;
  if(printLevel>1){
    cout<<Name<<" has "<<nEvents<<" events, normalized to "<<nEventsNormalized
        <<endl;
  }

  vector<double>*  s1=getColumn(S1_COLUMN);
  vector<double>*  s2=getColumn(S2_COLUMN);
  s1Max=VERY_SMALL;
  for(int i=0;i<nEvents;i++) s1Max=max(s1Max,(*s1)[i]);
  s2Max=VERY_SMALL;
  for(int i=0;i<nEvents;i++) s2Max=max(s2Max,(*s2)[i]);
  s2overs1Max=VERY_SMALL;
  for(int i=0;i<nEvents;i++) s2overs1Max=max(s2overs1Max,(*s2)[i]/(*s1)[i]);

  s1Min=VERY_LARGE;
  for(int i=0;i<nEvents;i++) s1Min=min(s1Min,(*s1)[i]);
  s2Min=VERY_LARGE;
  for(int i=0;i<nEvents;i++) s2Min=min(s2Min,(*s2)[i]);
  s2overs1Min=VERY_LARGE;
  for(int i=0;i<nEvents;i++) s2overs1Min=min(s2overs1Min,(*s2)[i]/(*s1)[i]);

  return true;
}

string S1S2Data::getTypeName() {return getTypeName(dataType);}
double S1S2Data::maxS1()       {return s1Max;}
double S1S2Data::maxS2()       {return s2Max;}
double S1S2Data::maxS2OverS1() {return s2overs1Max;}
double S1S2Data::minS1()       {return s1Min;}
double S1S2Data::minS2()       {return s2Min;}
double S1S2Data::minS2OverS1() {return s2overs1Min;}

double S1S2Data::minX(){
  int size=getNEvents();
  double m=VERY_LARGE;
  for(int i=0;i<size;i++) m=min(m,getX(i));
  return m;
}

double S1S2Data::minY(){
  int size=getNEvents();
  double m=VERY_LARGE;
  for(int i=0;i<size;i++) m=min(m,getY(i));
  return m;
}

double S1S2Data::maxX(){
  int size=getNEvents();
  double m=0.;
  for(int i=0;i<size;i++) m=max(m,getX(i));
  return m;
}

double S1S2Data::maxY(){
  int size=getNEvents();
  double m=0.;
  for(int i=0;i<size;i++) m=max(m,getY(i));
  return m;
}

TTree* S1S2Data::newTree(S1S2Bands* ref){
  bool withSB= ref!=NULL ;
  TTree *tree =new TTree(getNameChar(),getNameChar());
  tree->SetDirectory(NULL);
  tree->Branch("S1"             , &S1           , "S1/F");
  tree->Branch("S2"             , &S2           , "S2/F");
  tree->Branch("logS2S1"        , &logS2S1      , "logS2S1/F");
  tree->Branch("flattened"      , &flattened    , "flattened/F");
  if(withSB){
    tree->Branch("slice"        , &slice        , "slice/I");
    tree->Branch("band"         , &band         , "band/I");
  }
  for(int e=0;e<nEvents;e++) {
    S1=getS1(e);
    S2=getS2(e);
    logS2S1= S1==0.? UNDEFINED: log10(S2/S1);
    flattened=run->flatten(S1,S2);
    if(withSB){
      pair<int,int> sb=ref->whichSliceAndBand(S1,S2);
      slice=sb.first;      
      band=sb.second;
    }
    tree->Fill(); 
  }
  return tree;
}

void S1S2Data::add(double s1, int s, double s2, int b) {
  add(s1,s2);
  slices.push_back(s);
  bands.push_back(b);
}

void S1S2Data::add(double s1, double s2) {
  double p[2];
  p[S1_COLUMN]=s1;
  p[S2_COLUMN]=s2;
  addEntry(p);
}

bool   S1S2Data::passes(int i, XeSetOfCuts* cuts)  {
  double s1=getS1(i);
  double s2=getS2(i);
  return cuts->passes(s1,s2);
}

bool   S1S2Data::passes(int i,DarkMatterCut* cut) {
  return cut->passes(getS1(i),getS2(i));
}

int     S1S2Data::getNEvents()           {return nEvents;}
int     S1S2Data::getBand(int ev)        {return bands[ev];}
int     S1S2Data::getSlice(int ev)       {return slices[ev];}
double  S1S2Data::getNormalization()     {return normalization;}
double  S1S2Data::getNEventsNormalized() {return nEventsNormalized;}
double  S1S2Data::getS1(int ev)          {return getValue(ev,S1_COLUMN);}
double  S1S2Data::getS2(int ev)          {return getValue(ev,S2_COLUMN);}
double  S1S2Data::getX(int ev)           {return getS1(ev);}
double* S1S2Data::getS1()                {return getVector(S1_COLUMN);}
double* S1S2Data::getS2()                {return getVector(S2_COLUMN);}

double  S1S2Data::getY(int ev)  {return S1S2Display::getY(getS1(ev),getS2(ev));}
double  S1S2Data::getS2overS1(int ev)    {
  double s1=getS1(ev);
  return s1<=0.? UNDEFINED: getS2(ev)/s1;
}

void S1S2Data::fillAndMarkBandsAndSlices(S1S2Bands* s1s2b,bool mark){
  int ne=getNEvents();
  bands.resize(ne);
  slices.resize(ne);
  for(int e=0;e<ne;e++) {
    double s1=getS1(e); 
    double s2=getS2(e); 
    pair<int,int> sb=s1s2b->whichSliceAndBand(s1,s2);
    int s=sb.first; 
    int b=sb.second;
    if(b>=0 && s>=0) s1s2b->fill(s1,s,b);
    if(mark){
      bands[e]=b;
      slices[e]=s;
    }
  }
}

//--------------- a virtual class  for calculating pvalues ---------------

S1S2PValue::~S1S2PValue(){}

S1S2PValue::S1S2PValue(string n,int m,XeRun* r) : PValue(n,m), RunComponent(){
  setName(shortName);
  setRun(r);
  if(run!=NULL){
    Name += " for "+run->getName();
    if(debugLevel>0) {
      cout<<getName()<<" properly linked to "<<run->getName()<<endl;
    }
  }
}

void S1S2PValue::setInteraction(Interaction* inter){
  run->setInteraction(inter);
}

void   S1S2PValue::updateSigmaScale() {}

bool S1S2PValue::checkEverything(){
  if(run==NULL){
    cout<<"No S1S2 run for "<<getName()<<endl;
    return false;
  }
  run->checkTheCuts();
  if(run==NULL){
    cout<<"No run for "<<getName()<<endl;
    return false;
  }
  return run->checkIt();
}

void S1S2PValue::setWimpMass(double m) {
  if(printLevel>1) cout<<"S1S2PValue sets wimp mass to "<<m<<endl;
  run->setWimpMass(m);
}

double S1S2PValue::nSignalPerCm2(){
  return run->expectedSignal(0.)/run->getSigmaNucleon();
}

//---------------- Simple signal/bkg for CLs---------------------------------

S1S2CLs::~S1S2CLs(){}
S1S2CLs::S1S2CLs(XeRun* r) : S1S2PValue("CLs",CUTS_ANALYSIS,r){}

double S1S2CLs::pValueS(double sigma){
  MethodCounter::count("pValue for CLS");
  run->setSigmaNucleon(DEFAULT_SIGMA_NUCLEON);
  run->computeSignal();
  double signal=sigma*nSignalPerCm2();
  double mu=bkg+signal;
  double p=PoissonDist::belowEqual(nEvents,mu);
  if(printLevel>2){
    cout<<"S1S2CLs::pValueS(sigma="<<sigma<<"), signal="<<signal
        <<",  mu="<<mu<<", p="<<p<<endl;
  }  
  return p;
}

void S1S2CLs::printFlagsAndParameters() {cout<<endl<<getName()<<endl<<endl;}

bool S1S2CLs::update() {
  run->update();
  return checkPValue();
}

bool S1S2CLs::checkPValue(){
  if(!checkEverything()) return false;
  if(printLevel>2){
    cout<<" Initialization/Update of "<<getName()<<endl;
  }
  nEvents=(int)(run->getS1S2Data(DM_CUT_DATA)->getNEvents());
  bkg=run->getS1S2Data(E_GAMMA_CUT_DATA)->getNEventsNormalized()
     +run->getS1S2Data(AM_BE_CUT_DATA)->getNEventsNormalized();
  if(printLevel>0){
    cout<<"DM events after cuts  : "<<nEvents<<endl
        <<"Normalized background : "<<bkg<<endl;
  }
  return true;
}


//--------Slices and Bands --------------------------------------------

string S1Slice::getTheName(int s,S1S2Bands *s1s2b){
  char title[200];
  string t=s1s2b->getName();
  XeBins *s1Bins=s1s2b->getS1Bins();
  double s1min=s1Bins->getLowerEdge(s);
  double s1max=s1Bins->getUpperEdge(s);
  sprintf(title,"%s slice %2d, %5.2f<S1<%5.2f",t.c_str(),s,s1min,s1max);
  return title;
}

S1Slice::~S1Slice() {}
S1Slice::S1Slice() : DetectorComponent(), S1S2Object(), RunComponent() {} 
S1Slice::S1Slice(int s,S1S2Bands *s1s2b) : DetectorComponent()
                                        , S1S2Object(), RunComponent() {
  setName(getTheName(s,s1s2b));
  if(debugLevel>1){
    cout<<"Building "<<getName()<<endl;
  }
  mother=s1s2b;
  setRun(mother->getRun());
  XeBins *s1Bins=mother->getS1Bins();
  sequence=s;
  S1Min=s1Bins->getLowerEdge(sequence);
  S1Max=s1Bins->getUpperEdge(sequence);
  nBands=mother->getNBands();
  count=0;
  for(int b=0;b<nBands;b++) counts.push_back(0);

  S2overS1Bins = new XeBins(); 
}

void S1Slice::reset(){
  S2overS1Values.clear();
  count=0.;
  for(int b=0;b<nBands;b++) counts[b]=0.;
}

void S1Slice::normalize(double ratio){
  for(int b=0;b<nBands;b++) counts[b]*=ratio;
  count*=ratio;
}

void S1Slice::reweightTheBands(double * weights){
  count=0.; 
  for(int b=0;b<nBands;b++) {
    counts[b]*=weights[b];
    count+=counts[b];
  }
}

void   S1Slice::setS1Min(double s)   {S1Min=s;}
void   S1Slice::setS1Max(double s)   {S1Max=s;}
double S1Slice::getS1Min()           {return S1Min;}
double S1Slice::getS1Max()           {return S1Max;}
double S1Slice::getS1Width()         {return S1Max-S1Min;}
double S1Slice::getCentralS1()       {return .5*(S1Max-S1Min);}
double S1Slice::getContent()         {return count;}
double S1Slice::getContent(int b)    {return b>=0 && b<nBands? counts[b]:0.;}
double S1Slice::getS2overS1Min()     {return S2overS1Bins->getLowerEdge();}
double S1Slice::getS2overS1Max()     {return S2overS1Bins->getUpperEdge();}
double S1Slice::minX()               {return getS1Min();}
double S1Slice::maxX()               {return getS1Max();}
double S1Slice::minY()               {return getS2overS1Min();}
double S1Slice::maxY()               {return getS2overS1Max();}

double S1Slice::getArea() {
  return getS1Width()*(log10(getS2overS1Max()/getS2overS1Min()));
}

void   S1Slice::fill(double s1, double s2,double w) {
  if(s1<=0.) return;
  fill(S2overS1Bins->getBin(s2/s1),w);
}

void   S1Slice::fill(int bin,double w){
  if(bin>=0) {
    count+=w;
    counts[bin]+=w;
  }
}

double  S1Slice::getMaxCount(){
  double m=0.;
  for(int b=0;b<nBands;b++) m=max(m,counts[b]);
  return m;
}

double  S1Slice::getMinCount(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) m=min(m,counts[b]);
  if(m==VERY_LARGE) m=0.;
  return m;
}

bool S1Slice::printIt(int level){
  if(level<1) return true;
  cout<<getName()<<", count: "<<count<<endl;
  if(level>1){
    printVector(counts);
    if(level>2){
      printVector(S2overS1Bins->getLowerEdges(),nBands,"Min S2/S1");
      printVector(S2overS1Bins->getUpperEdges(),nBands,"Max S2/S1");
    }
    cout<<endl;
  }
  return true;
}

bool S1Slice::printItLog(int level){
  if(level<1) return true;
  cout<<getName()<<", count: "<<count<<endl;
  if(level>1){
    printVector(counts);
    if(level>2){
      printLogVector(S2overS1Bins->getLowerEdges(),nBands,"Min Log10(S2/S1)");
      printLogVector(S2overS1Bins->getUpperEdges(),nBands,"Max Log10(S2/S1)");
    }
    cout<<endl;
  }
  return true;
}

void S1Slice::draw(){
  if(S1S2Display::getMode()<0) return;
  double x[2];
  double y[2];
  y[0]=S1S2Display::getYmin();
  y[1]=S1S2Display::getYmax();
  for(int edge=0;edge<2;edge++){
    x[0]=x[1]=edge==0? S1Min:S1Max;
    XeGraph* g= new XeGraph("",2,x,y,LABEL_S1,LABEL_S2);
    g->setLineWidth(1);
    g->Draw("L");
  }
}


void S1Slice::fillForStudy(double s1,double s2) {
  if(s1>=0.) S2overS1Values.push_back(s2/s1);
}
void S1Slice::setS2overS1Bins(XeBins* s2s1)  {*S2overS1Bins= *s2s1;} //FIXED: SHOULD NOT initialize class member from a pointer to local variable 
vector<double>* S1Slice::getS2overS1Vector() {return &S2overS1Values;}
XeBins*         S1Slice::getS2overS1Bins()   {return S2overS1Bins;}


//--- The bands, which contain other info vital for analysis -----------------
S2Band::S2Band(): S1S2Object(), RunComponent() {}

S2Band::S2Band(int b,S1S2Bands *s1s2b): S1S2Object(), RunComponent() {
   mother=s1s2b;
   setName(getTheName(b,s1s2b));
   setRun(mother->getRun());
   s1Bins=mother->getS1Bins();
   sequence=b;
   nSlices=s1Bins->getNBins();
   count=0;

   // initialization of S1 bins relative of the band
   for(int i=0;i<nSlices;i++) {
	counts.push_back(0.);
	S2overS1Min.push_back(0.);
	S2overS1Max.push_back(0.);
   }
}

string S2Band::getTheName(int b,S1S2Bands *s1s2b){
  char title[200];
  string t=s1s2b->getName();
  sprintf(title,"%s band %2d",t.c_str(),b);
  return title;
}

void S2Band::fill(double s1,double w){
  fill(s1Bins->getBin(s1),w);
  if(w==1.) S1.push_back(s1);
}

vector<double>* S2Band::getS1() {return &S1;}

vector<int> S2Band::getS1InBin() {

	vector <int> corresponding_s1_bins;

	if (S1.size() == 0 ) return corresponding_s1_bins;
	
	for(unsigned int i =0; i < S1.size(); i++) {
		corresponding_s1_bins.push_back( s1Bins->getBin(S1[i]) );
	}
	return corresponding_s1_bins;
}


void S2Band::fill(int bin, double w) {if(bin>=0) {count+=w; counts[bin]+=w;} }

bool S2Band::printIt(int level){
  if(level<1) return true;
  cout<<getName()<<", count:"<<getContent()<<endl;  
  if(level>1) {
    printVector(counts);
    if(level>2){
      printVector(S2overS1Min,"Min S2/S1");
      printVector(S2overS1Max,"Max S2/S1");
    }
    cout<<endl;
  }
  return true;
}

XeSpectrum* S2Band::newSliceXeSpectrum(){
  XeSpectrum *sp=new XeSpectrum("S1 in "+getName(),nSlices,0.,nSlices,counts);
  attach(sp);
  return sp;
}
  

TabulatedDist* S2Band::newSliceXeDist(){
  XeSpectrum* sp=newSliceXeSpectrum();
  TabulatedDist *d= new TabulatedDist(sp);
  delete sp;
  attach(d);
  return d;
}

TH1F* S2Band::newHistogramOfSliceSpectrum(int plot){
  XeSpectrum *sp=newSliceXeSpectrum();
  TH1F* h=sp->newHistogram("Slice in "+getName());
  drawHist(h,plot,LABEL_SLICE,LABEL_EVT);
  delete sp;
  return h;
}

TH1F* S2Band::newHistogramOfSliceDistribution(int plot){
  TabulatedDist *dist=newSliceXeDist();
  TH1F* h=dist->newHistogram("Slice distribution in "+getName());
  drawHist(h,plot,LABEL_SLICE,LABEL_NORMTO1);
  delete dist;
  return h;
}

XeSpectrum* S2Band::newS1XeSpectrum(int nb, double s1Min,double s1Max){
  if(s1Min<0.) s1Min=getS1Min();
  if(s1Max<0.) s1Max=getS1Max();
  if(nb<=0) nb=(int)(s1Max-s1Min+.00001);
  XeSpectrum* sp= new XeSpectrum("S1 in "+getName(),nb,s1Min,s1Max,S1,VALUES);
  attach(sp);
  return sp;
}
  

TabulatedDist* S2Band::newS1XeDist(int nb,double s1Min,double s1Max){
  XeSpectrum* sp=newS1XeSpectrum(nb,s1Min,s1Max);
  TabulatedDist *d= new TabulatedDist(sp);
  delete sp;
  attach(d);
  return d;
}

TH1F* S2Band::newHistogramOfS1Spectrum(int pl,int nb,double s1Min,double s1Max){
  XeSpectrum *sp=newS1XeSpectrum(nb,s1Min,s1Max);
  TH1F* h=sp->newHistogram("S1 in "+getName());
  drawHist(h,pl,LABEL_S1);
  delete sp;
  return h;
}

TH1F* S2Band::newHistogramOfS1Distribution(int plot, int nb
              ,double s1Min,double s1Max){
  TabulatedDist *dist=newS1XeDist(nb,s1Min,s1Max);
  TH1F* h=dist->newHistogram("S1 distribution in "+getName());
  drawHist(h,plot,LABEL_S1,LABEL_NORMTO1);
  delete dist;
  return h;
}

void S2Band::printS1(int maxEvents){ 
  //  print the details with at most "maxEvents" entries
  int n=S1.size();
  cout<<getName()<<" has "<<n<<" entries"<<endl;
  if(n<=maxEvents ){
    for(int i=0;i<n;i++) {
      printf("%7.2f ",S1[i]);
      if(i%10==9) cout<<endl;
    }
    if(n%10!=0) cout<<endl;
  }
}

void S2Band::drawCellContent(double cMin,double cMax,int zMode){

  double x0=0.;
  double y0=0.;
  double x1=0.;
  double y1=0.;

  for(int s=0;s<nSlices;s++){
   if(getS1S2Mode()==BAND_VS_SLICE){
     y0=sequence;
     y1=y0+1.;
     x0=s;
     x1=x0+1.;
   }
   else { 
     x0=s1Bins->getLowerEdge(s);
     x1=s1Bins->getUpperEdge(s);
     y0=S2overS1Min[s];
     y1=S2overS1Max[s];
   }
   int c=indexInterval(getContent(s),cMin,cMax,gStyle->GetNumberOfColors(),zMode);
   if(c>=0) S1S2Display::drawBox(x0,x1,y0,y1,c,true);
  }
}

void S2Band::extendS1Range(double s1Min, double s1Max){
  if(s1Min>getS1Min()){
    cout<<"Can't extend lower range of "<<getName()<<" to "<<s1Min
             <<" while actual value is "<<getS1Min()<<endl;
    return;
  }
  if(s1Max<getS1Max()){
    cout<<"Can't extend upper range of "<<getName()<<" to "<<s1Max
             <<" while actual value is "<<getS1Max()<<endl;
    return;
  }
  s1Bins->extend(s1Min,s1Max);
}


void S2Band::normalize(double f) {
  for(int s=0;s<nSlices;s++) {
    counts[s]=f*counts[s];
  }
  count*=f;
}


void S2Band::drawBandContent(double cMin,double cMax,int zMode){
  int c=indexInterval(count,cMin,cMax,gStyle->GetNumberOfColors(),zMode);
  for(int s=0;s<nSlices;s++){
    double x0=s1Bins->getLowerEdge(s);
    double x1=s1Bins->getUpperEdge(s);
    double y0=S2overS1Min[s];
    double y1=S2overS1Max[s];
    if(c>=0) S1S2Display::drawBox(x0,x1,y0,y1,c,true);
  }
}

void S2Band::draw() {draw(UPPER_EDGE);draw(LOWER_EDGE);}

void S2Band::draw(int lu) {
  vector<double> vx;
  vector<double> vy;
  for(int s=0;s<nSlices;s++){
    double x0=s1Bins->getLowerEdge(s);
    double x1=s1Bins->getLowerEdge(s+1);
    double y= lu==UPPER_EDGE ? S2overS1Max[s]:S2overS1Min[s];
    for(int i=0;i<5;i++){
      float x=x0+(x1-x0)*i*.25;
      vx.push_back(x);
      vy.push_back(x*y);
    }
  } 
  int n=vx.size();
  XeGraph *g=new XeGraph("",n,&(vx[0]),&(vy[0]));
  g->setLineWidth(1);
  g->setLineColor(kGray);
  g->drawS1S2("L");
}

int    S2Band::getSequence()         {return sequence;}
double S2Band::getContent(double s1) {return getContent(s1Bins->getBin(s1));}
double S2Band::getContent(int s)     {return s>=0 && s<nSlices? counts[s]:0.;}
double S2Band::getContent()          {return count;}
vector<double>* S2Band::getContents(){return &counts;}

double  S2Band::getMaxCount(){
  double m=0.;
  for(int s=0;s<nSlices;s++) m=max(m,counts[s]);
  return m;
}

double  S2Band::getMinCount(){
  double m=VERY_LARGE;
  for(int s=0;s<nSlices;s++) m=min(m,counts[s]);
  if(m==VERY_LARGE) m=0;
  return m;
}

double  S2Band::getMinCountButZero(){
  double m=VERY_LARGE;
  for(int s=0;s<nSlices;s++) if(counts[s]>0.) m=min(m,counts[s]);
  if(m==VERY_LARGE) m=0;
  return m;
}

double  S2Band::minX()              {return s1Bins->getLowerEdge();}
double  S2Band::maxX()              {return s1Bins->getUpperEdge();}
double  S2Band::getS1Min()          {return s1Bins->getLowerEdge();}
double  S2Band::getS1Max()          {return s1Bins->getUpperEdge();}
double  S2Band::minY()              {return getLowerEdge();}
double  S2Band::maxY()              {return getUpperEdge();} 
double  S2Band::getUpperEdge()      {return S2overS1Max[nSlices-1];}
double  S2Band::getUpperEdge(int s) {return S2overS1Max[s];}
double  S2Band::getLowerEdge()      {return S2overS1Min[0];}
double  S2Band::getLowerEdge(int s) {return S2overS1Min[s];}
vector <double> S2Band::getUpperEdges()     {return S2overS1Max;}
vector <double> S2Band::getLowerEdges()     {return S2overS1Min;}
double  S2Band::getWidth()          {
  return s1Bins->getUpperEdge()-s1Bins->getLowerEdge();
}

void S2Band::setLowerEdge(int s, double lowEdge) {
   if(s < nSlices && s >= 0)
	S2overS1Min[s] = lowEdge;
   else cout << "S2Band::setLowerEdge  - ERROR : slice number " << s << " is out of range, cannot set edge." << endl;
}
void S2Band::setUpperEdge(int s, double upEdge)  {
   if(s < nSlices && s >= 0)
   	S2overS1Max[s] = upEdge;
   else cout << "S2Band::setUpperEdge  - ERROR : slice number " << s << " is out of range, cannot set edge." << endl;

}

double  S2Band::getWidth(int s)     {
  return s1Bins->getUpperEdge(s)-s1Bins->getLowerEdge(s);
}

void S2Band::fillUniformly(double w){
  reset();
  double density=w/getWidth();
  for(int s=0;s<nSlices;s++) fill(s,density*getWidth(s)); 
}

void S2Band::fillAccordingToS1Dist(double w, double *dist){
  reset();
  for(int s=0;s<nSlices;s++) fill(s,w*dist[s]); 
}

void S2Band::reset() {
  count=0.;
  for(int s=0;s<nSlices;s++) counts[s]=0.;
  S1.clear();
}

void S2Band::setS2overS1Limits(vector<double>* S2overS1Mi
                              ,vector<double>* S2overS1Ma){
  S2overS1Min=*S2overS1Mi;
  S2overS1Max=*S2overS1Ma;
}

S2Band::~S2Band(){}

//------------------ All bands together ---------------------------------

S1S2Bands::~S1S2Bands() {
  deleteWithPointer(s1Bins);
  int nb=bands.size();
  for(int b=0;b<nb;b++)  deleteWithPointer(bands[b]);
  int ns=slices.size();
  for(int s=0;s<ns;s++)  deleteWithPointer(slices[s]);
  counts.clear();
}

S1S2Bands::S1S2Bands() : DetectorComponent(), S1S2Object(),RunComponent()  {}
S1S2Bands::S1S2Bands(string nam, XeRun* r,int nB, int nS1,double s1Min
            ,double s1Max) : DetectorComponent(), S1S2Object() ,RunComponent(){
  if(nam=="") nam=getTheName(nB,nS1);
  setName(nam);
  setRun(r);
  s1Bins=new EquidistantBins("S1",nS1,s1Min,s1Max);
  instantiateBandsAndSlices(nam,nB);
}

S1S2Bands::S1S2Bands(string nam, XeRun* r,int nB,S1S2Data* s1s2,int nS1) 
  :  S1S2Object(), RunComponent() {
  if(nam=="") nam="Bands of "+s1s2->getName();
  setName(nam);
  setRun(r);
  s1Bins=s1s2->newEquiContentS1Bins(nS1);
  instantiateBandsAndSlices(nam,nB);
}

S1S2Bands::S1S2Bands(string nam, XeRun* r,int nB, XeBins* s1B) 
  :  S1S2Object(nam) {
  if(nam=="") nam=getTheName(nB,s1B->getNBins());
  run=r;
  s1Bins=s1B;
  instantiateBandsAndSlices(nam,nB);
}

S1S2Bands::S1S2Bands(string nam, S1S2Bands* bds, bool copy) : S1S2Object(nam) {
  run=bds->getRun();
  nBands=bds->getNBands();
  nSlices=bds->getNSlices();
  s1Bins=new XeBins((const XeBins&)*(bds->getS1Bins()));
  for(int b=0;b<nBands;b++){
    bands.push_back(new S2Band((const S2Band&)*(bds->getBand(b))));
  }
  for(int s=0;s<nSlices;s++){
    slices.push_back(new S1Slice((const S1Slice&)*(bds->getSlice(s))));
  }
  if(nam=="") nam=bds->getName();
  setNames(nam);
  establishPointersToCounts();
  if(!copy) reset();
}

S1S2Bands::S1S2Bands(string nam, S1S2Bands bds, bool copy) : S1S2Object(nam) {
  run=bds.getRun();
  nBands=bds.getNBands();
  nSlices=bds.getNSlices();
  s1Bins=new XeBins((const XeBins&)*(bds.getS1Bins()));
  for(int b=0;b<nBands;b++){
    bands.push_back(new S2Band((const S2Band&)*(bds.getBand(b))));
  }
  for(int s=0;s<nSlices;s++){
    slices.push_back(new S1Slice((const S1Slice&)*(bds.getSlice(s))));
  }
  if(nam=="") nam=bds.getName();
  setNames(nam);
  establishPointersToCounts();
  if(!copy) reset();
}

void S1S2Bands::establishPointersToCounts(){
  counts.clear();
  for(int b=0;b<nBands;b++) counts.push_back(bands[b]->getContents());
}

bool S1S2Bands::isDrawable(int ) {return true;}

void S1S2Bands::cumulate(S1S2Bands* original){
  // *this and *original assumed of same definition
  for(int b=0;b<nBands;b++){
    S2Band *bo=original->getBand(b);
    for(int s=0;s<nSlices;s++){
      double w=bo->getContent(s);
      if(b>0) w += bands[b-1]->getContent(s);
      fill(s,b,w);
    }
  }
}

XeSpectrum* S1S2Bands::newBandXeSpectrum(){
  XeSpectrum* spectrum = new XeSpectrum("Band in "+getName(),nBands,.0,nBands);
  double *v=spectrum->getSpectrum();
  for(int b=0;b<nBands;b++) v[b]=getBandContent(b);
  attach(spectrum);
  return spectrum;
}

TH1F* S1S2Bands::newHistogramOfBandSpectrum(int plot){
  XeSpectrum* sp=newBandXeSpectrum();
  TH1F* h=sp->newHistogram();
  delete sp;
  drawHist(h,plot,LABEL_BAND);
  return h;
}

TH1F* S1S2Bands::newHistogramOfBandDistribution(int plot){
  XeSpectrum* sp=newBandXeDist();
  TH1F* h=sp->newHistogram();
  delete sp;
  drawHist(h,plot,LABEL_BAND,LABEL_NORMTO1);
  return h;
}

TabulatedDist* S1S2Bands::newBandXeDist(){
  XeSpectrum* sp=newBandXeSpectrum();
  TabulatedDist* dist=new TabulatedDist(sp);
  delete sp;
  attach(dist);
  return dist;
}


string S1S2Bands::bandName(int b) {
  if(b>=0 && b<nBands) return "band"+formatI(b,3);
  else if(b==ALL) return "all bands";
  return "band "+UNDEFINED_STRING;
}

pair<int,int>S1S2Bands::whichBands(int band){
  if(band>=0 && band<nBands) return pair<int,int>(band,band);
  else if(band==ALL)         return pair<int,int>(0,nBands-1);
  cout<<"Invalid band index:"<<band<<endl;
  return pair<int,int> (NONE,NONE-1);
}

XeSpectrum* S1S2Bands::newSliceXeSpectrum(int band){
  trace(ENTRY,"S1S2Bands::newSliceXeSpectrum");
  pair<int,int> p=whichBands(band); 
  XeSpectrum* spectrum = new XeSpectrum("Slice in "+bandName(band)
                                       ,nSlices,.0,nSlices);
  for(int b=p.first;b<=p.second;b++) {
    XeSpectrum* sp=bands[b]->newSliceXeSpectrum();
    spectrum->add(sp);
    delete sp;
  }
  trace(EXIT,"S1S2Bands::newSliceXeSpectrum");
  attach(spectrum);
  return spectrum;
}

TH1F* S1S2Bands::newHistogramOfSliceSpectrum(int band,int plot){
  XeSpectrum* sp=newSliceXeSpectrum(band);
  TH1F* h=sp->newHistogram();
  delete sp;
  drawHist(h,plot,LABEL_SLICE);
  return h;
}

TH1F* S1S2Bands::newHistogramOfSliceDistribution(int band,int plot){
  XeSpectrum* sp=newSliceXeDist(band);
  TH1F* h=sp->newHistogram();
  delete sp;
  drawHist(h,plot,LABEL_SLICE,LABEL_NORMTO1);
  return h;
}

TabulatedDist* S1S2Bands::newSliceXeDist(int band){
  XeSpectrum* sp=newSliceXeSpectrum(band);
  TabulatedDist* dist=new TabulatedDist(sp);
  delete sp;
  attach(dist);
  return dist;
}

XeSpectrum* S1S2Bands::newS1XeSpectrum(int band,int nb,double s1Mi,double s1Ma){
  pair<int,int> p=whichBands(band); 
  if(s1Mi==AUTOMATIC) s1Mi=getS1LowerEdge();
  if(s1Ma==AUTOMATIC) s1Ma=getS1UpperEdge();
  if(nb==AUTO) nb=(int)(s1Ma-s1Mi+.00001);
  XeSpectrum* spectrum = new XeSpectrum("S1 in "+getName(),nb,s1Mi,s1Ma);
  for(int b=p.first;b<=p.second;b++) {
    XeSpectrum* sp=bands[b]->newS1XeSpectrum(nb,s1Mi,s1Ma);
    spectrum->add(sp);
    delete sp;
  }
  attach(spectrum);
  return spectrum;
}

TH1F S1S2Bands::getHistogramOfS1Distro(int from_band, int to_band){
   //simple implementation that returns histo filled with S1 distro for the integral of the specified bands
  	
    if(to_band == 999) to_band = from_band;

    TString name_and_title = "S1 distro Band "+ formatI(from_band)+ "-"+ formatI(to_band) + " "+ getName();

    TH1F h(name_and_title , name_and_title,getNSlices(),getS1LowerEdge(), getS1UpperEdge());

    if(from_band > to_band || from_band < 0 || from_band > getNBands() || to_band > getNBands() || to_band < 0) return h;


    for(int b=from_band ; b <=to_band; b++){
       for(int s=0; s < getNSlices(); s++) {
	  double cell_content  = getContent(s,b);
	  h.Fill(getS1Center(s), cell_content);	
	  if(debugLevel > 1) cout << "Band " << b << " Slice " << s << " content " << cell_content << endl;
       }
    }

    h.GetXaxis()->SetTitle("S1 [p.e]");
    h.GetYaxis()->SetTitle("Entries / PE");
    
    //NOTE: filling with weight calls automatically the storage of w^2 for errors.
    //Errors here are set as sqrt of counts, which is good only for data.
    h.Sumw2(kFALSE);
    
    return h;

}

TH1F* S1S2Bands::newHistogramOfS1Spectrum(int band,int plot,int nb,double s1Min
                                         , double s1Max) {
  XeSpectrum* sp=newS1XeSpectrum(band,nb,s1Min,s1Max);
  TH1F* h=sp->newHistogram();
  delete sp;
  drawHist(h,plot,LABEL_S1);
  return h;
}

TH1F* S1S2Bands::newHistogramOfS1Distribution(int band,int plot,int nb
                            ,double s1Min, double s1Max){
  XeSpectrum* sp=newS1XeDist(band,nb,s1Min,s1Max);
  TH1F* h=sp->newHistogram();
  delete sp;
  drawHist(h,plot,LABEL_S1,LABEL_NORMTO1);
  return h;
}

TabulatedDist* S1S2Bands::newS1XeDist(int band,int n,double s1Min,double s1Max){
  XeSpectrum* sp=newS1XeSpectrum(band,n,s1Min,s1Max);
  TabulatedDist* dist=new TabulatedDist(sp);
  delete sp;
  attach(dist);
  return dist;
}

XeMultiGraph* S1S2Bands::newMultiGraphOfSlice(string nam,XeSpectrum** spectra
                                             ,int plot){
  if(!run->checkAll()) return NULL;
  XeMultiGraph *mg=new XeMultiGraph(nam);
  for(int b=0;b<nBands;b++){
    XeGraph* g=spectra[b]->newGraph(); 
    g->setBandLegend(b);
    mg->add(g);
  }
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

XeMultiGraph* S1S2Bands::newMultiGraphOfSlice(string nam,TabulatedDist** dist
                                             ,int plot){
  if(!run->checkAll()) return NULL;
  XeMultiGraph *mg=new XeMultiGraph(nam);
  for(int b=0;b<nBands;b++){
    XeGraph* g=dist[b]->newGraph(); 
    g->setLegend(bandName(b));
    mg->add(g);
  }
  mg->setRainbowColors(ASYMMETRIC);
  mg->drawIt(plot);
  return mg;
}

XeMultiGraph* S1S2Bands::newMultiGraphOfSliceSpectrum(int plot){
  vector<XeSpectrum*> spectra;
  for(int b=0;b<nBands;b++){
    XeSpectrum *sp=bands[b]->newSliceXeSpectrum();
    spectra.push_back(sp);
  } 
  XeMultiGraph *mg=newMultiGraphOfSlice(getName()+"- Slice spectrum"
                                           , &(spectra[0]));
  for(int b=0;b<nBands;b++) delete spectra[b];
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

vector<TGraph> S1S2Bands::getTGraphOfBands(bool isS2_vs_S1){


   // the vector filled with all the bands
   vector<TGraph> bands_vector;

   for(int b=0;b<nBands;b++) {
	// temporary band holder
	TGraph tempBand; 

	for(int s=0; s<nSlices; s++){
		if(isS2_vs_S1) tempBand.SetPoint(s, getS1Center(s), getS2overS1LowerEdge(b,s) * getS1Center(s) );
		else tempBand.SetPoint(s, getS1Center(s), log10(getS2overS1LowerEdge(b,s)) );
	}	

	bands_vector.push_back(tempBand);
   }
   
   // last band upper edge
   TGraph tempBand;
	for(int s=0; s<nSlices; s++){
		if(isS2_vs_S1) tempBand.SetPoint(s, getS1Center(s), getS2overS1UpperEdge(nBands-1,s) * getS1Center(s) );
		else tempBand.SetPoint(s, getS1Center(s), log10(getS2overS1UpperEdge(nBands-1,s)) );
	}	

   bands_vector.push_back(tempBand);
   

   return bands_vector;

}

XeMultiGraph* S1S2Bands::newMultiGraphOfSliceDistribution(int plot){
  vector<TabulatedDist*> dist;
  for(int b=0;b<nBands;b++){
    TabulatedDist *sp=bands[b]->newSliceXeDist();
    dist.push_back(sp);
  } 
  XeMultiGraph *mg=newMultiGraphOfSlice(getName()+"- Slice distribution"
                                           ,&(dist[0]));
  for(int b=0;b<nBands;b++) delete dist[b];
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

TH2F* S1S2Bands::new2DHistogram(int plot){
  TH2F* h=new TH2F(Name.c_str(),Name.c_str(),nSlices,0.,nSlices
                  ,nBands,0.,nBands);
  for(int b=0;b<nBands;b++){
    S2Band* band=getBand(b);
    double y=b+.1;
    for(int s=0;s<nSlices;s++) h->Fill(s+.1,y,band->getContent(s));
  }
  drawHist(h,LABEL_SLICE,LABEL_BAND,plot);
  return h;
}

void S1S2Bands::instantiateBandsAndSlices(string nam,int nB){

  int runNumber = getRunNumber();

  fBottom=NULL;

  if (runNumber==12)  fBottom=new TF1("fBottom","-2.0e-02 * x^3 + 8.8e-01 * x^2 + 8.8e+00 * x + 3.9e+01",3,30);
  else if (runNumber==10)  fBottom=new TF1("fBottom","-2.3e-02 * x^3 + 1.0e+00 * x^2 + 1.1e+01 * x + 4.4e+01",3,30);
  else if (runNumber==8)   fBottom=new TF1("fBottom","-4.1e-02 * x^3 + 2.2e+00 * x^2 + -1.4e+01 * x + 1.9e+02",3,30);

  nBands=nB;
  nSlices=s1Bins->getNBins(); 
  if(debugLevel>0){
    cout<<"Establishing "<<nBands<<" bands and "<<nSlices<<" slices for "
        <<getName()<<endl;
  }
  for(int b=0;b<nBands;b++) bands.push_back(new S2Band(b,this));
  for(int s=0;s<nSlices;s++) slices.push_back(new S1Slice(s,this));
  establishPointersToCounts();
  setNames(nam);
}

int       S1S2Bands::getNBands()            {return nBands;}
int       S1S2Bands::getNSlices()           {return nSlices;}
double    S1S2Bands::getS1UpperEdge()       {return s1Bins->getUpperEdge();}
double    S1S2Bands::getS1UpperEdge(int s)  {return s1Bins->getUpperEdge(s);}
double    S1S2Bands::getS1LowerEdge()       {return s1Bins->getLowerEdge();}
double    S1S2Bands::getS1LowerEdge(int s)  {return s1Bins->getLowerEdge(s);}
double    S1S2Bands::getContent(int s,int b){return bands[b]->getContent(s);}
double    S1S2Bands::getBandContent(int b)  {return bands[b]->getContent(); }
double*   S1S2Bands::getS1UpperEdges()      {return s1Bins->getUpperEdges();}
double*   S1S2Bands::getS1LowerEdges()      {return s1Bins->getLowerEdges();}
XeBins*   S1S2Bands::getS1Bins()            {return s1Bins;}
S1Slice*  S1S2Bands::getSlice(int s)        {return slices[s];}
S2Band*   S1S2Bands::getBand(int b)         {return bands[b];}

double S1S2Bands::getS1Width(int s){return getS1UpperEdge(s)-getS1LowerEdge(s);}
double S1S2Bands::getS1Width()     {return getS1UpperEdge() -getS1LowerEdge() ;}

double S1S2Bands::getS1Center(int s){
  return 0.5*(getS1LowerEdge(s)+getS1UpperEdge(s));
}

double S1S2Bands::getCumulatedBandContent(int b){
  double c=0.;
  for(int i=0;i<=b;i++) c += getBandContent(i);
  return c;
}

double S1S2Bands::getS2overS1UpperEdge(int b)     {
  return bands[b]->getUpperEdge();
}

double S1S2Bands::getS2overS1UpperEdge(int b,int s){
  return bands[b]->getUpperEdge(s);
}
double  S1S2Bands::getS2overS1LowerEdge(int b)     {
  return bands[b]->getLowerEdge();
}
double S1S2Bands::getS2overS1LowerEdge(int b,int s){
  return bands[b]->getLowerEdge(s);
}
vector<double> S1S2Bands::getS2overS1UpperEdges(int b)    {
  return bands[b]->getUpperEdges();
}
vector<double> S1S2Bands::getS2overS1LowerEdges(int b)    {
  return bands[b]->getLowerEdges();
}

void S1S2Bands::reset() {
  for(int b=0;b<nBands;b++) bands[b]->reset();
  for(int s=0;s<nSlices;s++) slices[s]->reset();
}

string  S1S2Bands::getTheName(int nB, int nS){
  char title[50];
  sprintf(title,"[%d*%d ]",nB,nS);
  return string(title);
}

void S1S2Bands::setNames(string nam){
  setName(nam);
  for(int b=0;b<nBands;b++)  bands[b]->setName(S2Band::getTheName(b,this));
  for(int s=0;s<nSlices;s++) slices[s]->setName(S1Slice::getTheName(s,this));
}

void S1S2Bands::extendS1Range(double s1Min, double s1Max){
  if(s1Min>getS1UpperEdge()){
    cout<<"Can't extend lower range of "<<getName()<<" to "<<s1Min
             <<" while actual value is "<<getS1UpperEdge()<<endl;
    return;
  }
  if(s1Max<getS1LowerEdge()){
    cout<<"Can't extend upper range of "<<getName()<<" to "<<s1Max
             <<" while actual value is "<<getS1UpperEdge()<<endl;
    return;
  }
  getSlice(0)->setS1Min(s1Min);
  getSlice(getNSlices()-1)->setS1Max(s1Max);
  for(int b=0;b<nBands;b++) getBand(b)->extendS1Range(s1Min,s1Max);
}

void S1S2Bands::add(S1S2Bands *toBeAdded){
	//NOTE: since S1S2Bands doesn't know how to make a copy (there is a copy constructor)
        //one cannot pass input by value but only by reference.

	//safety checks
	if(toBeAdded->getNBands() != nBands  || toBeAdded->getNSlices() != nSlices) 
		cout <<"S1S2Bands::add() - ERROR : Bands '" << getName() << "' cannot be added with '" 
		     << toBeAdded->getName() << "' since they are different" << endl;

	else {
		for(int s=0; s < nSlices; s++){
	   	   for(int b=0; b < nBands; b++){
			fill(s,b, toBeAdded->getContent(s,b) ); 
		   }
		}
	}
        establishPointersToCounts();
}


void S1S2Bands::compute(TH2F h2d_reference, bool isLog){

  //----------------- NOTES -------------------------------- //
  // *) Bands are defined in S2/S1  VS S1 space.     	     //
  // *) the histogram is meant to be in log(s2/s1)           //
  //    versus S1 or S2 VS s1.				     //
  //---------------------------------------------------------//
  
  // NOTE: there is a problem here, with hagar histos lower band is always cutting away signal
  // due to the hard low band cut defined in s2/s1, this can be fixed defining a more suitable
  // lowest band cut   FIXME

  if(debugLevel>0){
    cout<<"Computing bands for "<<getName()<<endl;
  }

  double min_s2_over_s1 = run->getLowestS2OverS1();  // default is 20   
  double max_s2_over_s1 = run->getHighestS2OverS1(); // default is 500 

  //check histogram consistency with XEPHYR standards
  if( h2d_reference.GetYaxis()->GetBinWidth(1) > 0.001) cout << "S1S2Bands::compute() - WARNING : Histo may have too large Y binning." << endl;
	
  // Loop over S1 slices
  for(int slice = 0; slice < nSlices ; slice++)
  {
	double s1    = getS1Center(slice);		// s1 value of the current slice
	int    bin_x = h2d_reference.GetXaxis()->FindBin(s1);	
	
	//check histogram consistency with XEPHYR standards
	if( fabs(h2d_reference.GetXaxis()->GetBinLowEdge(bin_x) - getS1LowerEdge(slice) )  > 1.E-10   //there is some rounding issue here between TAxis::GetBinLowEdge and getS1LowerEdge()  FIXME
	    || fabs(h2d_reference.GetXaxis()->GetBinUpEdge(bin_x) - getS1UpperEdge(slice)) > 1.E-10 ) 
	{
	   cout << "S1S2Bands::compute() - ERROR : TH2F has not the standard s1 binning. BAD BANDS!"<< endl;
	   cout << "S1S2Bands::compute() - ERROR : Xephyr S1 low Edge " << getS1LowerEdge(slice)  << "    TH2F S1 low  "<< h2d_reference.GetXaxis()->GetBinLowEdge(bin_x)<< endl;
	   cout << "S1S2Bands::compute() - ERROR : Xephyr S1 up  Edge " << getS1UpperEdge(slice)  << "    TH2F S1 low  "<< h2d_reference.GetXaxis()->GetBinUpEdge(bin_x)<< endl;
	   cout << "S1S2Bands::compute() - ERROR : difference  "  << getS1UpperEdge(slice) - h2d_reference.GetXaxis()->GetBinUpEdge(bin_x) << endl;
	   exit(100);
	   break;
	}	
		
	//-------  Computing bands edges for given slice  -------//
    	  //getting y min and max, the histo is meant to be in log(s2/s1) or in S2 VS S1
	  int    bin_y_min = h2d_reference.GetYaxis()->FindBin(log10(min_s2_over_s1));	
		if(!isLog)   {
            // Override minimum S2 
            // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:hagar:pl_news_run12#pl_news_run12_26_5_2016_3nsigma_cut
	    double s2 = min_s2_over_s1 * s1;
            if(fBottom != NULL)  s2 = fBottom->Eval(s1);
            bin_y_min = h2d_reference.GetYaxis()->FindBin(s2);
        }
	  int    bin_y_max = h2d_reference.GetYaxis()->FindBin(log10(max_s2_over_s1));
		if(!isLog)   bin_y_max = h2d_reference.GetYaxis()->FindBin(max_s2_over_s1 * s1);
	 
	  int histoBinMax =  h2d_reference.GetYaxis()->GetNbins(); 
          
	  //sanity checks, NOTE: in case required Y range is out of histo bin_y_max = histoBinMax+1
	  if(bin_y_min == 0 || bin_y_max >= histoBinMax ) cout <<"S1S2Bands::compute() - WARNING required Y range is larger than histo Y range, possible problem. S1 slice " << slice << endl;
	  if(bin_y_min >= histoBinMax || bin_y_max == 0) { cout <<"S1S2Bands::compute() - ERROR Histo out of range, BAD BANDS." << endl; break;}

	  double equo_density = h2d_reference.Integral(bin_x, bin_x, bin_y_min, bin_y_max) / (double)nBands;

	  if(debugLevel>0) cout << "S1S2Bands::compute() - DEBUG : for slice S1= " << s1 <<" bin_y_min " << bin_y_min << " bin_y_max " << bin_y_max << " density " << equo_density << endl;
	  
          vector <double> low_edges;
	  if(isLog) low_edges.push_back(pow(10,h2d_reference.GetYaxis()->GetBinLowEdge(bin_y_min))); //remember histo is in log(s2/s1), bands are in s2/s1.
	  else      low_edges.push_back(h2d_reference.GetYaxis()->GetBinLowEdge(bin_y_min) / s1); //case histo is in S2 VS s1

          double temp_density = 0;
	  int counterBands = 0;
	  //Loop over Y bins to get low and up edges of equal density bands
	    for(int bin_itr = bin_y_min; bin_itr <= bin_y_max; bin_itr++)
	    {
		temp_density += h2d_reference.GetBinContent(bin_x , bin_itr);

		//safety check
		if(h2d_reference.GetBinContent(bin_x , bin_itr) / equo_density > 0.5 ) cout << "S1S2Bands::compute()  - WARNING : Histo " << h2d_reference.GetName() 
				 							    << " Y axis has too large bin width ( > 50% density) for bin X "<< bin_x << " Y " << bin_itr<< endl;

		//to minimize the jitter due to finite size of the Y axis bin we integrate Y until Y_bins_integral >= i * density
		if(temp_density >= equo_density + equo_density * counterBands) {
		   //fill edge
		   if(isLog) low_edges.push_back(pow(10,h2d_reference.GetYaxis()->GetBinUpEdge(bin_itr)));
		   else low_edges.push_back(h2d_reference.GetYaxis()->GetBinUpEdge(bin_itr) / s1);

		   if(debugLevel>0) cout << "S1S2Bands::compute()  -  DEBUG : Slice S1= " << s1 << " Band " << counterBands << " Density check: computed " << equo_density + equo_density * counterBands<< "  obtained " << temp_density  << endl;
		   counterBands++;
		}
		if(counterBands >= nBands -1) break;  // we are computing the lower edges of bands
						      // we already arrived to the low edge of the last band.
						      // note that here I fill the band lower edge as the 
						      // previous band upper edge.
	    }

	   unsigned int U_nBands = nBands; // just to fix usual compiler complaint
	   if(low_edges.size() != U_nBands) { cout <<"S1S2Bands::compute() - ERROR : wrong number of final bands" << endl; break;}
	
	  //extend the vector of edges with the last band edge (used in initializzation of S1slice)
	  vector <double> edges = low_edges;
	  if(isLog) edges.push_back(max_s2_over_s1);	// the maximum edge is extended to the maximum default
	  else      edges.push_back(max_s2_over_s1); 
	
	  // vector used in initializzation of S2Band
	  vector <double> up_edges;
	  for(int e=0; e < nBands -1 ; e++) up_edges.push_back(low_edges[e+1]);
	  if(isLog) up_edges.push_back(max_s2_over_s1);  // the maximum edge is extended to the maximum default
	  else      up_edges.push_back( max_s2_over_s1);

	//--------------------------------------------------------//

 	//initializing the S2bins of  member "slices" of S1S2Bands. 
	XeBins s2overs1Bins(getName() + " slice " , edges);
	slices[slice]->setS2overS1Bins(&s2overs1Bins);
  	if(debugLevel>0){
		cout << "--------------------- S1S2Bands::compute() - DEBUG : PRINTING SLICE " << slice << "  --------------------------------" << endl;
		slices[slice]->printIt(10);
	}	
	
	//Initializing "S2band" member of S1S2Bands
	for(int band_itr=0; band_itr < nBands; band_itr++) {
		bands[band_itr]->setUpperEdge(slice,  up_edges[band_itr]);
		bands[band_itr]->setLowerEdge(slice, low_edges[band_itr]);
	}	 
   }

}


void S1S2Bands::compute(S1S2Data* expectedSignal){
  if(debugLevel>0){
    cout<<"Computing signal for "<<getName()<<endl;
  }
  vector<vector<double>*> edges;    // accessed as "edges[band][slice]"
  
  for(int b=0;b<=nBands;b++) {
    if(b<nBands) bands[b]->reset();
    vector<double>* vd=new vector<double>(nSlices);
    edges.push_back(vd);
  }
  for(int s=0;s<nSlices;s++) slices[s]->reset();
  int ne=expectedSignal->getNEvents();
  if(debugLevel>0){
    cout<<"Looping over "<<ne<<" in "<<expectedSignal->getName()<<endl;
  }
  for(int e=0;e<ne;e++){
    double s1=expectedSignal->getS1(e);
    int s=s1Bins->getBin(s1);
    if(s>=0) slices[s]->fillForStudy(s1,expectedSignal->getS2(e));
  }
  if(debugLevel>0){
    cout<<"Looping over "<<nSlices<<" slices "<<endl;
  }
  double rMin=run->getLowestS2OverS1();
  double rMax=run->getHighestS2OverS1();
  for(int s=0;s<nSlices;s++) {
    vector<double>* values=slices[s]->getS2overS1Vector();
    if(debugLevel>0){
      cout<<"Finding the "<<nBands<<" equicontent bins for "
          <<slices[s]->getName()<<", which has "<<values->size()<<" entries"
          <<endl;
    }
    EquiContentBins* S2overS1bins=new EquiContentBins(
                           "S2 in "+expectedSignal->getName(),nBands,values); 
    S2overS1bins->extend(rMin,rMax);
    if(debugLevel>0){
      cout<<"Filling the edge tables"<<endl;
    }
   
    for(int b=0;b<=nBands;b++) (*edges[b])[s]=S2overS1bins->getLowerEdge(b);
    slices[s]->setS2overS1Bins(S2overS1bins);
  }
  
  for(int b=0;b<nBands;b++) {
    bands[b]->setS2overS1Limits(edges[b],edges[b+1]);
    delete edges[b];
  }
  delete edges[nBands];
}

double S1S2Bands::getMinimumBandSize( bool isLogRequired) {

  double minimum_band_size = VERY_LARGE;

  for(int s=0;s<nSlices;s++) {
	
	double s1 = getS1Center(s);
  	for(int b=0;b<nBands;b++){
	  double low_edge = getS2overS1LowerEdge(b,s) * s1;
	  double up_edge  = getS2overS1UpperEdge(b,s) * s1;
	
	  if(isLogRequired) {
		low_edge = log10(low_edge / s1);
		up_edge  = log10(up_edge / s1 );
	  }

	  if(minimum_band_size > (up_edge - low_edge))  minimum_band_size = (up_edge - low_edge) ;
	}
  }
  return minimum_band_size;

}

void S1S2Bands::reweightTheBands(double * weights){ 
  vector<double> coefs(nBands);
  double sumw=0.;
  for(int b=0;b<nBands;b++) sumw += weights[b];
  double sumo=getTotalContent();
  for(int b=0;b<nBands;b++) {
    coefs[b]=weights[b]/sumw*sumo/getBand(b)->getContent();
    getBand(b)->normalize(coefs[b]);
  }
  for(int s=0;s<nSlices;s++)  getSlice(s)->reweightTheBands(&(coefs[0]));
}

void S1S2Bands::drawSlices(){
  for(int s=0;s<nSlices;s++) slices[s]->draw();
}

void S1S2Bands::draw(){draw(LOG);}
void S1S2Bands::draw(int zMode){
  if(S1S2Display::getMode()==BAND_VS_SLICE) new2DHistogram(zMode);
  else {
    drawBands();
    drawCellContent(zMode);
  }
}

void S1S2Bands::drawTheFrame(){drawCellContentFrame(LOG);}

void S1S2Bands::drawBands(){
  bands[0]->draw(LOWER_EDGE);
  for(int b=0;b<nBands;b++) bands[b]->draw(UPPER_EDGE);
}

void S1S2Bands::drawCellContentFrame(int zMode){
  drawFrameWithZ(getMinCountButZero(),getMaxCount(),zMode);
  gPad->SetLogy(S1S2Display::getDefaultLogMode());
}

void S1S2Bands::drawCellContent(int zMode){
  double cMin=getMinCountButZero();
  double cMax=getMaxCount();
  for(int b=0;b<nBands;b++) bands[b]->drawCellContent(cMin,cMax,zMode);
  drawBands();
}

void S1S2Bands::drawCellContentWithFrame(int zMode){
  drawCellContentFrame(zMode);
  drawCellContent(zMode);
}


void S1S2Bands::produceGraphicOverview(){
  setS1S2Mode(BAND_VS_SLICE);
  openCanvas(3,3);
  subCanvas(0.00,0.97,.35,1.);
  drawCellContentWithFrame(LOG);
  subCanvas(6);
  newHistogramOfSliceSpectrum(ALL,LINEAR);
  nextSubCanvas(); 
  newMultiGraphOfSliceDistribution()->drawWithFrame();
  nextSubCanvas(); 
  newHistogramOfBandSpectrum(LOG);
}

void S1S2Bands::drawBandContentFrame(int zMode){
  drawFrameWithZ(getMinBandContentButZero(),getMaxBandContent(),zMode);
  gPad->SetLogy(S1S2Display::getDefaultLogMode());
}

void S1S2Bands::drawBandContent(int zMode){
  double cMin=getMinBandContentButZero();
  double cMax=getMaxBandContent();
  for(int b=0;b<nBands;b++) bands[b]->drawBandContent(cMin,cMax,zMode);
  drawBands();
}

void S1S2Bands::drawBandContentWithFrame(int zMode){
  drawBandContentFrame(zMode);
  drawBandContent(zMode);
}

bool S1S2Bands::printIt(int level){
  if(!run->checkAll()) return false;
  if(level<1) return true;
  cout<<Name<<" has "<<nBands<<" bands and "<<nSlices<<" slices; total count :"
      <<getTotalContent()<<endl;
  if(level<2) return true;
  cout<<getName()<<endl;
  printBands(level-1);
  printSlices(level-1);
  return true;
}

void S1S2Bands::printSlices(int level){
  if(level<1) return;
  if(!run->checkAll()) return;
  cout<<endl;
  double total=0;
  for(int s=0;s<nSlices;s++) {
    slices[s]->printIt(level);
    if(level > 3) 
	{	
		cout << "--------  Printing also the log10(s2/s1) min/max -----------" << endl;
		slices[s]->printItLog(level);  // print also the log10(s2/s1) min/max 
		cout << "--------  -----------------------------   -----------" << endl;
        }
    total += slices[s]->getContent(); 
  }
  cout<<"Total count :"<<total<<endl;
}

void S1S2Bands::printBands(int level){
  if(level<1) return;
  if(!run->checkAll()) return;
  cout<<endl;
  double total=0;
  for(int b=0;b<nBands;b++) {
    bands[b]->printIt(level);
    total += bands[b]->getContent();
  }
  cout<<"Total count :"<<total<<endl;
}

void S1S2Bands::printS1InBands(int maxEvents){
  cout<<endl;
  for(int b=0;b<nBands;b++) bands[b]->printS1(maxEvents);
}



S1Slice* S1S2Bands::getTheSlice(double s1){
  int s=s1Bins->getBin(s1);
  if(s>=0) return slices[s];
  return NULL;
}

pair<int,int> S1S2Bands::whichSliceAndBand(double s1,double s2){
  if(s1<=0.) return pair<int,int>(-1,-1);
  int b=-1;
  int s=s1Bins->getBin(s1);
  if(s>=0) b=slices[s]->getS2overS1Bins()->getBin(s2/s1);
  return pair<int,int>(s,b);
}

void S1S2Bands::copyBandsToSlices(){
  for(int s=0;s<nSlices;s++){
    S1Slice *slice=slices[s];
    slice->reset();
    for(int b=0;b<nBands;b++) slice->fill(b,bands[b]->getContent(s));
  }
}

void S1S2Bands::fill(double s1,int s,int b,double w){
  if(s>=0) {
    slices[s]->fill(b,w);
    if(b>=0) bands[b]->fill(s1,w);
  }
}

void S1S2Bands::fill(int s,int b,double w){
  if(s>=0) {
    slices[s]->fill(b,w);
    if(b>=0) bands[b]->fill(s,w);
  }
}

void S1S2Bands::fillUniformly(vector<double>& weights){
  if((int)weights.size()!=nBands){
    cout<<"Size mismatch in S1S2Bands::fillUniform;y: "<<weights.size()
        <<" vs "<<nBands<<endl;
    return;
  }
  for(int b=0;b<nBands;b++) bands[b]->fillUniformly(weights[b]);
  copyBandsToSlices();
}

void S1S2Bands::fillAccordingToS1Dist(vector<double>& weights, double* dist){
  if((int)weights.size()!=nBands){
    cout<<"Size mismatch in S1S2Bands::fillUniformly: "<<weights.size()
        <<" vs "<<nBands<<endl;
    return;
  }
  for(int b=0;b<nBands;b++) bands[b]->fillAccordingToS1Dist(weights[b],dist);
  copyBandsToSlices();
}

void S1S2Bands::fillExponential(vector<double>& weights,double slope){
  if((int)weights.size()!=nBands){
    cout<<"Size mismatch in S1S2Bands::fillExponential: "<<weights.size()
        <<" vs "<<nBands<<endl;
    return;
  }
  double overall=exp(-slope*getS1LowerEdge())-exp(-slope*getS1UpperEdge());
  for(int b=0;b<nBands;b++){
    double d=weights[b]/overall;
    for(int s=0;s<nSlices;s++) {
      double c=d*(exp(-slope*getS1LowerEdge(s))-exp(-slope*getS1UpperEdge(s)));
      fill(s,b,c);
    }
  }
}

void S1S2Bands::dispatchS1InEqualBands(double s1, double w){
  if(s1Bins==NULL) {
    cout<<"Major inconsistency "<<endl;
    return;
  }
  int s=s1Bins->getBin(s1);
  if(s<0) return;
  double wn=w/nBands;
  for(int b=0;b<nBands;b++) fill(s1,s,b,wn);
}

double S1S2Bands::getTotalContent(){
  double m=0.;
  for(int b=0;b<nBands;b++) m+=bands[b]->getContent();
  return m;
}

double S1S2Bands::getS1Density(double s1){
  S1Slice* slice=getTheSlice(s1);
  if(slice==NULL) return 0.;
  return slice->getContent()/slice->getS1Width();
}

double S1S2Bands::getContentBelowS1(double s1){
  if(s1<=getS1LowerEdge()) return 0;
  if(s1>=getS1UpperEdge()) return getTotalContent();
  pair<int,double> bf=s1Bins->getBinAndFraction(s1);
  int b=bf.first;
  double count=bf.second*slices[b]->getContent();
  for(int s=0;s<b;s++) count+= slices[s]->getContent();
  return count;
}

void S1S2Bands::normalizeToEvents(double events){
  double ratio=events/getTotalContent();
  normalize(ratio);
}

void S1S2Bands::normalize(double ratio){
  for(int b=0;b<nBands;b++) bands[b]->normalize(ratio);
  for(int s=0;s<nSlices;s++) slices[s]->normalize(ratio);
}

double S1S2Bands::getMaxCount(){
  double m=0.;
  for(int b=0;b<nBands;b++) m=max(m,bands[b]->getMaxCount());
  return m;
}

double S1S2Bands::getMinCount(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) m=min(m,bands[b]->getMinCount());
  if(m==VERY_LARGE) m=.0;
  return m;
}

double S1S2Bands::getMinCountButZero(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) m=min(m,bands[b]->getMinCountButZero());
  if(m==VERY_LARGE) m=1.;
  return m;
}

double S1S2Bands::getMaxBandContent(){
  double m=0.;
  for(int b=0;b<nBands;b++) m=max(m,bands[b]->getContent());
  return m;
}

double S1S2Bands::getMinBandContent(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) m=min(m,bands[b]->getContent());
  if(m==VERY_LARGE) m=0.;
  return m;
}

double S1S2Bands::getMinBandContentButZero(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) {
    double c=bands[b]->getContent();
    if(c>0.) m=min(m,c);
  }
  if(m==VERY_LARGE) m=1.;
  return m;
}


double S1S2Bands::maxX(){
  double m=0.;
  for(int b=0;b<nBands;b++) m=max(m,bands[b]->maxX());
  return m;
}

double S1S2Bands::minX(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) m=min(m,bands[b]->minX());
  if(m==VERY_LARGE) m=0.;
  return m;
}

double S1S2Bands::maxY(){
  double m=0.;
  for(int b=0;b<nBands;b++) m=max(m,bands[b]->maxY());
  return m;
}

double S1S2Bands::minY(){
  double m=VERY_LARGE;
  for(int b=0;b<nBands;b++) m=min(m,bands[b]->minY());
  if(m==VERY_LARGE) m=0.;
  return m;
}

double S1S2Bands::getArea(){
  double area=0.;
  for(int s=0;s<nSlices;s++) area+= slices[s]->getArea();
  return area;
}

//--------------  XeRun ------------------------------------------

int    XeRun::nBands  = DEFAULT_N_BANDS;
int    XeRun::nSlices = DEFAULT_N_SLICES;

void   XeRun::setNBands(int nb)  { 
  if(nBands>N_MAX_BANDS){
    cout<<"Can't set "<<nb<<" bands; max is "<<N_MAX_BANDS<<endl;
    return;
  }
  nBands=nb;  
  S1S2Display::setNBands(nBands); 
}
void   XeRun::setNSlices(int ns) { nSlices=ns; S1S2Display::setNSlices(ns); }

void   XeRun::traceTheFlags(){
  STATIC_CONST  int TAB=10;
  displayFlag("Changed status"     , TAB, changed              , false);
  displayFlag("Object status"      , TAB, getStatusName()      , isCheckedOK());
  displayFlag("forcedSignalERecoil", TAB, forcedSignalERecoilSpectrum , false);
  displayFlag("SignalTabulated"    , TAB, signalTabulated      , true);
  displayFlag("DataBandsAreBuilt"  , TAB, dataBandsAreBuilt    , true);
  displayFlag("BackgroundComputed" , TAB, backgroundsComputed  ,true);
}

XeRun::~XeRun(){
  for(int t=0;t<N_CUT_DATA_TYPES;t++) deleteWithPointer(s1s2[t]);
  if(referenceBands!=NULL) deleteWithPointer(referenceBands);
  deleteTheBands();

  
    for(unsigned int i=0 ; i < SignalBand.size(); i++){
	for(unsigned int j=0 ; j < SignalBand[i].size(); j++) delete SignalBand[i][j];
    }
    SignalBand.clear();   
}

void XeRun::resetAllFlags(){
  signalTabulated            = false;
  dataBandsAreBuilt          = false;
  backgroundsComputed        = false;
  markUnchecked();
}


string XeRun::getDefaultDirectory(int run){
  return "data/run"+format0I(run,2);
}

string XeRun::getDefaultFileName(int run){
  return getDefaultDirectory(run)+"/bands.root";
}

void XeRun::initialize() {
  Se                            = DEFAULT_SE; 
  Sr                            = DEFAULT_SR; 
  SigmaPMT                      = DEFAULT_SIGMA_PMT;
  LightYield                    = DEFAULT_LIGHT_YIELD_RUN_08;
  smearingPMT                   = DEFAULT_SMEARING_PMT;
  equallyFilledS1Slices         = DEFAULT_EQUALLY_FILLED_S1_SLICES;
  DoMCSignalBand                = DEFAULT_DOMCSIGNALBAND;
  S2OverS1Min                   = DEFAULT_S2_OVER_S1_MIN;
  S2OverS1Max                   = DEFAULT_S2_OVER_S1_MAX;
  SelectionUS1Min               = DEFAULT_US1_MIN ;
  sigmaForReports       = DEFAULT_SIGMA_NUCLEON;
  sigmaFactor                   = 1.;
  SelectionS1Min                = UNDEFINED;
  SelectionS1Max                = UNDEFINED;
  Exposure                      = UNDEFINED;
  FiducialMass                  = UNDEFINED;
  typicalErMin                  = UNDEFINED;
  DarkMatterNormalization       = UNDEFINED;
  AmBeNormalizationToEvents     = UNDEFINED;
  EGammaNormalizationToEvents   = UNDEFINED;
  AnalysisS1Min                 = UNDEFINED; 
  AnalysisS1Max                 = UNDEFINED; 
  firstAnalysisS1Bin            = 0;
  lastAnalysisS1Bin             = N_PE_POINTS -1 ;
  selectionCuts                 = NULL;
  darkMatterCuts                = NULL;
  referenceBands                = NULL;
  ebm                           = NULL;
  nbm                           = NULL;
  leff                          = NULL;
  pValue                        = NULL;
  dataFile			= NULL;
  forcedSignalERecoilSpectrum   = false;
  resetAllFlags();
}

XeRun::XeRun() {initialize();dataType=UNDEFINED_INT;}

XeRun::XeRun(int run, int dt, string file, string sFile , string dm, string em, string ne) 
            : S1S2Object(getTheName(run)){
  initialize(); 
  dataType= dt;
  switch(dataType) {
    case SIMULATED_DATA : 
    case REAL_DATA      : 
         if(file.size()==0) {
           file=getDefaultFileName(run);
           cout<<"... no ROOT file name given for run, assume "<<file<<endl;
         }
         if(file.find('/')==string::npos){
           file=getDefaultDirectory(run)+"/"+file;
           cout<<"... no ROOT directory given, assume "<<file<<endl;
         }
         break;
    case NO_DATA        : 
         if(file.size()>0) {
           cout <<"Run "<<run<<" has no data, ignore file name"<<endl;
         }
         break;
    default  : 
         cout<<"Unknown run data type "<<dt<<endl;
         return;       
  }

  runNumber  = run;
  setAnalysisS1Min(s1sTotCut::defaultS1Min(runNumber)); 
  setAnalysisS1Max(s1sTotCut::defaultS1Max(runNumber));
  if(runNumber == 8) setNSlices(DEFAULT_N_SLICES_RUN_8);

  string ns=getNameSpace();

  SignalERecoilXeSpectrum=new XeSpectrum( ns+ "Signal ERecoil",N_ERS_POINTS
                                , 0., ERS_MAX);
  SignalERecoilSpectrum=SignalERecoilXeSpectrum->getSpectrum();
  attach(SignalERecoilXeSpectrum);

  SignalERecoilXeDist=new TabulatedDist( ns+"Signal ERecoil"
                                       , N_ERS_POINTS, 0., ERS_MAX);
  SignalERecoilDistribution=SignalERecoilXeDist->getSpectrum();
  attach(SignalERecoilXeDist);

  for(int b=FIRST_BACKGROUND; b<=LAST_BACKGROUND; b++){
    string nam=ns+S1S2Data::getTypeName(b);
    
    BackgroundS1XeDist[b]= new TabulatedDist(nam,N_PE_POINTS, 0., PE_MAX);
    BackgroundS1SDistribution[b]=BackgroundS1XeDist[b]->getSpectrum();
    attach(BackgroundS1XeDist[b]);
    
    BackgroundS1XeSpectrum[b]= new XeSpectrum(nam,N_PE_POINTS, 0., PE_MAX);
    BackgroundS1Spectrum[b]=BackgroundS1XeSpectrum[b]->getSpectrum();
    attach(BackgroundS1XeSpectrum[b]);
  }


  for(int t=0;t<N_LEFF_TABULATED;t++) {
    string namt="S1 "+LEff::tabulatedName(t);

    SignalS1XeSpectrum[t]=new XeSpectrum(namt,N_PE_POINTS,0.,PE_MAX);
    SignalS1Spectrum[t]=SignalS1XeSpectrum[t]->getSpectrum();
    attach(SignalS1XeSpectrum[t]);

    SignalS1XeDist[t]=new TabulatedDist(namt, N_PE_POINTS,0.,PE_MAX);  
    SignalS1Distribution[t]=SignalS1XeDist[t]->getSpectrum();
    attach(SignalS1XeDist[t]);
  }

  switch(runNumber){
    case  RUN_08 : setReference("PRL 107, 131302 (2011)");
                   LightYield  = DEFAULT_LIGHT_YIELD_RUN_08; 
                   Exposure=100.9*day;
                   FiducialMass=48.*kg;
                   typicalErMin=8.8;
                   break;
    case RUN_10 : setXenon100Reference("run10:pldata");
                   LightYield  = DEFAULT_LIGHT_YIELD_RUN_10;
                   Exposure=224.56*day;
                   FiducialMass=34.*kg;
                   typicalErMin=5.0;
                   break;
    case RUN_12 : setXenon100Reference("Private communication");
                   LightYield  = DEFAULT_LIGHT_YIELD_RUN_10;
                   Exposure=153.8 *day;
                   FiducialMass=34.*kg;
                   typicalErMin=5.0;
                   break;
     default : 
              cout<<"Don't know the properties of Xenon detector for run "
                  <<runNumber<<endl;
              return;
  }


  setNormalization(); // to default values
  setTarget(getTarget(XE100_MIXTURE));
 flattener= new RunFlattener(this);
  attach(flattener);
  S1S2Display::setCurrentFlattener(flattener);

  for(int t=0;t<N_CUT_DATA_TYPES;t++) s1s2[t]=NULL;
  for(int t=0;t<N_ALL_DATA_TYPES;t++) bands[t]=NULL;
  
  if(withoutData()) return;

  signalFile =new TFile(sFile.c_str(),"OLD");
  if(signalFile==NULL){
    cout<<"Problem: can't access signal file "<<sFile<<endl;
    exit(100);
  }

  dataFile =new TFile(file.c_str(),"OLD");
  if(dataFile==NULL){
    cout<<"Problem: can't access file "<<file<<endl;
    exit(100);
  }

  TGraph2D *gdm2d = (TGraph2D*) dataFile->Get(dm.c_str());
  if(gdm2d==NULL){
    cout<<"Problem: TGraph " + dm  + " does not exist in file " << file <<endl;
    exit(100);
  }
  
  TGraph gdm = from2dto1dGraph(*gdm2d) ;


  s1s2[DM_DATA]=new S1S2Data(Name+" DM",DM_DATA,this,&gdm);
  if(debugLevel>0){
    cout<<"original TGraph has "<<gdm.GetN() <<" events, "
        <<s1s2[DM_DATA]->getName()
        <<" has "<<s1s2[DM_DATA]->getNEvents()<<" events"<<endl;
	s1s2[DM_DATA]->printIt(3);
  } 

  TGraph *gEGamma = (TGraph*) dataFile->Get(em.c_str());
  if(gEGamma==NULL){
    cout<<"Problem: TGraph " + em  + " does not exist in file " << file <<endl;
    exit(100);
  }
  s1s2[E_GAMMA_DATA]=new S1S2Data(Name+" e-Gamma",E_GAMMA_DATA,this,gEGamma);
  if(debugLevel>0){
    cout<<"original TGraph has "<<gEGamma->GetN() <<" events, "
        <<s1s2[E_GAMMA_DATA]->getName()
        <<" has "<<s1s2[E_GAMMA_DATA]->getNEvents()<<" events"<<endl;
  }

  TGraph *gAmBe = (TGraph*) dataFile->Get(ne.c_str());
  if(gAmBe==NULL){
    cout<<"Problem: TGraph " + ne  + " does not exist in file " << file <<endl;
    exit(100);
  }
  s1s2[AM_BE_DATA]=new S1S2Data(Name+" Am-Be",AM_BE_DATA,this,gAmBe);
  if(debugLevel>0){
    cout<<"original TGraph has "<<gAmBe->GetN() <<" events, "
        <<s1s2[AM_BE_DATA]->getName()
        <<" has "<<s1s2[AM_BE_DATA]->getNEvents()<<" events"<<endl;
  }


    //Defining default Leff and Qy
    setLEff(newDefaultLEff()); 
    qy = new Qy();

  // initializing the Signal Handler to Null, so it can be chosen later
  signalHandler = NULL;  

  markAsChanged();
}
    

TGraph XeRun::from2dto1dGraph(TGraph2D gr2d) {

  int nPoints = gr2d.GetN();

  TGraph gr1d;
  Double_t *s1 = gr2d.GetX();
  Double_t *s2 = gr2d.GetY();
   

  for(int p = 0 ; p < nPoints; p++){
	gr1d.SetPoint(p, s1[p], log10(s2[p]/s1[p])); 

  }

  return gr1d;

}
 
void XeRun::setPublishedElectronBackgroundModel(){
  switch(runNumber){
    case RUN_08 : ebm= new PublishedElectronBackgroundRun8(this)  ; break;
    case RUN_10 : ebm= new PublishedElectronBackgroundRun10(this) ; break;
    default :
      cout<<"No PublishedElectronBackground for run "<<runNumber<<endl;
      return;
  }
  cout<<"... assume published electron background model for "<<getName()<<endl;
  setElectronBackgroundModel(ebm);
}


void   XeRun::setElectronBackgroundModel(ElectronBackgroundModel* m) {
  detach(ebm);
  ebm=m;
  attach(ebm);
} 
void   XeRun::setNeutronBackgroundModel(NeutronBackgroundModel* m)   {
  detach(nbm);
  nbm=m;
  attach(nbm);
} 

void   XeRun::setLowestS2OverS1(double r) {S2OverS1Min=r;}
void   XeRun::setHighestS2OverS1(double r){S2OverS1Max=r;}
void   XeRun::setAnalysisS1Min(double s1) {
  AnalysisS1Min=s1;
  firstAnalysisS1Bin=s1/PE_STEP;
}
void   XeRun::setAnalysisS1Max(double s1) {
  AnalysisS1Max=s1;
  lastAnalysisS1Bin=s1/PE_STEP;
}
void   XeRun::setLightYield(double ly)    {LightYield=ly;}
void   XeRun::setPValue(PValue *pv)       {pValue=pv;}

int    XeRun::getNumber()                 {return runNumber;}
int    XeRun::getNSelectionCuts()         {return selectionCuts->getNCuts();}
int    XeRun::getDataType()               {return dataType;}
bool   XeRun::update()                    {fillDataBands(); return true;}
bool   XeRun::withData()                  {return dataType!=NO_DATA;}
bool   XeRun::withoutData()               {return dataType==NO_DATA;}
bool   XeRun::withBands()                 {return withData() && XeStat::isPL();}
bool   XeRun::withoutBands()              {return !withBands();}

double XeRun::getLowestS2OverS1()         {return S2OverS1Min;}
double XeRun::getHighestS2OverS1()        {return S2OverS1Max;}
double XeRun::getSigmaPMT()               {return SigmaPMT;}
double XeRun::getSelectionUS1Min()        {return SelectionUS1Min;}
double XeRun::getSelectionS1Min()         {return SelectionS1Min;}
double XeRun::getSelectionS1Max()         {return SelectionS1Max;}
double XeRun::getAnalysisS1Min()          {return AnalysisS1Min;}
double XeRun::getAnalysisS1Max()          {return AnalysisS1Max;}
double XeRun::getSe()                     {return Se;}
double XeRun::getSr()                     {return Sr;}
double XeRun::getLightYield()             {return LightYield;}
double XeRun::getTypicalErMin()           {return typicalErMin;}
double XeRun::getExposure()               {return Exposure;}
double XeRun::getWimpMass()               {return target->getWimpMass();}
double XeRun::getFiducialMass()           {return FiducialMass;}
double XeRun::ErMax()                     {return target->ErMax();}
double XeRun::getLEffTValue()             {return leff->getTValue();}
double XeRun::getLEffErMin()              {return leff->getErMin();}
string XeRun::getDataTypeName()           {return getDataTypeName(dataType);}

string XeRun::getShortTypeName(int type){
  if(type==DEFAULT_SIGNAL) return formatG(sigmaForReports,8,2);
  return S1S2Data::getShortTypeName(type);
}

PValue*       XeRun::getPValue()          {return pValue;}
Interaction*  XeRun::getInteraction()     {return target->getInteraction();}
Target*       XeRun::getTarget()          {return target;}
Wimp*         XeRun::getWimp()            {return target->getWimp();}
GalaxyModel*  XeRun::getGalaxyModel()     {return target->getGalaxyModel();}
LEff*         XeRun::getLEff()            {return leff;}
XeSetOfCuts*  XeRun::getAllCuts()         {return &allCuts;}
RunFlattener* XeRun::getFlattener()       {return flattener;}
XeSetOfCuts*  XeRun::getDarkMatterCuts()  {return darkMatterCuts;}
XeSetOfSelectionCuts* XeRun::getSelectionCuts() {return selectionCuts;}

LEff*  XeRun::newDefaultLEff()       {return LEff::newDefault(runNumber);}

XeSetOfSelectionCuts* XeRun::newDefaultXeSetOfSelectionCuts() {
  return XeSetOfSelectionCuts::newDefault(this);
}

XeSetOfCuts* XeRun::newDefaultXeSetOfDarkMatterCuts(){
  return XeSetOfCuts::newDefaultXeSetOfDarkMatterCuts(this);
}


ElectronBackgroundModel* XeRun::getElectronBackgroundModel() {return ebm;}
NeutronBackgroundModel*  XeRun::getNeutronBackgroundModel()  {return nbm;}

ElectronBackgroundModel* XeRun::newDefaultElectronBackgroundModel() {
  return ElectronBackgroundModel::newDefault(this);
}

NeutronBackgroundModel* XeRun::newDefaultNeutronBackgroundModel() {
  return NeutronBackgroundModel::newDefault(this);
}

void   XeRun::setLEffTValue(double t){
  if(leff==NULL){
    cout<<"No LEff for run "<<getName()<<", can't set LEff T-value"<<endl;
    return;
  }
  leff->setTValue(t);
  doNotMarkAsChanged(); // because it is tabulated
}  

void   XeRun::setLEffErMin(double e) {
  if(leff==NULL){
    cout<<"No LEff for run "<<getName()<<", can't set LEff ErMin"<<endl;
    return;
  }
  LEffErMin=e;
  leff->setErMin(e);
  markAsChanged();
}

void   XeRun::setWimpMass(double m)  {
  forcedSignalERecoilSpectrum=m<0.;
  target->setWimpMass(forcedSignalERecoilSpectrum? 1000.: m );
  markAsChanged();
  signalTabulated=false;
}

void   XeRun::setWimp(Wimp *w)       {
  target->setWimp(w);
  markAsChanged();
  signalTabulated=false;
}

void   XeRun::setSelectionCutTValue(int i,double t) {
  selectionCuts->setTValue(i,t);
  markAsChanged();
}

double XeRun::getSigmaNucleon() {
  Interaction * inter= getInteraction();
  if(inter==NULL){
    cout<<"No interaction defined in run!"<<endl;
    return UNDEFINED;
  }
  return inter->getSigmaNucleon();
}

double XeRun::dRate(double E)   {
  return target->dRate(E)*Exposure*FiducialMass;
}

double XeRun::nSignalPerCm2(double lt){
  return expectedSignal(lt)/getSigmaNucleon();
}

SelectionCut*       XeRun::getSelectionCutBySequence(int i)  {
  return selectionCuts->getSelectionCutBySequence(i);
}

string XeRun::getDataTypeName(int dt) {
  switch(dt) {
    case NO_DATA        : return " no data";
    case REAL_DATA      : return "";
    case SIMULATED_DATA : return " simulated data";
  }
  cout<<"Invalid data type "<<dt<<endl;
  return "???";
}

string XeRun::getName(){return Name+getDataTypeName();}

string XeRun::getTheName(int r){return "Run"+format0I(r,2);}

double XeRun::flatten(double S1, double S2){
  if(flattener==NULL) {
    cout<<"No flattener defined for "<<getName()<<endl;
    return UNDEFINED;
  }
  return flattener->flatten(S1,S2);
}

double XeRun::expectedSignalInErWindow(double E1,double E2) {
  if(!checkIt()) return 0.;
  double ErMi= E1<0. ? getTypicalErMin() : E1;
  double ErMa= E2<=0.? target->ErMax()   : E2;
  tabulateSignalERecoilSpectrum();
  if(ErMi>=ErMa) return 0.;
  return SignalERecoilXeSpectrum->integral(ErMi,ErMa);
}

void XeRun::tabulateSignalERecoilSpectrum(){
  if(forcedSignalERecoilSpectrum) return;
  trace(ENTRY,"XeRun::tabulateSignalERecoilSpectrum");
  if(!isChanged() && signalTabulated) return;
  stopper *stop=NULL; 
  if(debugLevel>0) stop=new stopper("tabulate SignalERecoilSpectrum");
  for(int e=0;e<N_ERS_POINTS;e++){
    SignalERecoilSpectrum[e]=integrate(RATE,e*ERS_STEP,(e+1)*ERS_STEP); 
  }
  SignalERecoilXeDist->importSpectrum(SignalERecoilSpectrum);
  markAsRecomputed();
  MethodCounter::count("Tabulate ER spectrum");
  if(debugLevel>0)  {stop->print(); delete stop;}
  trace(EXIT,"XeRun::tabulateSignalERecoilSpectrum");
}

void XeRun::printSignalERecoilSpectrum(){
  if(!checkIt()) return;
  SignalERecoilXeSpectrum->printIt(2);
}

void XeRun::setEquallyFilledS1Slices(bool eq)   {equallyFilledS1Slices=eq;}
void XeRun::setDarkMatterCuts(XeSetOfCuts* c)   {
  if(darkMatterCuts==c) return;
  darkMatterCuts=c;
  signalTabulated=false;
}

void XeRun::setSelectionCuts(XeSetOfSelectionCuts* c) {
  if(selectionCuts==c) return;
  selectionCuts=c;
  attach(selectionCuts);
  SelectionS1Min=selectionCuts->getMinS1();
  SelectionS1Max=selectionCuts->getMaxS1();
  markAsChanged();
  signalTabulated=false;
}

void XeRun::setTarget(Target *t) {
  target=t; 
  attach(target);
  markAsChanged(); 
  signalTabulated=false;
}

void XeRun::setLEff(LEff *lef) {
  leff=lef;
  LEffErMin=lef->getErMin();
  attach(leff);
  markAsChanged(); 
  signalTabulated=false;
}

void XeRun::smearPMT(bool smear) {
  smearingPMT=smear;
  markAsChanged();
}

void XeRun::setSigmaNucleon(double s) {getInteraction()->setSigmaNucleon(s);}

void XeRun::simulate(double sigma) {

  dataBandsAreBuilt=false;
  if(withoutData()){
    cout<<"Can't simulate a run without data"<<endl;
    return;
  }
  if(!checkIt() || sigma<0.) {
    cout<<"Can't simulate "<<getName()<<" with sigma="<<sigma<<endl;
    return;
  }

  computeBackgrounds();

  // reset
  if(bands[DM_DATA]==NULL || bands[DM_CUT_DATA]==NULL ||
     s1s2[DM_DATA]==NULL  || s1s2[DM_CUT_DATA]==NULL ) {
    cout<<"Houston, we have a serious problem"<<endl;
    return;
  }
  dataType=SIMULATED_DATA;
  instantiateSignalBandsIfNeeded();
  bands[DM_DATA]->reset();
  bands[DM_CUT_DATA]->reset();
  s1s2[DM_DATA]->reset();
  s1s2[DM_CUT_DATA]->reset();
  PoissonDist poisson(1.);
  UniformDist uniform(minS1(),maxS1());
  
  // fill signal
  if(sigma>0.) {
    for(int b=0;b<nBands;b++){
      double sn = bands[DEFAULT_SIGNAL]->getBand(b)->getContent()
                * sigma/DEFAULT_SIGMA_NUCLEON;
      if(sn>0.) {
        poisson.setMu(sn);
        int ns=(int)(.001+poisson.generate());
        for(int e=0;e<ns;e++){
          double s1=0;
          do {s1=SignalS1XeDist[SIDE_LEFF_TABULATED]->generate();} 
             while(s1>maxS1()||s1<minS1());
          if(!fillSimulatedDataAndBands(b,s1)){
            cout<<"Invalid s1 for simulated signal:"<<s1<<endl;
          }
        }
      }
    }
  }

  // fill gaussian ER  background; no fluctuations for upper band

  for(int b=0;b<nBands;b++){
    double neg= bands[ER_BACKGROUND]->getBand(b)->getContent();
    if(neg>0.) {
      int ne=neg;
      if(b<nBands-1){
        poisson.setMu(neg);
        ne=(int)(.001+poisson.generate());
      }
      if(debugLevel>2) {
        cout<<"Simulating gaussian ER in band ="<<b<<" ne="<<ne<<endl;
      }
      for(int e=0;e<ne;e++){
        double s1=BackgroundS1XeDist[ER_BACKGROUND]->generate();
        if(debugLevel>3) cout<<"    S1="<<s1<<endl;
        if(!fillSimulatedDataAndBands(b,s1)){
          cout<<"Invalid s1 for simulated gaussian ER:"<<s1<<endl;
        }
      } 
    }
  }
  
  // fill Anomalous ER  background; 

  for(int b=0;b<nBands;b++){
    double nea= bands[ER_LEAKAGE]->getBand(b)->getContent();
    if(nea>0.) {
      poisson.setMu(nea);
      int ne=(int)(.001+poisson.generate());
      if(debugLevel>2) {
        cout<<"Simulating anomalous leakage ER in band ="<<b<<" ne="<<ne<<endl;
      }
      for(int e=0;e<ne;e++){
        double s1=BackgroundS1XeDist[ER_LEAKAGE]->generate();
        if(debugLevel>3) cout<<"    S1="<<s1<<endl;
        if(!fillSimulatedDataAndBands(b,s1)){
          cout<<"Invalid s1 for simulated anomalous ER:"<<s1<<endl;
        }
      } 
    }
  }


  // fill neutron background
 
  for(int b=0;b<nBands;b++){
    double nneut= bands[NR_BACKGROUND]->getBand(b)->getContent();
    if(nneut>0.) {
      poisson.setMu(nneut);
      int nn=(int)(.001+poisson.generate());
      if(debugLevel>2) cout<<"Simulating neutron in band "<<b<<" nn="<<nn<<endl;
      for(int e=0;e<nn;e++){
        double s1=BackgroundS1XeDist[NR_BACKGROUND]->generate();
        if(debugLevel>3) cout<<"    S1="<<s1<<endl;
        if(!fillSimulatedDataAndBands(b,s1)){
          cout<<"Invalid s1 for simulated neutron:"<<s1<<endl;
        }
      } 
    }

  }

  fillDataBands();
  s1s2[DM_DATA]->update();
  s1s2[DM_CUT_DATA]->update();
}

bool XeRun::fillSimulatedDataAndBands(int band,double simulatedS1){
  int slice=bands[DM_CUT_DATA]->getS1Bins()->getBin(simulatedS1);
  if(slice<0) return false;
  double s2Mi=bands[DM_CUT_DATA]->getS2overS1LowerEdge(band,slice)*simulatedS1;
  double s2Ma=bands[DM_CUT_DATA]->getS2overS1UpperEdge(band,slice)*simulatedS1;
  double simulatedS2=s2Mi+XeDist::rndm()*(s2Ma-s2Mi);
  s1s2[DM_DATA]->add(simulatedS1,slice,simulatedS2,band);
  bands[DM_DATA]->fill(simulatedS1,slice,band,1.);
  s1s2[DM_CUT_DATA]->add(simulatedS1,slice,simulatedS2,band);
  bands[DM_CUT_DATA]->fill(simulatedS1,slice,band,1.);
  return true;
}

bool XeRun::recomputeEverything(){
  resetAllFlags();
  return checkAll();
}

bool XeRun::checkIt(bool recomputeIfNeeded){
  trace(ENTRY,"XeRun::checkIt()");
  if(traceLevel>0) cout<<"RecomputeIfNeeeded :"<<recomputeIfNeeded<<endl;
  if(dataType==UNDEFINED_INT){
    cout<<"This run isn't defined at all !!!"<<endl;
    return false;
  }
  if(changed) {
    if(recomputeIfNeeded) markUnchecked();
    markAsRecomputed();
  }
  else if(isError()) {
    if(traceLevel>0) cout<<"return false because in error"<<endl;
    return false; 
  }
  if(isCheckedOK() || !recomputeIfNeeded) {
    if(traceLevel>0) {
      cout<<"return true because";
      if(isCheckedOK()) cout<<" checkedOK";
      if(!recomputeIfNeeded) cout<<" no need to recompute";
      cout<<endl;
    }
    return true;
  }

  markOK();
  changed=false;
  bool compatible=true;
  if(!XeStat::isAnalysisDefined(false)){
    cout<<"... no analysis mode defined, assume Profile Likelihood"<<endl;
    XeStat::setAnalysisMode(PL_ANALYSIS);
  }
  if(FiducialMass<=0.) {
    cout<<"Invalid fiducial mass for run "<<getName()<<endl;
    markError();
  }
 
  if(target==NULL){
    cout<<"No target defined for run "<<getName()<<endl;
    markError();
  }
  else if(!target->checkIt())  markError();

  if(leff==NULL) markError();

  if(qy == NULL) markError();

  if(selectionCuts==NULL){
    cout<<"... no set of selection cuts defined for run "<<getName()
        <<", assume default"<<endl;
    setSelectionCuts(newDefaultXeSetOfSelectionCuts());
    //if(selectionCuts==NULL) markError();
  }
  if(!selectionCuts->checkRunCompatibility()) compatible=false;

  if(darkMatterCuts==NULL){
    cout<<"... no dark matter cuts defined for run "<<getName()
        <<", assume default"<<endl;
    setDarkMatterCuts(newDefaultXeSetOfDarkMatterCuts());
    //if(darkMatterCuts==NULL) markError();
  }
  if(!darkMatterCuts->checkRunCompatibility()) compatible=false;

  if(ebm==NULL){
    cout<<"... no electron background model for run "<<getName()
        <<", assume default"<<endl;
    setElectronBackgroundModel(newDefaultElectronBackgroundModel());
    if(ebm==NULL) markError();
  }
  if(!ebm->checkRunCompatibility()) compatible=false;

  if(nbm==NULL){
    cout<<"... no neutron background model for run "<<getName()
        <<", assume default"<<endl;
    setNeutronBackgroundModel(newDefaultNeutronBackgroundModel());
    if(nbm==NULL) markError();
  }
  if(!nbm->checkRunCompatibility()) compatible=false;
 
  if(SelectionS1Min==UNDEFINED || SelectionS1Max==UNDEFINED) {
    cout<<"No S1Min/S1Max cut defined for run "<<getName()<<endl;
    markError();
  }
  
  if(signalHandler==NULL){  
        signalHandler = new MCSignalModel(this);  // just for test to be harmonized FIXME
  }

  if(withData()){
    if(!checkS1LimitsConsistency()) compatible=false;
    changed=false;
    if(isCheckedOK()) {
      fillDataBands();
      changed=false;
      computeBackgrounds();
    }
  }
   
  computeSignal();
  if(isError()){
    cout<<endl
        <<"---------------------------------------------------"<<endl
        <<"  There is an error in the definition of the run   "<<endl
        <<"---------------------------------------------------"<<endl
        <<endl;
  }
  else if(!compatible){
    cout<<endl
        <<"-----------------------------------------------------------"<<endl
        <<"  Warning! Incompatiblity between items. At your own risk  "<<endl
        <<"-----------------------------------------------------------"<<endl
        <<endl;
  }

  trace(EXIT,"XeRun::checkIt()");
  return isCheckedOK();
}

bool XeRun::checkS1LimitsConsistency(){
  trace(ENTRY,"XeRun::checkS1LimitsConsistency()");
  bool ok=true;
  if(AnalysisS1Min<SelectionS1Min){
    cout<<"S1Min for Analysis ("<<AnalysisS1Min
        <<") can't be smaller than S1Min for selection cuts ("
        <<SelectionS1Min<<")"<<endl;
    ok=false;
  }
  else if(AnalysisS1Min>SelectionS1Min){
    cout<<"S1Min for Analysis ("<<AnalysisS1Min
        <<")  larger than S1Min for selection cuts ("
        <<SelectionS1Min<<"); beware!"<<endl;
  }

  if(AnalysisS1Max>SelectionS1Max){
    cout<<"S1Max for Analysis ("<<AnalysisS1Max
        <<") can't be larger than S1Max for selection cuts ("
        <<SelectionS1Max<<")"<<endl;
    ok=false;
  }
  else if(AnalysisS1Max<SelectionS1Max){
    cout<<"S1Max for Analysis ("<<AnalysisS1Max
        <<")  smaller than S1Max for selection cuts ("
        <<SelectionS1Max<<"); beware!"<<endl;
  }
  XeTolerance tol(1.);
  if(!tol.check(SelectionS1Min ,s1s2[DM_DATA]->minS1()
               ,s1s2[E_GAMMA_DATA]->minS1() ,s1s2[AM_BE_DATA]->minS1(),false)){
    ok=false;
  }
  if(!tol.check(SelectionS1Max ,s1s2[DM_DATA]->maxS1()
               ,s1s2[E_GAMMA_DATA]->maxS1() ,s1s2[AM_BE_DATA]->maxS1(),false)){
    ok=false;
  }
  if(!ok) {
    cout<<endl
        <<"Warning ! Inconsistency in S1 cuts"<<endl
        <<"           S1Min S1Max" <<endl
        <<"Selection "<<formatF(SelectionS1Min,6,1)
                      <<formatF(SelectionS1Max,6,1) <<endl
        <<"Analysis  "<<formatF(AnalysisS1Min,6,1)
                      <<formatF(AnalysisS1Max,6,1) <<endl
        <<"DM data   "<<formatF(s1s2[DM_DATA]->minS1(),6,1)
                      <<formatF(s1s2[DM_DATA]->maxS1(),6,1)<<endl
        <<"E/gamma   "<<formatF(s1s2[E_GAMMA_DATA]->minS1(),6,1)
                      <<formatF(s1s2[E_GAMMA_DATA]->maxS1(),6,1)<<endl
        <<"Am-Be     "<<formatF(s1s2[AM_BE_DATA]->minS1(),6,1)
                      <<formatF(s1s2[AM_BE_DATA]->maxS1(),6,1)<<endl;
  }
  trace(EXIT,"XeRun::checkS1LimitsConsistency");
  return ok;
}

void XeRun::setInteraction(Interaction* inter)  {
  target->setInteraction(inter);
}

void XeRun::setGalaxyModel(GalaxyModel* gal) {target->setGalaxyModel(gal);}
bool XeRun::printCandidates(int maxB)    {return printEvents(DM_CUT_DATA,maxB);}
bool XeRun::printEvents(int dt, int maxB) {
  if(!checkIt()) return false;
  s1s2[dt]->printComputedBands(1,maxB);
  cout<<endl;
  return true;
}

bool XeRun::printIt(int level){

  STATIC_CONST  int tab=17;
  STATIC_CONST  int tab2=20;
  if(!checkAll()) return false;
  if(level<1) return true;
  checkTheCuts();
  cout<<endl<<"==================== Data and cuts for "<<getName()
            <<" ===================="<<endl<<endl;
  printReference();
  cout<<endl;
  target->printInteraction();
  cout<<endl
      <<"Target          : "<<getFiducialMass()<<" kg of "<<getName()
      <<" ("<<target->getName()<<")"<<endl;
  target->printReference(tab+1);
  if(level>1) target->printIt();
  cout<<endl
      <<"Exposure        : "<<Exposure<<" days"<<endl
      <<"LEffective      : "<<leff->getName()<<endl
      <<"<NPhotons>      : "<<setprecision(3)
      <<setw(4)<<ErToNPhotons(3.) <<" @ 3KeV  "
      <<setw(7)<<ErToNPhotons(5.) <<" @ 5KeV  "
      <<setw(7)<<ErToNPhotons(10.)<<" @ 10KeV "
      <<endl
      <<"                : "
      <<setw(4)<<ErToNPhotons(30.) <<" @ 30KeV "
      <<setw(7)<<ErToNPhotons(50.) <<" @ 50KeV "
      <<setw(7)<<ErToNPhotons(100.)<<" @ 100KeV"
      <<endl<<endl<<RESET_PRECISION;
  double eMax=ErMax();
  double eMin=max(LEffErMin,getTypicalErMin());
  double n0=expectedSignalInErWindow(0.);
  double n1=expectedSignalInErWindow(LEffErMin);
  double n2=expectedSignalInErWindow(eMin);
  cout<<"Er Spectrum     : LEff_ErMin="<<LEffErMin
      <<setprecision(5)<<", Typical_ErMin="<<eMin
      <<",  ErMax="<<eMax<<" KeV"<<endl<<setprecision(4)
      <<"N events        : Total:"<<n0
      <<", >LEff_ErMin:"<<n1;
  if(eMax>eMin && eMin>LEffErMin) cout<<", >Typical_ErMin:"<<n2;
  cout<<endl<<endl<<RESET_PRECISION
      <<"S1 Spectrum     : flags for computation"<<endl;
  displayFlag("Smear PMT response",tab
                        ,smearingPMT,DEFAULT_SMEARING_PMT);
  displayFlag("Minimum uS1",tab,SelectionUS1Min,DEFAULT_US1_MIN);

  vector<double> remaining;
  fillSelectionCutsBreakdown(remaining);
  cout<<"Selection Cuts  : ";
  selectionCuts->printIt(level+1,&remaining);
  cout<<endl<<"Dark matter cuts: ";
  darkMatterCuts->printIt(level);
  double overall=darkMatterCuts->getOverallAcceptance();
  double nExpected=expectedSignal(0.);
  cout<<"                  Overall acceptance : "
      <<setprecision(4)<<overall<<endl
      <<"Expected events : "<<nExpected<<" after all cuts"
      <<RESET_PRECISION<<endl<<endl;
  if(withoutData()) {
    cout<<endl<<"No data associated to XeRun "<<getName()<<"!"<<endl;
    return true;
  }
  cout<<"Electron backg. : "<<ebm->getName()<<endl;
  ebm->printIt(level);
  cout<<"Neutron  backg. : "<<nbm->getName()<<endl;
  nbm->printIt(level);
  cout<<endl<<"Analysis options: "<<XeStat::getAnalysisModeName()<<endl;
  displayFlag("Minimum S1 for analysis",tab2,AnalysisS1Min
                        ,s1sTotCut::defaultS1Min(runNumber));
  displayFlag("Maximum S1 for analysis",tab2,AnalysisS1Max
                        ,s1sTotCut::defaultS1Max(runNumber));
  displayFlag("Lowest  S2/S1 in bands",tab2,S2OverS1Min
                     ,DEFAULT_S2_OVER_S1_MIN);
  displayFlag("Highest S2/S1 in bands",tab2,S2OverS1Max
                     ,DEFAULT_S2_OVER_S1_MAX);
  displayFlag("S1slices equally filled",tab2
                     ,equallyFilledS1Slices,DEFAULT_EQUALLY_FILLED_S1_SLICES);
  cout<<endl;
  printDataSummary();
  return true;
}
 
void XeRun::printDataSummary(){
  if(withoutData()) return;
  S1S2Data::printSummaryHeader();
  s1s2[DM_DATA]->printSummary("Dark Matter");
  s1s2[E_GAMMA_DATA]->printSummary("E-Gamma    ");
  s1s2[AM_BE_DATA]->printSummary("Am-Be      ");
  cout<<endl<<endl;
  if(withBands()) printBandContent();
}

double XeRun::getValue(int what, double x) {
  switch(what){
    case RATE: return dRate(x);
  }
  cout<<"Invalid quantity :"<<what<<endl;
  return 0.;
}

XeSpectrum* XeRun::getSignalS1XeSpectrum(double leffT) {
  int bin=LEff::getBinAndFraction(leffT).first;
  return getSignalS1XeSpectrum(bin);
}

XeSpectrum* XeRun::getSignalS1XeSpectrum(int tIndex) {
  if(!checkIt()) return NULL;
  if(tIndex<0 ||tIndex>N_LEFF_TABULATED) return NULL; 
  return SignalS1XeSpectrum[tIndex];
}

TabulatedDist* XeRun::getSignalS1XeDist(double leffT) {
  if(!checkIt()) return NULL;
  int bin=LEff::getBinAndFraction(leffT).first;
  return SignalS1XeDist[bin];
}


double* XeRun::getSignalS1Spectrum(int tIndex) {
  if(!checkIt()) return NULL;
  if(tIndex<0 ||tIndex>N_LEFF_TABULATED) return NULL; 
  return SignalS1Spectrum[tIndex];
}
double* XeRun::getSignalS1Spectrum(double leffT) {
  int bin=LEff::getBinAndFraction(leffT).first;
  return getSignalS1Spectrum(bin);
}

double* XeRun::getSignalS1Distribution(double leffT) {
  if(!checkIt()) return NULL;
  int bin=LEff::getBinAndFraction(leffT).first;
  return SignalS1Distribution[bin];
}

double* XeRun::getBackgroundS1Distribution(int type) {
  if(! (checkIt() && S1S2Data::isBackground(type,true) ) ) return NULL;
  computeBackgrounds();
  return BackgroundS1SDistribution[type];
}

double* XeRun::getBackgroundS1Spectrum(int type) {
  if(! (checkIt() && S1S2Data::isBackground(type,true) ) ) return NULL;
  computeBackgrounds();
  return BackgroundS1Spectrum[type];
}


TabulatedDist* XeRun::getBackgroundS1XeDist(int type) {
  if(! (checkIt() && S1S2Data::isBackground(type,true) ) ) return NULL;
  computeBackgrounds();
  return BackgroundS1XeDist[type];
}

XeSpectrum* XeRun::getBackgroundS1XeSpectrum(int type) {
  if(! (checkIt() && S1S2Data::isBackground(type,true) ) ) return NULL;
  computeBackgrounds();
  return BackgroundS1XeSpectrum[type];
}

double XeRun::getTotalContent(int type){
  if(!checkIt()) return UNDEFINED;
  return TotalContent[type];
}

double XeRun::getBandContent(int type, int band){
  if(!checkBand(band)) return UNDEFINED;
  else if(band==ALL)   return getTotalContent(type);
  return getS1S2Bands(type)->getBandContent(band);
}

double  XeRun::expectedSignal(double leffT)  {
  computeSignal();
  pair<int,double> bf=LEff::getBinAndFraction(leffT);
  int bin=bf.first;
  double f=bf.second;
  if(f==0.) return SignalS1Total[bin];
  return SignalS1Total[bin]*(1.-f)+SignalS1Total[bin+1]*f;
}

bool XeRun::computeSignal()   {
  trace(ENTRY,"XeRun::computeSignal()");
  bool b=isChanged() || !signalTabulated;
  if(b){ 
    if(debugLevel>1) {
      cout<<"Need to tabulate "<<getName()<<endl; 
      if(debugLevel>2) printTree();
    }
    tabulateSignalERecoilSpectrum();
    //tabulateSignalS1Spectrum();
    sigmaForReports=getInteraction()->getSigmaNucleon();
    sigmaFactor=1.;
    if(withBands()) computeSignalBands(0.);
    if(debugLevel>2) printTree();
  }
  else {
    if(debugLevel>1) cout<<"No need no tabulate again "<<getName()<<endl; 
  }
  signalTabulated=true;
  trace(EXIT,"XeRun::computeSignal()");
  return b;
}

double XeRun::getSigmaForReports(){
  return sigmaForReports;
}

void   XeRun::setSigmaForReports(double sigma, int unit){
  if(!checkIt()) return;
  switch(unit) {
    case SIGMA_UNIT : 
         sigmaForReports=sigma; 
         sigmaFactor=sigmaForReports/getSigmaNucleon();
         return;
    case EVENT_UNIT : 
         sigmaFactor=sigma/nSignalPerCm2();
         sigmaForReports=getSigmaNucleon()*sigmaFactor; 
         return;
  }
  cout<<"Invalid unit in setSigmaForReports: "<<unit<<endl;
}

void XeRun::forceSignalERecoilSpectrum(double *Er){
  if(Er!=NULL) {
    for(int e=0;e<N_ERS_POINTS;e++) SignalERecoilSpectrum[e]=Er[e];
  }
  setWimpMass(-1.);
  SignalERecoilXeDist->importSpectrum(SignalERecoilSpectrum);
  //tabulateSignalS1Spectrum();                                  // THIS IS ANYWAY NOT NEEDED
  if(withBands()) computeSignalBands(0.);
  signalTabulated=true;
  markAsRecomputed();
}
 
void XeRun::tabulateSignalS1Spectrum() {
  trace(ENTRY,"XeRun::computeSignalS1Spectrum()");
  stopper *stop=NULL;
  if(debugLevel>0) {
    stop=new stopper("Tabulate S1 Spectra, mass="+formatF(getWimpMass(),6,2));
  }
  for(int index=0;index< N_LEFF_TABULATED;index++) {
    tabulateSignalS1Spectrum(index);
  }
  if(debugLevel>0)  {stop->print(); delete stop;}
  trace(EXIT,"XeRun::computeSignalS1Spectrum()");
}

void XeRun::tabulateSignalS1Spectrum(int T_INDEX) {
  MethodCounter::count("Tabulate S1 spectrum");
  // Get Leff value
  double le=LEff::tabulated(T_INDEX);
  leff->setTValue(le);

  // Initialize arrays
  double unsmeared[UNSMEARED_PE_MAX];
  for(int s=0;s<N_PE_POINTS;s++)      SignalS1Spectrum[T_INDEX][s]=0.;
  for(int s=0;s<UNSMEARED_PE_MAX;s++) unsmeared[s]=0.;
  if(debugLevel>1) {
    cout<<"Tabulating S1 spectrum of "<<getName()<<", LEff-t="<<le<<endl;
  }
  leff->markAsRecomputed();
  markAsRecomputed();

  // define minimum and maximum bin for E recoil spectrum
  int bMinQ=SignalERecoilXeDist->quantileBin(1.E-6);
  int bMax=SignalERecoilXeDist->quantileBin(1.-1.E-6);
  int bLEffErMin=LEffErMin/ERS_STEP-.5;
  int bMin=max(bLEffErMin,bMinQ);
  if(bMin>bMax || bMin<0 || bMin>N_ERS_POINTS || bMax<0 || bMax>N_ERS_POINTS){
    cout<<"Strange problem in XeRun::computeSignal()S1Spectrum: bErmin="
        <<bLEffErMin
        <<",bMinQ="<<bMinQ<<",bMin="<<bMin<<", bMax="<<bMax<<endl;
    return;
  }
  if(debugLevel>1) {
    cout<<"  ... bMinQ="<<bMin<<", bLEffErMin"<<bLEffErMin
        <<", bMin="<<bMin<<", bMax="<<bMax<<endl;
  }

  
  //From E recoil gets the smeared cS1 distribution
  //and applies the first acceptance for cuts on S2.
  for(int b=bMin;b<=bMax;b++){
    double ec=(b+.5)*ERS_STEP;
    double npe=ErToNPhotons(ec);   // Get the cS1 from E recoil
    //double npe_s1=ErToNPhotons(ec)/LCE;   // to get the S1 to be fixed /!\ Ale
    if(npe>=SelectionUS1Min) {  // hard cut on unsmeared S1 !!! 
      double r=SignalERecoilSpectrum[b];
      double a=selectionCuts->getAcceptance(npe,false); //These are acceptances on cut on S2
							// which has to be applied at the cS1 level.
      for(int pe=UNSMEARED_PE_MIN;pe<UNSMEARED_PE_MAX;pe++){
        double q=ROOT::Math::poisson_pdf(pe,npe);
        if(pe>npe && q<1.E-8) break;
        unsmeared[pe]+=r*q*a; 
      }
    }
  }



  if(debugLevel>2) {
    printVector(unsmeared,UNSMEARED_PE_MAX,"Unsmeared S1 in p.e.");
  }

  int pmin=(int)(UNSMEARED_PE_MIN/PE_STEP);
  if(smearingPMT){
    // Further gaussian smear to take into account pmt resolution
    // Acceptance of all s1 cuts is applied at this stage
    // /!\ gaussian smear should be applied to S1 not to cS1, to be fixed.  Ale.
    // fill SignalS1Spectrum with the computed spectrum
    for(int s=UNSMEARED_PE_MIN;s<UNSMEARED_PE_MAX;s++) {
      if(unsmeared[s]==0.) continue;
      double sig=SigmaPMT*sqrt((double)s);
      if(sig==0.){
        int b=(int)(s/PE_STEP);
        double a=selectionCuts->getAcceptance((b+.5)*PE_STEP,true);
        SignalS1Spectrum[T_INDEX][b]+=unsmeared[s]*a;  
      }
      else {
        int b0=max(pmin,(int)((s-5.*sig)/PE_STEP));
        b0=max(b0,firstAnalysisS1Bin);
        int b1=min(lastAnalysisS1Bin,(int)((s+5.*sig)/PE_STEP));
        for(int b=b0;b<=b1;b++){
          double s0=(b*PE_STEP-s)/sig;
          double s1=s0+PE_STEP/sig;
          double a=selectionCuts->getAcceptance((b+.5)*PE_STEP,true);
          double g=GaussianDist::inside(s0,s1);  
          SignalS1Spectrum[T_INDEX][b]+=unsmeared[s]*a*g;  
        }
      }
    } 
  }
  else { 
    for(int s=UNSMEARED_PE_MIN;s<UNSMEARED_PE_MAX;s++) {
      int bin=(int)(s/PE_STEP);
      if(bin<N_PE_POINTS) SignalS1Spectrum[T_INDEX][bin]=unsmeared[s];
    }
  }


  //Compute the Total signal
  SignalS1Total[T_INDEX]=0.;
  for(int s=0;s<N_PE_POINTS;s++) {
    double delta=SignalS1Spectrum[T_INDEX][s];
    if(delta<0.){
      cout<<"Problem in "<<getInteraction()->getName()<<", mass="
           <<getWimpMass()<<" LEff-t="<<le<<", bin "<<": S1 <0!"<<endl;
    }  
     else SignalS1Total[T_INDEX]+=delta;
  }
  SignalS1XeDist[T_INDEX]->importSpectrum(SignalS1Spectrum[T_INDEX]);
  if(XeStat::isCutsBased()){
     SignalS1Total[T_INDEX] *=darkMatterCuts->getOverallAcceptance(); 
  }
  if(debugLevel>1) cout<<"S1 tabulated, s1Total="<<SignalS1Total[T_INDEX]<<endl;


}

bool XeRun::checkBand(int band) {
  if(band==ALL || (band>=0 && band<nBands)) return true;
  cout<<"Invalid band number: "<<band<<endl;
  return false;
}

double* XeRun::getAllBackgroundsInBandsS1Spectrum(int band){
  if(!checkBand(band)) return NULL;
  if(band==ALL) return getBackgroundS1Spectrum(ALL_BACKGROUNDS);
  return AllBackgroundsInBandsS1Spectrum[band];
}

double* XeRun::getAllBackgroundsInBandsS1Distribution(int band){
  if(!checkBand(band)) return NULL;
  if(band==ALL) return getBackgroundS1Distribution(ALL_BACKGROUNDS);
  return AllBackgroundsInBandsS1Distribution[band];
}

XeSpectrum* XeRun::getAllBackgroundsInBandsS1XeSpectrum(int band){
  if(!checkBand(band)) return NULL;
  if(band==ALL) return getBackgroundS1XeSpectrum(ALL_BACKGROUNDS);
  return AllBackgroundsInBandsS1XeSpectrum[band];
}

TabulatedDist* XeRun::getAllBackgroundsInBandsS1XeDist(int band){
  if(!checkBand(band)) return NULL;
  if(band==ALL) return getBackgroundS1XeDist(ALL_BACKGROUNDS);
  return AllBackgroundsInBandsS1XeDist[band];
}

XeGraph* XeRun::newGraphOfAllBackgroundsInBandsS1Distribution(int band
          ,int plot){
  if(!checkBand(band)) return NULL;
  TabulatedDist* dist=getAllBackgroundsInBandsS1XeDist(band);
  XeGraph *g= dist->newGraph(LABEL_S1,LABEL_NORMTO1);
  attach(g);
  g->drawIt(plot);
  return g;
}

XeGraph* XeRun::newGraphOfAllBackgroundsInBandsS1Spectrum(int band, int plot){
  if(!checkBand(band)) return NULL;
  XeSpectrum* sp=getAllBackgroundsInBandsS1XeSpectrum(band);
  XeGraph *g= sp->newGraph(LABEL_S1,LABEL_EVT);
  attach(g);
  g->drawIt(plot);
  return g;
}

XeMultiGraph* XeRun::newMultiGraphOfAllBackgroundsInBandsS1Distribution(int pl){
  if(!checkIt()) return NULL;
  XeMultiGraph *mg=new XeMultiGraph("S1 in background",LABEL_S1,LABEL_NORMTO1
                                   ,LABEL_SLICE);
  for(int b=0;b<nBands;b++){
    XeGraph *g=newGraphOfAllBackgroundsInBandsS1Distribution(b);
    g->setBandLegend(b);
    mg->add(g,b);
  }
  mg->setRainbowColors(ASYMMETRIC);
  attach(mg);
  mg->drawIt(pl);
  return mg;
}

XeMultiGraph* XeRun::newMultiGraphOfAllBackgroundsInBandsS1Spectrum(int plot){
  if(!checkIt()) return NULL;
  XeMultiGraph *mg=new XeMultiGraph("S1 in background",LABEL_S1,LABEL_EVT
                                   ,LABEL_SLICE);
  for(int b=0;b<nBands;b++){
    XeGraph *g=newGraphOfAllBackgroundsInBandsS1Spectrum(b);
    g->setBandLegend(b);
    mg->add(g,b);
  }
  mg->setRainbowColors(ASYMMETRIC);
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

XeMultiGraph* XeRun::newMultiGraphOfSignalS1Spectrum(int step,int plot){

  if(!checkIt()) return NULL;
  XeMultiGraph *mg=new XeMultiGraph("S1 dependence on LEff",LABEL_S1,LABEL_EVT
                                   ,LABEL_LEFF_TVALUE);
  for(int i=0;i<N_LEFF_TABULATED;i+=step){
    double le=LEff::tabulated(i);
    XeGraph *g=newGraphOfSignalS1Spectrum(i);
    g->setTValueLegend(le);
    mg->add(g,le);
  }
  mg->setRainbowColors(SYMMETRIC);
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

XeGraph* XeRun::newGraphOfSignalS1Spectrum(int plot, int t) {
  if(!checkIt()) return NULL;
  XeStyle* style=getInteraction()->getStyleFromFF();
  string leg=getInteraction()->getLegendFromFF();
  XeGraph *g=SignalS1XeSpectrum[t]->newGraph(LABEL_S1,LABEL_EVT_PE,style,leg);
  g->setDefaultYScale(LOG);
  attach(g);
  g->drawIt(plot);
  return g;
} 

TH1F* XeRun::newHistogramOfSignalS1Spectrum(int plot,int t) {
  if(!checkIt()) return NULL;
  string nam="S1 spectrum of signal, "+LEff::tabulatedName(t);
  TH1F* h= SignalS1XeSpectrum[t]->newHistogram(nam);
  drawHist(h,plot,LABEL_S1);
  return h;
} 

TH2F* XeRun::new2DHistogramOfSignalS1Spectrum(XeRange* mr,int plot,int t) {
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  TH2F* h = NULL;
  for(int i=0;i<n;i++){
    double mas=mr->getValue(i);
    if(debugLevel>0) cout<<"Filling for mass "<<mas<<endl;
    setWimpMass(mas);
    computeSignal();
    if(i==0) h=SignalS1XeSpectrum[t]->newHistogram(n);
    SignalS1XeSpectrum[t]->fillHistogram(h,i);
  }
  drawHist(h,plot,LABEL_S1);
  return h;
}


XeGraph* XeRun::newGraphOfSignalERecoilSpectrum(int plot) {
  if(!checkIt()) return NULL;
  XeStyle* style=getInteraction()->getStyleFromFF();
  string leg=getInteraction()->getLegendFromFF();
  XeGraph *g=SignalERecoilXeSpectrum->newGraph(LABEL_ER,LABEL_EVT_KEV,style,leg);
  g->setDefaultYScale(LOG);
  attach(g);
  g->drawIt(plot);
  return g;
} 

XeGraph* XeRun::newGraphOfSignalERecoilDistribution(int plot) {
  if(!checkIt()) return NULL;
  XeStyle* style=getInteraction()->getStyleFromFF();
  string leg=getInteraction()->getLegendFromFF();
  XeGraph *g=SignalERecoilXeDist->newGraph(LABEL_ER,LABEL_NORMTO1,style,leg);
  g->setDefaultYScale(LOG);
  attach(g);
  g->drawIt(plot);
  return g;
} 

TH1F* XeRun::newHistogramOfSignalERecoilSpectrum(int plot) {
  if(!checkIt()) return NULL;
  return SignalERecoilXeSpectrum->newHistogram("E Recoil",plot);
}

TH1F* XeRun::newHistogramOfSignalERecoilDisribution(int plot) {
  if(!checkIt()) return NULL;
  return SignalERecoilXeDist->newHistogram("E Recoil, normalized to 1",plot);
}

XeMultiGraph* XeRun::newMultiGraphOfSignalERecoilSpectrum(XeRange* range
        , int plot) {
  if(!checkIt()) return NULL;
  if(range==NULL) range=XeRange::getDefaultMassRange();
  int n=range->getNPoints();
  XeMultiGraph * mg=new XeMultiGraph("Er for various masses"
                   ,LABEL_ER,LABEL_EVT_KEV,LABEL_MASS);
  for(int i=0;i<n;i++){
    double mass=range->getValue(i);
    setWimpMass(mass);
    mg->add(newGraphOfSignalERecoilSpectrum(),mass);
  }
  mg->setReferenceGraph(range->getReferenceIndex());
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

XeMultiGraph* XeRun::newMultiGraphOfSignalERecoilDistribution(XeRange* range
               , int plot) {
  if(!checkIt()) return NULL;
  if(range==NULL) range=XeRange::getDefaultMassRange();
  int n=range->getNPoints();
  XeMultiGraph * mg=new XeMultiGraph("Er for various masses"
                   ,LABEL_ER,LABEL_NORMTO1,LABEL_MASS);
  for(int i=0;i<n;i++){
    double mass=range->getValue(i);
    setWimpMass(mass);
    mg->add(newGraphOfSignalERecoilDistribution(),mass);
  }
  mg->setReferenceGraph(range->getReferenceIndex());
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

TH2F* XeRun::new2DHistogramOfSignalERecoilSpectrum(XeRange* range,int plot) {
  if(!checkIt()) return NULL;
  if(range==NULL) range=XeRange::getDefaultMassRange();
  int n=range->getNPoints();
  TH2F* h = NULL;
  for(int i=0;i<n;i++){
    double mas=range->getValue(i);
    if(debugLevel>0) cout<<"Filling for mass "<<mas<<endl;
    setWimpMass(mas);
    computeSignal();
    if(i==0) h=SignalERecoilXeSpectrum->newHistogram(n);
    SignalERecoilXeSpectrum->fillHistogram(h,i);
  }
  drawHist(h,plot);
  return h;
}

TH2F* XeRun::new2DHistogramOfSignalERecoilDistribution(XeRange* range,int plot){
  if(!checkIt()) return NULL;
  if(range==NULL) range=XeRange::getDefaultMassRange();
  int n=range->getNPoints();
  TH2F* h =NULL;
  for(int i=0;i<n;i++){
    double mas=range->getValue(i);
    if(debugLevel>0) cout<<"Filling for mass "<<mas<<endl;
    setWimpMass(mas);
    computeSignal();
    if(i==0) {
      h=SignalERecoilXeDist->newHistogram(n);
    }
    SignalERecoilXeDist->fillHistogram(h,i);
  }
  drawHist(h,plot);
  return h;
}

double* XeRun::getSignalERecoilDistribution() {
  if(!checkIt()) return NULL;
  return SignalERecoilDistribution;
}

double* XeRun::getSignalERecoilSpectrum() {
  if(!checkIt()) return NULL;
  return SignalERecoilSpectrum;
}

XeSpectrum* XeRun::getSignalERecoilXeSpectrum() {
  checkIt();
  return SignalERecoilXeSpectrum;
}

TabulatedDist* XeRun::getSignalERecoilXeDist() {
  checkIt();
  return SignalERecoilXeDist;
}


double XeRun::ErToNPhotons(double Er){
   return leff->getLEff(Er)*Sr/Se*LightYield*Er;
}

double XeRun::NPhotonsToEr(double nP){
  double eff=leff->getLEff(10.);
  double Er=0.;
  for(int i=0;i<5;i++){
    Er=nP/eff*Se/Sr/LightYield;
    eff=leff->getLEff(Er);
  }
  return Er;
}

XeGraph* XeRun::newGraphOfExpectedSignalInErWindow(XeRange *mr,int plot){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  XeStyle* style=getInteraction()->getStyleFromFF();
  string leg=getInteraction()->getLegendFromFF();
  XeGraph* gr=new XeGraph("Expected events in "+getName(),n
                          ,LABEL_MASS,LABEL_EVT,style,leg);
  for(int i=0;i<n;i++){
    double mas=mr->getValue(i);
    setWimpMass(mas);
    if(debugLevel>0){
      cout<<"Mass for newGraphOfExpectedSignal set to "<<mas<<"..."<<flush;
    } 
    double e=expectedSignalInErWindow();
    gr->setPoint(i,mas,e);
    if(debugLevel>0){
      cout<<"    result is "<<e<<endl;
    } 
  }
  gr->setDefaultXScale(LOG);
  gr->setDefaultYScale(LOG);
  attach(gr);
  gr->drawIt(plot);
  return gr;
}

XeGraph* XeRun::newGraphOfExpectedSignal(XeRange *mr,double tl,int plot){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  XeStyle* style=getInteraction()->getStyleFromFF();
  string leg=getInteraction()->getLegendFromFF();
  XeGraph* gr=new XeGraph("Expected Signal in "+getName(),n
             ,LABEL_MASS,LABEL_EVT_AFTER_SEL,style,leg);
  for(int i=0;i<n;i++){
    double mas=mr->getValue(i);
    setWimpMass(mas);
    if(debugLevel>0){
      cout<<"Mass for newGraphOfExpectedSignal set to "<<mas
          <<endl;
    } 
    double e=expectedSignal(tl);
    gr->setPoint(i,mas,e);
  }
  gr->setDefaultYScale(LOG);
  gr->setDefaultXScale(LOG);
  attach(gr);
  gr->drawIt(plot);
  return gr;
}

XeGraph* XeRun::newGraphOfLEff(ErRange *er, double t,int plot){
  if(!checkIt()) return NULL;
  return leff->newGraphOfLEff(er,t,plot); 
}

XeGraph* XeRun::newGraphOfCurrentLEff(ErRange *er,int plot){
  if(!checkIt()) return NULL;
  return leff->newGraphOfLEff(er,plot); 
}

XeGraph* XeRun::newGraphOfSelectionCutsAcceptance(S1Range *s1, bool smeared
                , int plot){
  if(!checkIt()) return NULL;
  if(s1==NULL) s1=XeRange::getDefaultS1Range();
  int n=s1->getNPoints();
  XeGraph* g=new XeGraph("Selection cuts acceptance",n,LABEL_S1);
  for(int i=0;i<n;i++){
    double p=s1->getValue(i);
    double a=selectionCuts->getAcceptance(p,smeared);
    g->setPoint(i,p,a);
  }
  attach(g);
  g->drawIt(plot);
  return g;
}

XeMultiGraph* XeRun::newMultiGraphOfSignalS1Distribution(XeRange *mr,int plot
                 , int t){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  XeMultiGraph* mg=new XeMultiGraph("Normalized S1",LABEL_S1,LABEL_NORMTO1);
  for(int i=0;i<n;i++){
    double mas=mr->getValue(i);
    setWimpMass(mas);
    XeGraph*g=newGraphOfSignalS1Distribution(t);
    g->setMassLegend(mas);
    mg->add(g);
  }
  mg->setReferenceGraph(mr->getReferenceIndex());
  mg->setAsymmetricRainbowColors();
  attach(mg);
  mg->drawIt(plot);
  return mg;
}


XeGraph* XeRun::newGraphOfSignalS1Distribution(int plot,int lt){
  computeSignal();
  computeSignal();
  XeGraph *g=SignalS1XeDist[lt]->newGraph(LABEL_S1,LABEL_NORMTO1,NULL,"Signal");
  attach(g);
  g->drawIt(plot);
  return g;
}

TH1F* XeRun::newHistogramOfSignalS1Distribution(int plot,int lt){
  computeSignal();
  return SignalS1XeDist[lt]->newHistogram("Signal S1 distribution",plot);
}

TH2F* XeRun::new2DHistogramOfSignalS1Distribution(XeRange* mr,int plot, int lt){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  TH2F* h = NULL;
  for(int i=0;i<n;i++){
    double mas=mr->getValue(i);
    if(debugLevel>0) cout<<"Filling for mass "<<mas<<endl;
    setWimpMass(mas);
    computeSignal();
    if(i==0) {
      h=SignalS1XeDist[lt]->newHistogram(n);
    }
    SignalS1XeDist[lt]->fillHistogram(h,i);
  }
  drawHist(h,plot,LABEL_MASS,LABEL_S1);
  return h;
}

XeGraph* XeRun::newGraphOfBackgroundS1Distribution(int type, int plot){
  if(! (checkIt() && S1S2Data::isBackground(type,true) ) ) return NULL;
  computeBackgrounds();
  XeGraph *g=BackgroundS1XeDist[type]->newGraph(LABEL_S1,LABEL_NORMTO1,NULL
                      ,S1S2Data::getTypeName(type)+"S1 distribution");
  g->setDefaultYScale(LOG);
  attach(g);
  g->drawIt(plot);
  return g;
}

TH1F* XeRun::newHistogramOfBackgroundS1Distribution(int type, int plot){
  if(! (checkIt() && S1S2Data::isBackground(type,true) ) ) return NULL;
  if(!checkIt()) return NULL;
  computeBackgrounds();
  return BackgroundS1XeDist[type]->newHistogram(S1S2Data::getTypeName(type)
                                       +" S1 distribution",plot);
}

XeMultiGraph* XeRun::newMultiGraphOfS1Distributions(int plot){
  if(!checkIt()) return NULL;
  vector<XeGraph*> v;
  v.push_back(newGraphOfSignalS1Distribution());
  for(int type=FIRST_BACKGROUND;type<=LAST_BACKGROUND;type++){
    v.push_back(newGraphOfBackgroundS1Distribution(type));
  }
  XeMultiGraph* mg=new XeMultiGraph("S1 distribution",v,NULL
                                    ,LABEL_S1,LABEL_NORMTO1);
  mg->setRainbowColors(ASYMMETRIC);
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

XeSpectrum* XeRun::newSpectrumOfSignalAndBackgroundS1(int t){
  return new XeSpectrum(*getSignalS1XeSpectrum(t)
                       +*getBackgroundS1XeSpectrum(ALL_BACKGROUNDS));
}

TabulatedDist* XeRun::newDistributionOfSignalAndBackgroundS1(int t){
  XeSpectrum    *spect= newSpectrumOfSignalAndBackgroundS1(t);
  TabulatedDist *dist = new TabulatedDist(spect);
  delete spect;
  return dist;
}

XeGraph* XeRun::newGraphOfSignalAndBackgroundS1Spectrum(int plot, int t){
  XeSpectrum *spect = newSpectrumOfSignalAndBackgroundS1(t);
  XeGraph    *graph = spect->newGraph(LABEL_S1,LABEL_EVT);
  graph->drawIt(plot);
  delete spect;
  return graph;
}

XeGraph* XeRun::newGraphOfSignalAndBackgroundS1Distribution(int plot, int t){
  TabulatedDist *dist = newDistributionOfSignalAndBackgroundS1(t);
  XeGraph       *graph= dist->newGraph(LABEL_S1,LABEL_EVT);
  graph->drawIt(plot);
  delete dist;
  return graph;
}


XeGraph* XeRun::newGraphOfTheoreticalUpperSigma(XeRange *mr,double cl
        , bool afterCuts, int plot){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  XeStyle* style=getInteraction()->getStyleFromFF();
  string leg=getInteraction()->getLegendFromFF();
  XeGraph* gr=new XeGraph(
                     "\\sigma excluded at"+formatF(100.*cl,4)+" \% C.L"
                     ,n,LABEL_MASS,LABEL_SIGMA,style,leg);
  for(int i=0;i<n;i++){
    double mas=mr->getValue(i);
    setWimpMass(mas);
    double upper=getTheoreticalUpperSigma(cl,afterCuts);
    gr->setPoint(i,mas,upper);
    if(getDebugLevel()>0) {
      cout<<"Mass="<<mas<<", upper sigma="<<upper<<endl;
    }
  }
  gr->setDefaultXScale(LOG);
  gr->setDefaultYScale(LOG);
  attach(gr);
  gr->drawIt(plot);
  return gr;
}

XeMultiGraph* XeRun::newMultiGraphOfExpectedSignalByLEff(XeRange *mr
                                  , tValueRange* tr, int plot){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  if(tr==NULL) tr=XeRange::getDefaultTValueRange();
  int nm=mr->getNPoints();
  int nt=tr->getNPoints();
  vector<XeGraph*> vg(nt);
  vector<double> vz;
  for(int t=0;t<nt;t++){
    double lef=tr->getValue(t);
    vz.push_back(lef);
    string leg="LEff t-value="+formatF(lef,4,1);
    vg[t]=new XeGraph("Expected signal, "+leg,nm,LABEL_MASS,LABEL_EVT,NULL,leg);
    vg[t]->setDefaultXScale(LOG);
  }
  for(int m=0;m<nm;m++){ 
    double mass=mr->getValue(m);
    setWimpMass(mass);
    computeSignal();
    for(int t=0;t<nt;t++){
      double tv=tr->getValue(t);
      double s=expectedSignal(tv);
      vg[t]->setPoint(m,mass,s);
    }
  } 
  XeMultiGraph *mg=new XeMultiGraph("Expected number of events after cuts",vg
                                   ,&vz,LABEL_MASS,LABEL_EVT,LABEL_LEFF_TVALUE);
  mg->setRainbowColors(SYMMETRIC);
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

XeMultiGraph* XeRun::newMultiGraphOfExpectedSignalByVEsc( XeRange *mr
                                , VEscRange *vr, int plot){
  if(!checkIt()) return NULL;
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  if(vr==NULL) vr=XeRange::getDefaultVEscRange();
  int nv=vr->getNPoints();
  
  XeMultiGraph *mg=new XeMultiGraph("Expected number of events after cuts"
                                   ,LABEL_MASS,LABEL_EVT,LABEL_VESC);
  for(int v=0;v<nv;v++){ 
    double vesc=vr->getValue(v);
    getGalaxyModel()->setV_ESC(vesc);
    XeGraph* g=newGraphOfExpectedSignal(mr);
    string leg=LABEL_VESC+" = "+formatF(vesc,6,0);
    cout<<"Computed graph for "<<leg<<endl;
    g->setLegend(leg);
    mg->add(g,vesc);
  }
  mg->setRainbowColors(SYMMETRIC);
  attach(mg);
  mg->drawIt(plot);
  return mg;
}

double  XeRun::getTheoreticalUpperSigma(double cl,bool afterCuts){
  if(!checkIt()) return UNDEFINED;
  double evMax=-log(1-cl);
  double s=getInteraction()->getSigmaNucleon();
  double e=afterCuts ? expectedSignal(0.)
                     : expectedSignalInErWindow();
  if(e<0.) {
    cout<<"Negative number of events for "<<getInteraction()->getName()
        <<", mass="<<getWimpMass()<<", after cuts="<<afterCuts<<endl;
    return 0.;
  }
  return evMax*s/e;
}

void XeRun::printSelectionCutsBreakdown(double wmass,double sigma,double leffT){
  if(!checkIt()) return;
  if(wmass>0.) setWimpMass(wmass);
  if(sigma>0.) setSigmaNucleon(sigma);
  if(!checkIt()) return;
  vector<double> remaining;
  fillSelectionCutsBreakdown(remaining,leffT);
  cout<<leftJustify("No cut",30)<<" : "<<remaining[0]<<endl;
  int nc=selectionCuts->getNCuts();
  for(int c=0;c<nc;c++){
    string cn=selectionCuts->getCutNameBySequence(c+1);
    cout<<leftJustify("After "+cn,30)<<" : "<<remaining[c]<<endl;
  }
}


void XeRun::fillSelectionCutsBreakdown(vector<double> &remaining,double leffT){
  remaining.clear();
  selectionCuts->disableAll();
  //tabulateSignalS1Spectrum();
  remaining.push_back(expectedSignal(leffT));
  int nc=selectionCuts->getNCuts();
  for(int c=0;c<nc;c++){
    selectionCuts->enableFirst(1+c);
    //tabulateSignalS1Spectrum();
    remaining.push_back(expectedSignal(leffT));
  }
}

double XeRun::expectedNumberOfNeutrons(double peMin,double peMax){
  if(nbm==NULL) {
    cout<<"No neutron background model defined for "<<getName()<<endl;
    return UNDEFINED;
  }
  return nbm->expected(peMin,peMax);
}

Target* XeRun::getTarget(int a){
  Target *t = new Target(XENON,a);
  if(a==XE100_MIXTURE) {
    Target *natural  = new Target(XENON,NATURAL); 
    Target *depleted = new Target(XENON,DEPLETED); 
    double fractionOfNaturalXenon=0.869;
    double fractionOfDepletedXenon=1.-fractionOfNaturalXenon;
    t->add(natural,fractionOfNaturalXenon);
    t->add(depleted,fractionOfDepletedXenon);
  }
  return t;
}

void XeRun::setAutomaticLimits(bool var) {
  for(int t=0;t<N_CUT_DATA_TYPES;t++) s1s2[t]->setAutomaticLimits(var);
}

S1S2Bands* XeRun::getReferenceBands(){return referenceBands;}

S1S2Bands* XeRun::getS1S2Bands(int type){
  trace(ENTRY,"XeRun::getS1S2Bands");
  if(!checkIt(false)) return NULL;
  trace(EXIT,"XeRun::getS1S2Bands");
  if(type>=0 && type<N_ALL_DATA_TYPES) return bands[type];
  cout<<"Can't return bands of "<<S1S2Data::getTypeName(type)<<endl;
  return NULL;
}

S1S2Data* XeRun::getS1S2Data(int type){
  if(!checkIt(false)) return NULL;
  if(type>=0 && type<N_CUT_DATA_TYPES) {
    if(S1S2Data::isAfterCuts(type)) checkTheCuts();
    return s1s2[type];
  }
  cout<<"Can't return S1S2 data set of "<<S1S2Data::getTypeName(type)<<endl;
  return NULL;
}

TTree* XeRun::newTree(int type){
  if(!checkIt()) return NULL;
  S1S2Data* s12=getS1S2Data(type);
  if(s12==NULL) return NULL;
  return s12->newTree(referenceBands);
}

TH2F* XeRun::new2DHistogram(int type, int plot){
  S1S2Bands* bds=getS1S2Bands(type);
  if(bds==NULL) return NULL;
  return bds->new2DHistogram(plot);
}

void XeRun::setNormalization(double egamma , double ambe , double darkMatter) {
  DarkMatterNormalization=darkMatter;
  AmBeNormalizationToEvents=ambe;
  EGammaNormalizationToEvents=egamma;
}

double XeRun::defaultDarkMatterNormalization()  {return 1.;}

double XeRun::defaultAmBeNormalizationToEvents(){
  switch(runNumber){
    //This is used for the cutbased analysis. For PL will be overwritten. 
    case RUN_08 : return .119;  // Corresponds to 1.8 events total
    case RUN_10 : return .39 ;  // between 4 and 20
  }
  cout<<endl<<"=> Don't know how to normalize AmBe of run "<<runNumber<<endl;
  return 0.;
}

double XeRun::defaultEGammaNormalizationToEvents(){
  switch(runNumber){
    //This is used for the cutbased analysis. For PL will be overwritten. 
    case RUN_08 : return 1.14+.54;  // Gaussian and Anomalous
    case RUN_10 : return .79;       // all together
  }
  cout<<endl<<"=> Don't know how to normalize E-Gamma of run "<<runNumber<<endl;
  return 0.;
}

void XeRun::checkTheCuts(){
  if(withData() && s1s2[DM_CUT_DATA]==NULL) {
    cout<<"You forgot to apply the cuts and normalize, doing it for you"<<endl;
  }
}


double XeRun::maxS1(int type)  {return s1s2[type]->maxS1();}
double XeRun::maxS2(int type)  {return s1s2[type]->maxS2();}
double XeRun::minS1(int type)  {return s1s2[type]->minS1();}
double XeRun::minS2(int type)  {return s1s2[type]->minS2();}
double XeRun::maxX(int type)   {return s1s2[type]->maxX();}
double XeRun::maxY(int type)   {return s1s2[type]->maxY();}
double XeRun::minX(int type)   {return s1s2[type]->minX();}
double XeRun::minY(int type)   {return s1s2[type]->minY();}

double XeRun::maxS1(){
  double m=VERY_SMALL;
  for(int t=0;t<=DM_DATA;t++) m=max(m,maxS1(t));
  return m;
}

double XeRun::maxS2(){
  double m=VERY_SMALL;
  for(int t=0;t<=DM_DATA;t++) m=max(m,maxS2(t));
  return m;
}

double XeRun::maxX(){
  double m=VERY_SMALL;
  for(int t=0;t<=DM_DATA;t++) m=max(m,maxX(t));
  return m;
}

double XeRun::maxY(){
  double m=VERY_SMALL;
  for(int t=0;t<=DM_DATA;t++) m=max(m,maxY(t));
  return m;
}

double XeRun::minS1(){
  double m=VERY_LARGE;
  for(int t=0;t<=DM_DATA;t++) m=min(m,minS1(t));
  return m;
}

double XeRun::minS2(){
  double m=VERY_LARGE;
  for(int t=0;t<=DM_DATA;t++) m=min(m,minS2(t));
  return m;
}

double XeRun::minX(){
  double m=VERY_LARGE;
  for(int t=0;t<=DM_DATA;t++) m=min(m,minX(t));
  return m;
}

double XeRun::minY(){
  double m=VERY_LARGE;
  for(int t=0;t<=DM_DATA;t++) m=min(m,minY(t));
  return m;
}

void XeRun::draw(){
  if(debugLevel>0) cout<<"Drawing "<<getName()<<endl;
  s1s2[DM_DATA]->draw();
  s1s2[AM_BE_DATA]->draw();
  s1s2[E_GAMMA_DATA]->draw();
}

void XeRun::drawCuts(){getAllCuts()->draw();}

void XeRun::deleteTheBands(int first, int last){
  for(int t=first;t<=last;t++) deleteTheBand(t);
}

void XeRun::deleteTheBand(int t) {deleteWithPointer(bands[t]);};

void XeRun::instantiateSignalBandsIfNeeded(){
  trace(ENTRY,"XeRun::instantiateSignalBandsIfNeeded");
  if(bands[DEFAULT_SIGNAL]!=NULL) return;
  if(bands[DM_DATA]==NULL) {
   cout<<"Inconsistency in XeRun::instantiateSignalBandsIfNeeded()"<<endl;
  } 
  else bands[DEFAULT_SIGNAL]=
     new  S1S2Bands(getName()+" Expected signal",bands[DM_DATA]);
  trace(EXIT,"XeRun::instantiateSignalBandsIfNeeded");
}

void XeRun::computeSignalBands(double leffT){
  trace(ENTRY,"XeRun::computeSignalBands");
  deleteTheBand(DEFAULT_SIGNAL);
  instantiateSignalBandsIfNeeded();
  pair<int,double> bf=LEff::getBinAndFraction(leffT);
  int bin=bf.first;
  double f=bf.second;
  for(int i=0;i<N_PE_POINTS;i++) {
    double w=SignalS1Spectrum[bin][i];
    if(f>0.) w= w*(1.-f)+SignalS1Spectrum[bin+1][i]*f;
    bands[DEFAULT_SIGNAL]->dispatchS1InEqualBands(i*PE_STEP,w);
  }
  deleteTheBand(DEFAULT_SIGNAL_INTEGRATED);
  bands[DEFAULT_SIGNAL_INTEGRATED]=new S1S2Bands(
      getName()+" Total Signal integrated", bands[DEFAULT_SIGNAL]);
  bands[DEFAULT_SIGNAL_INTEGRATED]->cumulate(bands[DEFAULT_SIGNAL]);
  trace(EXIT,"XeRun::computeSignalBands");
}

void XeRun::tabulateSignalLeffQy(){

    if(signalHandler==NULL) { cout <<"XeRun::tabulateSignalLeffQy()  - FATAL ERROR : Signal Model not defined, exit" << endl;  exit(100);}

    if(debugLevel > 0) cout << "XeRun::tabulateSignalLeffQy() - DEBUG : entering tabulation" << endl;

    //empty stored value if any
    for(unsigned int i=0 ; i < SignalBand.size(); i++){
	for(unsigned int j=0 ; j < SignalBand[i].size(); j++) delete SignalBand[i][j];
    }
    SignalBand.clear();   

    //Get minimum and max Leff and Qy
    double Leff_min = leff->getMinimumTval();
    double Qy_min   = qy->getMinimumTval();
    double Leff_max = leff->getMaximumTval();
    double Qy_max   = qy->getMaximumTval();

    double LEff_step_size = leff->getStepSize();
    double Qy_step_size   = qy->getStepSize();

    if(debugLevel > 1) cout << "Leff_min " << Leff_min <<" Qy_min " << Qy_min  << " Leff_max " << Leff_max
			     << " Qy_max " << Qy_max <<" LEff_step_size " << LEff_step_size << " Qy_step_size " << Qy_step_size << endl;

    // Looping on all Histograms of Leff and Qy, all possible combinations, push them in bands
    // and store them into SignalBand vector. 
    for(double Leff = Leff_min; Leff <= Leff_max; Leff += LEff_step_size) {

	vector <S1S2Bands*> tempBands;

	    // this is just because of rounding in Double which wont get to zero perfectly but remains negative -1E-9
	    if(abs(Leff) < 0.0001) Leff = 0.; 

	for(double Qy = Qy_min; Qy <= Qy_max; Qy += Qy_step_size){

	    if(debugLevel > 0)  cout << "XeRun::tabulateSignalLeffQy() - DEBUG : Leff " << Leff << "  Qy " << Qy << endl;

	    // this is just because of rounding in Double which wont get to zero perfectly but remains negative -1E-9
	    if(abs(Qy) < 0.0001) Qy = 0.; 

	    TString Name = "SignalBands LEff "+ formatF(Leff)+ " Qy " +formatF(Qy);
	    S1S2Bands *signalBandsTemp = new S1S2Bands( Name.Data(), signalHandler->computeSignalBands(Leff, Qy),true);

  	   if(debugLevel > 1.) cout << "XeRun::tabulateSignalLeffQy() INFO : " << signalBandsTemp->getName() << " N signal events rescaled to " << signalBandsTemp->getTotalContent() << "  *  " << signalHandler->getXsecMultiplier() << "  =  " << signalBandsTemp->getTotalContent() * signalHandler->getXsecMultiplier() << endl;

	   tempBands.push_back(signalBandsTemp);
	}
 	SignalBand.push_back(tempBands);
    }
    
    cout << "XeRun::tabulateSignalLeffQy() - DEBUG : END of tabulation" << endl;

    //now that the signal bands are filled we can set the XsecMultiplier to the necessary value to normalize to 1 event
    if(signalHandler->isNormOne()) {
	    double normTemp = getSignalBand(0.,0.)->getTotalContent();
	    signalHandler->setXsecMultiplier(1. / normTemp );
    }


}

unsigned int XeRun::getTabulatedIndex(double Leff){
	//finding the zero
	unsigned int Leff_zero_bin = leff->getMaximumTval() / leff->getStepSize();

	unsigned int bin_number = fabs(Leff / leff->getStepSize());

	if(Leff >= 0.) bin_number += Leff_zero_bin;
	else 	       bin_number =  Leff_zero_bin - bin_number;

        return bin_number;
}

S1S2Bands XeRun::interpolateLeffQy(double Leff_val, double Qy_val){
   // The required Leff and Qy combination lies in a 2D-grid of Qy and Leff discrete values for wich histo are computed
   // four histograms has to be loaded and linearly interpolated. Note that linear interpolation is just
   // a simple approximation and wont give the exact result, the number of steps has to be properly calibrated
   // otherwise the result wont represent the effect.

   double Qy_step_size   = qy->getStepSize() ;
   double LEff_step_size = leff->getStepSize();
   int Leff_N_steps = Leff_val / LEff_step_size; 
   int Qy_N_steps   = Qy_val / Qy_step_size ; 

   double Leff_residual = Leff_val / LEff_step_size - Leff_N_steps; 
   double Qy_residual =   Qy_val   / Qy_step_size   - Qy_N_steps; 
   
   //NOTE: *_N_steps will underestimate in case Leff or Qy are positive and overestimate in case are negative,
   // so in case of negative Qy and Leff one has to subtract to obtain the "upper" edge.
   double Qy_plus = 0;
   if(Qy_N_steps >= 0) Qy_plus = (Qy_N_steps+1)*Qy_step_size;
   else 	       Qy_plus = (Qy_N_steps-1)*Qy_step_size;
   int Leff_plus = 0;
   if(Leff_N_steps >= 0) Leff_plus = (Leff_N_steps+1)*LEff_step_size;
   else 	       Leff_plus = (Leff_N_steps-1)*LEff_step_size;

   //protect against minimum/maximum value
   if(Leff_val <= leff->getMinimumTval()) Leff_plus = Leff_val;
   if(Leff_val >= leff->getMaximumTval()) Leff_plus = Leff_val;
   if(Qy_val <= qy->getMinimumTval())  Qy_plus = Qy_val;
   if(Qy_val >= qy->getMaximumTval())  Qy_plus = Qy_val;

/*   if(debugLevel > 1) cout << "S1S2PL::interpolateLeffQy() - DEBUG : required Leff " << Leff  << "  Qy  " << Qy  
			    << " retriving Leff_N_steps " << Leff_N_steps  << " Qy_N_steps " << Qy_N_steps << " Leff_residual " 
                            << Leff_residual << "  Qy_residual  "  << Qy_residual <<"\n\t\t HISTO: Leff = " << Leff_N_steps*LEff_step_size << " Leff+1 =" << Leff_N_steps*LEff_step_size + LEff_step_size
				<< "  Qy= " << Qy_N_steps * Qy_step_size << "  Qy+1= " << Qy_N_steps * Qy_step_size + Qy_step_size << endl;
*/
   // retrieve the three histograms lowQy_lowLeff, HighQy_lowLeff, lowQy_HigLeff
   S1S2Bands lowQy_lowLeff ( "TEST", SignalBand[getTabulatedIndex(Leff_val)][ getTabulatedIndex(Qy_val)],true);
   S1S2Bands highQy_lowLeff( "TEST", SignalBand[getTabulatedIndex(Leff_val)][ getTabulatedIndex(Qy_plus)],true);
   S1S2Bands lowQy_highLeff( "TEST", SignalBand[getTabulatedIndex(Leff_plus)][ getTabulatedIndex(Qy_val)],true);

  //Normalize the histos: interpolate separately along the two dimensions (Qy, Leff) and then average the result.
  lowQy_lowLeff.normalize((1.- Qy_residual)/2. + (1. - Leff_residual)/2.);
  highQy_lowLeff.normalize( Qy_residual/2. );
  lowQy_highLeff.normalize( Leff_residual/2.);

  //Add histo and return final combination
  lowQy_lowLeff.add(&highQy_lowLeff);
  lowQy_lowLeff.add(&lowQy_highLeff);


  return lowQy_lowLeff;
  
}

S1S2Bands * XeRun::getSignalBand(double t_Leff, double t_Qy){

    return SignalBand[getTabulatedIndex(t_Leff)][ getTabulatedIndex(t_Qy)];

}


void XeRun::computeBackgrounds(){
  trace(ENTRY,"XeRun::computeBackgrounds");
  if(backgroundsComputed && !changed) return;

  // all background bands are computed *BEFORE* acceptance!"

  deleteTheBand(ER_BACKGROUND);
  deleteTheBand(ER_LEAKAGE);
  deleteTheBand(NR_BACKGROUND);
  deleteTheBand(ALL_BACKGROUNDS);
  bands[ER_BACKGROUND] = ebm->computeBands();
  bands[ER_LEAKAGE]    = ebm->computeAnomalousBands();
  bands[NR_BACKGROUND] = nbm->computeBands();

  bands[ALL_BACKGROUNDS]=new S1S2Bands(getName()+" Total Background"
                                       ,bands[ER_BACKGROUND]);
 
  for(int b=0;b<nBands;b++){
    for(int s=0;s<nSlices;s++){
      double bkg=bands[ER_LEAKAGE]->getContent(s,b)
                +bands[ER_BACKGROUND]->getContent(s,b)
                +bands[NR_BACKGROUND]->getContent(s,b);
      bands[ALL_BACKGROUNDS]->fill(s,b,bkg);
    }
  } 
 
  deleteTheBand(ALL_BACKGROUNDS_INTEGRATED);
  bands[ALL_BACKGROUNDS_INTEGRATED]=new S1S2Bands(
             getName()+" TotalBackground integtrated",bands[ALL_BACKGROUNDS]);
  bands[ALL_BACKGROUNDS_INTEGRATED]->cumulate(bands[ALL_BACKGROUNDS]);


// fill S1 distributions of the various background
// Smaller grain than slice distributions in S1S2bands;
  for(int p=0;p<N_PE_POINTS;p++){
    double s1=p*PE_STEP;
    double s2=(p+1)*PE_STEP;
    if(s1<SelectionS1Max && s2>SelectionS1Min){

      BackgroundS1Spectrum[NR_BACKGROUND][p] = nbm->expected(s1,s2);

      BackgroundS1Spectrum[ER_BACKGROUND][p]
         = bands[ER_BACKGROUND]->getS1Density(s1);

      BackgroundS1Spectrum[ER_LEAKAGE][p] 
         = bands[ER_LEAKAGE]->getS1Density(s1);

      BackgroundS1Spectrum[ALL_BACKGROUNDS][p] 
         = BackgroundS1Spectrum[ER_LEAKAGE][p] 
         + BackgroundS1Spectrum[ER_BACKGROUND][p]
         + BackgroundS1Spectrum[NR_BACKGROUND][p]
         ;
    }
  } 

  for(int type=FIRST_BACKGROUND;type<=LAST_BACKGROUND;type++) {
    BackgroundS1XeSpectrum[type]->importSpectrum(BackgroundS1Spectrum[type]);
    BackgroundS1XeDist[type]->importSpectrum(BackgroundS1Spectrum[type]);
    TotalContent[type]=bands[type]->getTotalContent();
  }
  
// Recompute S1 spectrum for all combined backgrounds
// If needed, delete existing bands since nBands might have changed

  int size=AllBackgroundsInBandsS1XeSpectrum.size();
  for(int b=0;b<size;b++){
    deleteWithPointer(AllBackgroundsInBandsS1XeSpectrum[b]);
    deleteWithPointer(AllBackgroundsInBandsS1XeDist[b]);
  }
  AllBackgroundsInBandsS1XeSpectrum.resize(nBands);
  AllBackgroundsInBandsS1Spectrum.resize(nBands);
  AllBackgroundsInBandsS1XeDist.resize(nBands);
  AllBackgroundsInBandsS1Distribution.resize(nBands);

// the assumption is that the S1 spectrum is the same in all bands 

  for(int b=0;b<nBands;b++){
    string nam="All backgrounds in band "+formatI(b,3);
    AllBackgroundsInBandsS1XeSpectrum[b]=new XeSpectrum(nam,N_PE_POINTS
                                                       , 0., PE_MAX);
    AllBackgroundsInBandsS1Spectrum[b]=AllBackgroundsInBandsS1XeSpectrum[b]
                                      ->getSpectrum();
    double el=bands[ER_LEAKAGE]->getBandContent(b);
    double er=bands[ER_BACKGROUND]->getBandContent(b);
    double nr=bands[NR_BACKGROUND]->getBandContent(b);
    for(int p=0;p<N_PE_POINTS;p++){
      AllBackgroundsInBandsS1Spectrum[b][p]=
                      el*BackgroundS1SDistribution[ER_LEAKAGE][p]
                    + er*BackgroundS1SDistribution[ER_BACKGROUND][p]
                    + nr*BackgroundS1SDistribution[NR_BACKGROUND][p] 
                    ;
    } 

    AllBackgroundsInBandsS1XeDist[b]=new TabulatedDist(
                                         AllBackgroundsInBandsS1XeSpectrum[b]);
    AllBackgroundsInBandsS1Distribution[b]=AllBackgroundsInBandsS1XeDist[b]
                                        ->getSpectrum();
  }  

  backgroundsComputed=true;
  trace(EXIT,"XeRun::computeBackgrounds");
}

void XeRun::printBandContent(double leffT){
  trace(ENTRY,"XeRun::printBandContent");
  const int nTypes=9;
  const int types[nTypes]={ ER_BACKGROUND
                          , ER_LEAKAGE
                          , NR_BACKGROUND
                          , ALL_BACKGROUNDS
                          , ALL_BACKGROUNDS_INTEGRATED
                          , DM_CUT_DATA
                          , DM_CUT_DATA_INTEGRATED
                          , DEFAULT_SIGNAL
                          , DEFAULT_SIGNAL_INTEGRATED
                          } ;
  if(withoutBands()) return;
  computeSignalBands(leffT);
  double total[nTypes];
  string line(80,'-');
  cout<<line<<endl<<"Band  ";
  for(int ty=0;ty<nTypes;ty++) {
    total[ty]=0;
    int t=types[ty];
    if(t==DM_CUT_DATA) cout<<" ";
    cout<<rightJustify(getShortTypeName(t),8);
  }
  cout<<endl<<line<<endl;
  for(int b=0;b<nBands;b++){
    cout<<setw(3)<<b<<" ";
    for(int t=0;t<nTypes;t++) {
      int ty=types[t];
      S1S2Bands *band=bands[ty];
      if(band==NULL) {
        cout<<"*** Internal problem, ty="<<getShortTypeName(ty)<<endl;
        return;
      }
      double c=band->getBandContent(b);
      total[t]+=c;
      if(t==DM_CUT_DATA) cout<<"  ";
      cout<<S1S2Data::formatContent(ty,c,8);
    }
    cout<<endl;
  }
  cout<<line<<endl<<"Sum  ";
  for(int t=0;t<nTypes;t++) {
    if(t==DM_CUT_DATA) cout<<"  ";
    cout<<S1S2Data::formatContent(types[t],total[t],8,true);
  }
  cout<<endl<<endl;
  trace(EXIT,"XeRun::printBandContent");
}

void XeRun::printBandCLsCounting(int nb){
  cout<<endl<<"Band Signal Backg. Limit Overall"<<endl;
  for(int b=0;b<nb;b++){
    int s=(int)bands[DM_CUT_DATA]->getCumulatedBandContent(b);
    double bkg=bands[ALL_BACKGROUNDS]->getCumulatedBandContent(b); 
    double limit=PoissonCI::upperLimit(CLS_UP,s,bkg,DEFAULT_CL);
    printf("%4d %4d %7.2f %7.2f %7.2f \n",b,s,bkg,limit,limit*nBands/(b+1));
  }
  cout<<endl;
}

void XeRun::fillDataBands(){
  trace(ENTRY,"XeRun::fillDataBands");
  if(dataBandsAreBuilt) return;
  if(!XeStat::isPL()) return;

  // apply all cuts and normalize
  allCuts.clear();
  allCuts.addCuts(selectionCuts);
  allCuts.addCuts(darkMatterCuts);

  s1s2[DM_CUT_DATA]= new S1S2Data(s1s2[DM_DATA],&allCuts); 
  s1s2[E_GAMMA_CUT_DATA]=new S1S2Data(s1s2[E_GAMMA_DATA],&allCuts); 
  s1s2[AM_BE_CUT_DATA]=new S1S2Data(s1s2[AM_BE_DATA],&allCuts);

  double s=DarkMatterNormalization>0.? DarkMatterNormalization 
                                     : defaultDarkMatterNormalization();
  s1s2[DM_DATA]->normalize(s);
  s1s2[DM_CUT_DATA]->normalize(s);

  double a=AmBeNormalizationToEvents>0.? AmBeNormalizationToEvents
                                 : defaultAmBeNormalizationToEvents();
  s1s2[AM_BE_CUT_DATA]->normalizeToEvents(a);

  double n=s1s2[AM_BE_CUT_DATA]->getNormalization();
  s1s2[AM_BE_DATA]->normalize(n);

  double g=EGammaNormalizationToEvents>0.? EGammaNormalizationToEvents
                                 : defaultEGammaNormalizationToEvents();
  s1s2[E_GAMMA_CUT_DATA]->normalizeToEvents(g);
  n=s1s2[E_GAMMA_CUT_DATA]->getNormalization();
  s1s2[E_GAMMA_DATA]->normalize(n);

  // define reference bands
  if(equallyFilledS1Slices){
    referenceBands=new S1S2Bands(getNameSpace()+"Reference bands",this
                                ,nBands,s1s2[AM_BE_CUT_DATA],nSlices);
    referenceBands->extendS1Range(AnalysisS1Min,AnalysisS1Max);
  }
  else {
    referenceBands=new S1S2Bands(getNameSpace()+"Reference bands",this
                                ,nBands,nSlices,AnalysisS1Min,AnalysisS1Max);
  }

  if(DoMCSignalBand){
	cout <<"XeRun::fillDataBands() - INFO : Computing reference band from MC... " << endl;
        referenceBands->compute(signalHandler->getHisto(), false);
  	cout << "XeRun::fillDataBands() INFO : " << referenceBands->getName() << " N signal events rescaled to " << referenceBands->getTotalContent() * signalHandler->getXsecMultiplier() << endl;
	
	cout <<"XeRun::fillDataBands() - INFO : Tabulating signal bands from MC... " << endl;
	tabulateSignalLeffQy();	
   }
  else{ 
	// referenceBands->compute(s1s2[E_GAMMA_CUT_DATA]);
	 referenceBands->compute(s1s2[AM_BE_CUT_DATA]);
         s1s2[E_GAMMA_CUT_DATA]->fillAndMarkBandsAndSlices(referenceBands,false);
         s1s2[AM_BE_CUT_DATA]->fillAndMarkBandsAndSlices(referenceBands,false);
   }



  for(int t=0;t<N_CUT_DATA_TYPES;t++){
    deleteTheBand(t);
    bands[t]=new S1S2Bands(getNameSpace()+S1S2Data::getTypeName(t)
                          ,referenceBands);
    s1s2[t]->fillAndMarkBandsAndSlices(bands[t],true);
    TotalContent[t]=bands[t]->getTotalContent();
  }
  int t=DM_CUT_DATA_INTEGRATED;
  deleteTheBand(t);
  bands[t]=new S1S2Bands(getName()+" DM" ,bands[DM_CUT_DATA]);
  bands[t]->cumulate(bands[DM_CUT_DATA]);
  dataBandsAreBuilt=true;
  trace(EXIT,"XeRun::fillDataBands");
}
//----------------------------------------------------------//
//-----------------Generate Data according to BKG model -----///

void XeRun::generateData(double seed, double sigma){

  //you have to make sure that bkg is computed
  if(!checkAll()) {cout <<"XeRun::generateData()  - ERROR: something wron in setup, exit." <<endl;
		   exit(100);}

  //checks
  if(SignalBand.empty()) { cout << "\nXeRun::generateData - FATAL ERROR - Signal Bands have NOT been tabulated. \n\t.....Exiting....\n" << endl; 
			   exit(100);
  }

  //delete if is istanziated already 
  deleteTheBand(DM_SIMULATED_DATA);  
 
  TRandom3 rand;
  rand.SetSeed(seed);
  bands[DM_SIMULATED_DATA] = new S1S2Bands(getNameSpace()+ " DM Simulated Data sigma "+ formatF(sigma,2,2) ,  referenceBands); 
  bands[DM_SIMULATED_DATA]->printBands();

  //loop over bkg bands, for each fill data in generatedData bands
  for(int b=0; b < nBands ; b++){
	double bkgBandContent = bands[ALL_BACKGROUNDS]->getBandContent(b);
	double signalBandContent = sigma * signalHandler->getXsecMultiplier() *getSignalBand(0.,0.)->getBandContent(b); // signal for nominal Leff&Qy
	int GeneratedBandContent = rand.Poisson(bkgBandContent + signalBandContent);

	cout << "band " << b << " bkgBandContent " << bkgBandContent << " GeneratedBandContent " << GeneratedBandContent << endl;

	TH1F h_bkg    = bands[ALL_BACKGROUNDS]->getHistogramOfS1Distro(b);	
	TH1F h_signal = getSignalBand(0,0)->getHistogramOfS1Distro(b);
	h_signal.Scale(sigma * signalHandler->getXsecMultiplier());
	h_bkg.Add(&h_signal);
	// generate S1
	for(int j=0; j<GeneratedBandContent; j++){
	    double s1Point      = h_bkg.GetRandom();
	    int relativeS1Slice = h_bkg.FindBin(s1Point);
	    bands[DM_SIMULATED_DATA]->fill(s1Point,relativeS1Slice-1, b);
	}
  }
  cout << "AFTER GENERATION" << endl;
  bands[DM_SIMULATED_DATA]->printBands();

}


void XeRun::generateAsimovData(double sigma){

  //you have to make sure that bkg is computed
  if(!checkAll()) {cout <<"XeRun::generateAsimovData()  - ERROR: something wrong in setup, exit." <<endl;
		   exit(100);}

  //delete if is istanziated already 
  deleteTheBand(ASIMOV_DATA);  

  //checks
  if(SignalBand.empty()) { cout << "\nXeRun::generateAsimovData - FATAL ERROR - Signal Bands have NOT been tabulated. \n\t.....Exiting....\n" << endl; 
			   exit(100);
  }

  // copy the signal bands scaled to nominal value to the ASIMOV_BAND
  bands[ASIMOV_DATA] = new S1S2Bands(getNameSpace()+ " Asimov Data for sigma "+formatF(sigma,2,2) , getSignalBand(0.,0.)  , true); 

  // normalize to sigma times the nominal value	
  double norm = sigma * signalHandler->getXsecMultiplier();
  bands[ASIMOV_DATA]->normalize(norm);

  // add all the bkg
  bands[ASIMOV_DATA]->add(bands[ALL_BACKGROUNDS]);

}





//-------------------  Background Models ------------------------------

XeBackgroundModel::~XeBackgroundModel() {}
XeBackgroundModel::XeBackgroundModel() : XeObject(), RunComponent() {
  requestedAnalysisMode=NO_ANALYSIS;
}

XeBackgroundModel::XeBackgroundModel(string nam, XeRun* r,int req) : XeObject() 
                   , RunComponent(){
  setName(nam);
  setRun(r);
  dayKg = run->getFiducialMass()*run->getExposure();
  expected=0.;
  requestedAnalysisMode=req;
  checkRunCompatibility();
}

bool XeBackgroundModel::checkIt(bool){
  bool ok=XeStat::checkAnalysisMode(Name,requestedAnalysisMode);
  markFromCheck(ok);
  return ok;
}

double XeBackgroundModel::getExpected()             {return expected;}
void   XeBackgroundModel::setExpected(double e)     {expected=e;}
void   XeBackgroundModel::printExpected(){
  string tab(18,' ');
  cout<<tab<<"- Total background         : "<<formatF(expected,8,2)<<endl;
}


NeutronBackgroundModel::~NeutronBackgroundModel() {}
NeutronBackgroundModel::NeutronBackgroundModel() : XeBackgroundModel() {}
NeutronBackgroundModel::NeutronBackgroundModel(string n, XeRun* r, int req) 
  : XeBackgroundModel("Neutron background model "+n,r,req) {
  modelType=NR_BACKGROUND_MODEL;
}

NeutronBackgroundModel* NeutronBackgroundModel::newDefault(XeRun* r){
  switch(r->getNumber()) {
    case RUN_08 : return new NeutronBackgroundModelRun8(r);
    case RUN_10 : return new NeutronBackgroundModelRun10(r);
  }
  cout<<"No default neutron background model for run "<<r->getName()<<endl;
  return NULL;
}

double NeutronBackgroundModel::expected(double peMin,double peMax){
  return expectedBelow(peMax)-expectedBelow(peMin);
}
bool NeutronBackgroundModel::printIt(int level){
  if(level<1) return true;
  printExpected();
  return true;
}

NeutronBackgroundModelRun10::~NeutronBackgroundModelRun10(){}
NeutronBackgroundModelRun10::NeutronBackgroundModelRun10() 
  : NeutronBackgroundModel() {}
NeutronBackgroundModelRun10::NeutronBackgroundModelRun10(XeRun* r)
  : NeutronBackgroundModel("Standard for run 10",r,PL_ANALYSIS) {}
bool NeutronBackgroundModelRun10::isRunCompatible(int r) {return r==10;}

S1S2Bands* NeutronBackgroundModelRun10::computeBands(){
  trace(ENTRY,"NeutronBackgroundModelRun10::computeBands");
  S1S2Bands* abands=run->getS1S2Bands(AM_BE_CUT_DATA);
  if(abands==NULL) {
    cout<<"Error: Can't compute neutron bands!"<<endl;
    return NULL;
  }
  int nSlices=abands->getNSlices();
  int nBands=abands->getNBands();
  S1S2Bands* nbands=new S1S2Bands(
    run->getNameSpace()+S1S2Data::getTypeName(NR_BACKGROUND),abands,false);
  double s1Mi = abands->getS1LowerEdge();
  double s1Ma = abands->getS1UpperEdge();
  double total= expected(s1Mi,s1Ma);
  for(int s=0;s<nSlices;s++){
    s1Mi = abands->getS1LowerEdge(s);
    s1Ma = abands->getS1UpperEdge(s);
    //Now that we use MC the neutron expectation is not equal for all bands
    //IMPLEMENTATION: distribution along S2 from AmBe, while along S1 is computed
    double tot_AmBe_in_S1Slice = abands->getSlice(s)->getContent() ;
    double expected_alongS1 = expected(s1Mi,s1Ma);
    for(int b=0;b<nBands;b++) { 
		double w = expected_alongS1 * abands->getContent(s,b) / tot_AmBe_in_S1Slice;
		nbands->fill(s,b,w);
    }
  }
  setExpected(total);
  trace(EXIT,"NeutronBackgroundModelRun10::computeBands");
  return nbands;
}

double NeutronBackgroundModelRun10::expectedBelow(double pe){
  int index=(int)pe;
  if(index<0) return 0.;
  else if(index>=N_PE_NR_BACKGROUND-1) {
    return _integrated[N_PE_NR_BACKGROUND-1];
  }
  double x=pe-index;
  return dayKg*((1.-x)*_integrated[index] + x *_integrated[index+1]) ;
}

NeutronBackgroundModelRun8::~NeutronBackgroundModelRun8(){}
NeutronBackgroundModelRun8::NeutronBackgroundModelRun8() 
  : NeutronBackgroundModel() {}

NeutronBackgroundModelRun8::NeutronBackgroundModelRun8(XeRun* r)
 :  NeutronBackgroundModel("Standard for run 8",r,PL_ANALYSIS) {}

bool NeutronBackgroundModelRun8::isRunCompatible(int r) {return r==8;}

S1S2Bands* NeutronBackgroundModelRun8::computeBands(){
 ////////////////////////////////////////////////////////////
 //Neutron contribution
 // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat3
 // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:run8ubp:nrbg
 // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:code_implementation	
 ////////////////////////////////////////////////////////////

  trace(ENTRY,"NeutronBackgroundModelRun8::computeBands");

  //load the bands produced from the AmBe data graph after applying the 4 < S1 < 30 cut.
  S1S2Bands* abands=run->getS1S2Bands(AM_BE_CUT_DATA);

  if(abands==NULL) {
    cout<<"Error : Can't compute neutron bands!"<<endl;
    return NULL;
  }

  // producig a copy of abands (NOTE: the last option is TRUE, this means copy band content)
  S1S2Bands* nbands=new S1S2Bands(run->getNameSpace()+S1S2Data::getTypeName(NR_BACKGROUND),abands,true);

  double exp_neutron 	   = 0.31; // Number of predicted NR, error is 0.22  should be a NP!!! Ale. FIXME
  double nevents_4_30_pe   = 56529.;  // Number of events in AmBe dataset afer cuts.

  // sanity check with run8 notes
  if(nbands->getTotalContent() != nevents_4_30_pe) 
	{
	  cout << "\nFatal Error in NeutronBackgroundModelRun8::computeBands() - normalization of AmBe doesn't match, maybe you are using the wrong input file. \nQuitting..." << endl;
	  cout << nevents_4_30_pe << "   " << nbands->getTotalContent() << endl;
	  exit(100);
	}
	

  nbands->normalize(exp_neutron / nevents_4_30_pe);

  trace(EXIT,"NeutronBackgroundModelRun8::computeBands");
  return nbands;
}

double NeutronBackgroundModelRun8::expectedBelow(double pe){
  return run->getS1S2Bands(NR_BACKGROUND)->getContentBelowS1(pe);
}


// ------------------ Electron Background model ---------------------

ElectronBackgroundModel::~ElectronBackgroundModel() {}
ElectronBackgroundModel::ElectronBackgroundModel() : XeBackgroundModel() {}
ElectronBackgroundModel::ElectronBackgroundModel(string n, XeRun* r, int req) 
  : XeBackgroundModel("E-recoil background model "+n,r,req) {

  trace(ENTRY,"ElectronBackgroundModel::ElectronBackgroundModel");
  modelType=ELECTRON_BACKGROUND_MODEL;
  gaussianExpected=0.;
  anomalousExpected=0.;
  anomalousSlope=0.;
  ERGaussianFlatInS1=DEFAULT_ER_GAUSSIAN_FLAT_IN_S1;
  switch(runNumber){
    case RUN_08 : normCo60=1./43.3;
                  break;
    case RUN_10 : normCo60=0.027;  
                  break;
    default : 
      cout<<"Don't know (yet) how to compute the ER background for run"
          <<runNumber<<endl;
  }
  trace(EXIT,"ElectronBackgroundModel::ElectronBackgroundModel");
}

ElectronBackgroundModel* ElectronBackgroundModel::newDefault(XeRun* r){
  int rn=r->getNumber();
  switch(rn) {
    case RUN_08 : return new PublishedElectronBackgroundRun8(r) ;  
    case RUN_10 : return new SimpleERBackground(r);
  }
  cout<<"No default electron background model for run "<<r->getName()<<endl;
  return NULL;
}

bool ElectronBackgroundModel::printIt(int level){
  if(level<1) return true;
  string tab(18,' ');
  displayFlag("Gaussian ER flat in S1",17,ERGaussianFlatInS1
             ,DEFAULT_ER_GAUSSIAN_FLAT_IN_S1);
  cout<<tab<<"- Average of flattened     : " <<formatF(meanCo60,8,3)<<endl
      <<tab<<"- Sigma of flattened       : " <<formatF(sigmaCo60,8,3)<<endl
      <<tab<<"- Normalization to DM      : " <<formatF(normCo60,8,4)<<endl
      <<tab<<"- Gaussian  leakage        : " <<formatF(gaussianExpected,8,2) 
      <<endl
      <<tab<<"- Anomalous leakage        : " <<formatF(anomalousExpected,8,2)
      <<endl
      <<tab<<"- S1 slope of anomalous    : " <<formatF(anomalousSlope,8,3)
      <<endl;
  printExpected();
  return true;
 } 

void ElectronBackgroundModel::setERGaussianFlatInS1(bool flat){
  ERGaussianFlatInS1=flat;
}

bool ElectronBackgroundModel::isERGaussianFlatInS1(){return ERGaussianFlatInS1;}

double ElectronBackgroundModel::getExpectedGaussian()   {
  return gaussianExpected;
}

void   ElectronBackgroundModel::setExpectedGaussian(double e) {
  gaussianExpected=e;
  expected=gaussianExpected+anomalousExpected;
}

double ElectronBackgroundModel::getExpectedAnomalous() {
  return anomalousExpected;
}

double ElectronBackgroundModel::getAnomalousSlope() {
  return anomalousSlope;
}
void ElectronBackgroundModel::setAnomalousSlope(double slope) {
  anomalousSlope=slope;
}

void   ElectronBackgroundModel::setExpectedAnomalous(double e) {
  anomalousExpected=e;
  expected=gaussianExpected+anomalousExpected;
}

SimplisticERBackground::SimplisticERBackground() : ElectronBackgroundModel() {}
SimplisticERBackground::SimplisticERBackground(XeRun* r) :
  ElectronBackgroundModel("Simplistic, fixed gaussian leakage",r,PL_ANALYSIS) {
  doSpreadLeakageOverBands   = DEFAULT_SPREAD_LEAKAGE_OVER_BANDS;
  switch(runNumber){
    case RUN_08 : meanCo60=0.022242;   
                  sigmaCo60=0.135846;
                  break;
    case RUN_10 : meanCo60=0.000873;
                  sigmaCo60=0.136333;
                  break;
    default : 
      cout<<"Don't know (yet) how to compute the ER background for run"
          <<runNumber<<endl;
  }
}

SimplisticERBackground::~SimplisticERBackground(){}
void SimplisticERBackground::spreadLeakageOverBands(bool spread) {
  doSpreadLeakageOverBands=spread;
}
bool SimplisticERBackground::isLeakageSpreadOverBands() {
  return doSpreadLeakageOverBands;
}
bool SimplisticERBackground::printIt(int level){
  if(level<1) return true;
  ElectronBackgroundModel::printIt(level);
  displayFlag("Spread over bands",17,doSpreadLeakageOverBands
                        ,DEFAULT_SPREAD_LEAKAGE_OVER_BANDS);
  return true;
}

S1S2Bands* SimplisticERBackground::computeBands(){
  trace(ENTRY,"SimplisticERBackground::computeBands");
  S1S2Bands* ebands=run->getS1S2Bands(E_GAMMA_CUT_DATA);
  if(ebands==NULL) {
    cout<<"Error: can't compute gaussian leakage!"<<endl;
    return NULL;
  }
  int nSlices=ebands->getNSlices();
  int nBands=ebands->getNBands();
  S1S2Bands* gbands=new S1S2Bands(
             run->getNameSpace()+S1S2Data::getTypeName(ER_BACKGROUND),ebands);
  double npbg=run->getS1S2Data(DM_CUT_DATA)->getNEvents();
  double AnalysisS1Max=run->getAnalysisS1Max();
  double AnalysisS1Min=run->getAnalysisS1Min();
  for(int s=0;s<nSlices;s++){
    double s1Mi = ebands->getS1LowerEdge(s);
    double s1Ma = ebands->getS1UpperEdge(s);
    double s1    = .5*(s1Mi+s1Ma);
    for(int b=0;b<nBands;b++){
      double fMin=b==0? -10.: 
             run->flatten(s1,s1*ebands->getS2overS1LowerEdge(b,s)); 
      double fMax=b==nBands-1? 10.:
             run->flatten(s1,s1*gbands->getS2overS1UpperEdge(b,s));  
      double g=GaussianDist::inside(fMin-meanCo60,fMax-meanCo60,sigmaCo60)
              *npbg*(s1Ma-s1Mi)/(AnalysisS1Max-AnalysisS1Min);
      gbands->fill(s,b,g);
    }
  }
  setExpectedGaussian(npbg);
  trace(EXIT,"SimplisticERBackground::computeBands");
  return gbands;
}

S1S2Bands* SimplisticERBackground::computeAnomalousBands(){
  trace(ENTRY,"SimplisticERBackground::computeAnomalousBands");
  bool spreadOverBands=isLeakageSpreadOverBands();
  double npCo60=run->getS1S2Data(E_GAMMA_CUT_DATA)->getNEvents();
  double AnalysisS1Max=run->getAnalysisS1Max();
  double AnalysisS1Min=run->getAnalysisS1Min();

  S1S2Bands* ebands=run->getS1S2Bands(E_GAMMA_CUT_DATA);
  if(ebands==NULL) {
    cout<<"Error: can't compute anomalous leakage!"<<endl;
    return NULL;
  }
  int nSlices=ebands->getNSlices();
  int nBands=ebands->getNBands();
  vector<double> anomalous(nBands);
  S1S2Bands* abands=new S1S2Bands(
      run->getNameSpace()+S1S2Data::getTypeName(ER_LEAKAGE),ebands);
  double total=0.;
  for(int b=0;b<nBands;b++){
    double predicted=0.;
    for(int s=0;s<nSlices;s++){
      double s1Mi = ebands->getS1LowerEdge(s);
      double s1Ma = ebands->getS1UpperEdge(s);
      double s1    = .5*(s1Mi+s1Ma);
      double fMin=b==0? -10.: 
             run->flatten(s1,s1*ebands->getS2overS1LowerEdge(b,s)); 
      double fMax=b==nBands-1? 10.:
             run->flatten(s1,s1*ebands->getS2overS1UpperEdge(b,s));  
      predicted+=GaussianDist::inside(fMin,fMax,sigmaCo60)*(s1Ma-s1Mi);
    }
    predicted=predicted*npCo60/(AnalysisS1Max-AnalysisS1Min);
    double nLeak=ebands->getBandContent(b);
    if(nLeak>predicted) {
      anomalous[b]=(nLeak-predicted)*normCo60;
      total+= anomalous[b];
    }
  }
  if(spreadOverBands) equalizeVector(anomalous);
  abands->fillUniformly(anomalous);
  setExpectedAnomalous(total);
  trace(EXIT,"SimplisticERBackground::computeAnomalousBands");
  return abands;
}


SimpleERBackground::SimpleERBackground() : ElectronBackgroundModel() {}
SimpleERBackground::SimpleERBackground(XeRun* r) :
  ElectronBackgroundModel("Gaussian and Anomalous evaluated",r,PL_ANALYSIS) {
  flattenedMin          = DEFAULT_FLATTENED_MIN           ;
  flattenedMax          = DEFAULT_FLATTENED_MAX           ;
  flattenedFitMin       = DEFAULT_FLATTENED_FIT_MIN       ;
  flattenedFitMax       = DEFAULT_FLATTENED_FIT_MAX       ;
  flattenedAnomalousMin = DEFAULT_FLATTENED_ANOMALOUS_MIN ;
  flattenedAnomalousMax = DEFAULT_FLATTENED_ANOMALOUS_MAX ;
}

void SimpleERBackground::setFlattenedMin(double fMin){
  flattenedMin=fMin;
}

void SimpleERBackground::setFlattenedMax(double fMax){
  flattenedMax=fMax;
}

void SimpleERBackground::setFlattenedSymmetric(double fSymmetric){
  flattenedFitMin=-fSymmetric;
  flattenedFitMax= fSymmetric;
}

void SimpleERBackground::setFlattenedFitMin(double fMin){
  flattenedFitMin=fMin;
}

void SimpleERBackground::setFlattenedFitMax(double fMax){
  flattenedFitMax=fMax;
}

void SimpleERBackground::setFlattenedFitSymmetric(double fSymmetric){
  flattenedFitMin=-fSymmetric;
  flattenedFitMax= fSymmetric;
}

void SimpleERBackground::setFlattenedAnomalousMin(double fMin){
  flattenedAnomalousMin=fMin;
}

void SimpleERBackground::setFlattenedAnomalousMax(double fMax){
  flattenedAnomalousMax=fMax;
}

void SimpleERBackground::setFlattenedAnomalousSymmetric(double fSymmetric){
  flattenedAnomalousMin=-fSymmetric;
  flattenedAnomalousMax= fSymmetric;
}

SimpleERBackground::~SimpleERBackground(){}

bool SimpleERBackground::printIt(int level){
  if(level<1) return true;
  ElectronBackgroundModel::printIt(level);
  STATIC_CONST  int tab=17;
  displayFlag("Min flattened log S2/S1",tab,flattenedMin
                     ,DEFAULT_FLATTENED_MIN);
  displayFlag("Max flattened log S2/S1",tab,flattenedMax
                     ,DEFAULT_FLATTENED_MAX);
  displayFlag("Min flattened for fit  ",tab,flattenedFitMin
                     ,DEFAULT_FLATTENED_FIT_MIN);
  displayFlag("Max flattened for fit  ",tab,flattenedFitMax
                     ,DEFAULT_FLATTENED_FIT_MAX);
  displayFlag("Min flattened anomalous",tab,flattenedAnomalousMin
                     ,DEFAULT_FLATTENED_ANOMALOUS_MIN);
  displayFlag("Max flattened anomalous",tab,flattenedAnomalousMax
                     ,DEFAULT_FLATTENED_ANOMALOUS_MAX);
  return true;
}

S1S2Bands* SimpleERBackground::computeBands(){

  trace(ENTRY,"SimpleERBackground::computeBands");

  S1S2Bands* ebands=run->getS1S2Bands(E_GAMMA_CUT_DATA);
  if(ebands==NULL) {
    cout<<"Error: can't compute gaussian leakage!"<<endl;
    return NULL;
  }

  STATIC_CONST  int nBins=100;
  double bWidth=(flattenedMax-flattenedMin)/nBins;

  TTree* tree=run->newTree(E_GAMMA_CUT_DATA);
  TH1F egamma("egamma","egamma",nBins,flattenedMin,flattenedMax);
  tree->Project("egamma","flattened");
  delete tree;
  TF1 fitEGamma("fitEGamma", "gaus", flattenedFitMin, flattenedFitMax);
  egamma.Fit("fitEGamma","QNRL");
  meanCo60=fitEGamma.GetParameter(1);
  sigmaCo60=fitEGamma.GetParameter(2);
  double fittedCo60=fitEGamma.GetParameter(0)/bWidth*sqrt(2.*PI)*sigmaCo60;
  if(debugLevel>0) {
    cout<<"SimpleERBackground  fittedCo60="<<fittedCo60<<endl;
  }

  tree=run->newTree(DM_CUT_DATA);
  TH1F dm("dm","dm",nBins,flattenedMin,flattenedMax);
  tree->Project("dm","flattened");
  delete tree;
  TF1 fitDM("fitDM", "gaus", flattenedFitMin, flattenedFitMax);
  dm.Fit("fitDM","QNRL");
  double fittedDM=fitDM.GetParameter(0)/bWidth*sqrt(2.*PI)
                 *fitDM.GetParameter(2);
  if(debugLevel>0) {
    cout<<"SimpleERBackground  fittedDM="<<fittedDM<<endl;
  }

  int nSlices=ebands->getNSlices();
  int nBands=ebands->getNBands();
  S1S2Bands* gbands=new S1S2Bands(
      run->getNameSpace()+S1S2Data::getTypeName(ER_BACKGROUND),ebands);
  double AnalysisS1Max=run->getAnalysisS1Max();
  double AnalysisS1Min=run->getAnalysisS1Min();
  for(int s=0;s<nSlices;s++){
    double s1Mi = ebands->getS1LowerEdge(s);
    double s1Ma = ebands->getS1UpperEdge(s);
    double s1    = .5*(s1Mi+s1Ma);
    for(int b=0;b<nBands;b++){
      double fMin=run->flatten(s1,s1*ebands->getS2overS1LowerEdge(b,s)); 
      double fMax=run->flatten(s1,s1*gbands->getS2overS1UpperEdge(b,s));  
      double g=GaussianDist::inside(fMin-meanCo60,fMax-meanCo60,sigmaCo60)
              *fittedDM*(s1Ma-s1Mi)/(AnalysisS1Max-AnalysisS1Min);
      gbands->fill(s,b,g);
    }
  }
  setExpectedGaussian(fittedDM);

  normCo60=fittedDM/fittedCo60;
  int bin1=1+(flattenedAnomalousMin-flattenedMin)/bWidth; 
  int bin2=1+(flattenedAnomalousMax-flattenedMin)/bWidth; 
  float outside=egamma.GetSumOfWeights()-egamma.Integral(bin1,bin2);
  float anom=outside*nBins/(nBins-bin2+bin1-1);
  setExpectedAnomalous(normCo60*anom);
  if(debugLevel>0){
    
    cout<<"SimpleERBackground outside  = "<<outside<<endl
        <<"                   bin1     = "<<bin1<<endl
        <<"                   bin2     = "<<bin2<<endl
        <<"                   anom     = "<<anom<<endl
        <<"                   Expected = "<<anomalousExpected<<endl;
  }

  trace(EXIT,"SimpleERBackground::computeBands");
  return gbands;
}

void  SimpleERBackground::computeAnomalousSlope() {

  trace(ENTRY,"SimpleERBackground::computeAnomalousSlope");

  //  variable definition
  float  HighPE	    = 30;	    // Maximum  S1
  float  LowPE	    = 3;	    // Minimum  S1
  int    n_bins_s1  = HighPE-LowPE; 
  float  HighS1S2   = 1. ;	    // Maximum log(s2/s1)	
  float  LowS1S2    = -1.;	    // Minimum log(s2/s1)	
  int    NS1S2 	    = 200; 	    // Number of bins in log(S2/s1) flattened space
 
  TH2F *h2d_logS2S1_vs_S1 =new TH2F("logS2S1_vs_S1","",n_bins_s1,LowPE,HighPE,NS1S2,LowS1S2,HighS1S2); 
  TH1D *h1d_leakage_vs_S1 =new TH1D("leakage_vs_S1","",n_bins_s1,LowPE,HighPE);


  // Collecting calibration data and filling logS2S1_vs_S1 histo
  TTree* tree=run->newTree(E_GAMMA_CUT_DATA);
  tree->Project("logS2S1_vs_S1","flattened:S1","S1 > 3. && S1 < 30. && flattened > -1. && flattened < 1.");
  delete tree;


   TF1  *F_gaus; 	// pointer to the 'current' function
   TH1D *Proj_s2s1;     // pointer to the 'current' S1 projection
 

  for(int s1_bin = 1; s1_bin <= n_bins_s1; s1_bin++)
  // Loop over slice in S1 and fit with gauss each of them. 
  // Calculate the leakage and fill it in h1d_leakage_vs_S1.
  {
	// Project each S1 slice on log(S2/S1)
  	Proj_s2s1 = h2d_logS2S1_vs_S1->ProjectionY(Form("S1_proj_%d", s1_bin), s1_bin, s1_bin);

	F_gaus    = new TF1(Form("f_gaus_%d",s1_bin),"gaus",-1.,1.);

	// Likelihood fit in range [-0.4,0.4], N=do not plot, Q= don't show results
	Proj_s2s1->Fit(F_gaus, "LNQ","",-0.4 ,0.4);
	double g_mean  = F_gaus->GetParameter(1);
	double g_sigma = F_gaus->GetParameter(2);

	double leakageLimit = 3.5*g_sigma;

	double dataTop     = Proj_s2s1->Integral(Proj_s2s1->FindBin(g_mean+leakageLimit),Proj_s2s1->FindBin(1.));  //  /!\ This returns number of events and NOT area.
	double modelTop    = F_gaus->Integral(g_mean+leakageLimit, 1.) / Proj_s2s1->GetBinWidth(1);  		//  /!\ TF1::Integral return area --> This is number of event under area.

	double dataBottom  = Proj_s2s1->Integral(Proj_s2s1->FindBin(-1.), Proj_s2s1->FindBin(g_mean - leakageLimit));
	double modelBottom = F_gaus->Integral(-1. , g_mean-leakageLimit) / Proj_s2s1->GetBinWidth(1);

	double norm = 2. / (1. - leakageLimit);  // Leakage supposed flat in log(s2/s1)

	double tot_leak_s1_slice = ( dataTop - modelTop + dataBottom - modelBottom ) * norm / 2. ; // calculate the average leakage

	h1d_leakage_vs_S1->SetBinContent(s1_bin, tot_leak_s1_slice);

	//cleaning up
	delete Proj_s2s1;
	delete F_gaus;
  }
      

   //-- Exponetial fit	
   TF1 *F_exp  = new TF1("f_exp","expo",LowPE,HighPE);
   h1d_leakage_vs_S1->Fit(F_exp,"LRNQ");
   double slope = -1.* F_exp->GetParameter(1); // the slope is used in other methods as a positive parameter
  
   setAnomalousSlope(slope);

   //--- some plotting ---//
//   h1d_leakage_vs_S1->Draw();
//   cout<< "Anomalous slope " << slope << endl;

   //cleaning up
   delete h1d_leakage_vs_S1;
   delete h2d_logS2S1_vs_S1;
   delete F_exp;

}



S1S2Bands* SimpleERBackground::computeAnomalousBands(){

  trace(ENTRY,"SimpleERBackground::computeAnomalousBands");

  computeAnomalousSlope();  // Ale new

  S1S2Bands* ebands=run->getS1S2Bands(E_GAMMA_CUT_DATA);
  S1S2Bands* abands=new S1S2Bands(
      run->getNameSpace()+S1S2Data::getTypeName(ER_LEAKAGE),ebands);
  int nSlices=abands->getNSlices();
  int nBands=abands->getNBands();
  double w=anomalousExpected/abands->getArea();
  
  double wtag=1.;
  if(anomalousSlope>0.) {
    double S1Mi = ebands->getS1LowerEdge();
    double S1Ma = ebands->getS1UpperEdge();
    // /!\ Ale: I don't get this!!! 
    wtag=(S1Ma-S1Mi)/(exp(-anomalousSlope*S1Mi)-exp(-anomalousSlope*S1Ma));
  }  

  for(int s=0;s<nSlices;s++){
    double s1Mi = ebands->getS1LowerEdge(s);
    double s1Ma = ebands->getS1UpperEdge(s);
    double s1    = .5*(s1Mi+s1Ma);
    for(int b=0;b<nBands;b++){
      double fMin= run->flatten(s1,s1*ebands->getS2overS1LowerEdge(b,s)); 
      double fMax= run->flatten(s1,s1*ebands->getS2overS1UpperEdge(b,s));  
      double e=1.;
      if(anomalousSlope>0.) {
        e=exp(-anomalousSlope*s1Mi)-exp(-anomalousSlope*s1Ma);  
      }
      abands->fill(s,b,w*wtag*e*(s1Ma-s1Mi)*(fMax-fMin));
    }
  }
  trace(ENTRY,"SimpleERBackground::computeAnomalousBands");
  return abands;
}


/*PublishedElectronBackground for Run 8 */

bool PublishedElectronBackgroundRun8::isRunCompatible(int r) {return r==8;}
PublishedElectronBackgroundRun8::~PublishedElectronBackgroundRun8()
  {
//	delete gausBands;  // somehow if I add this lines it crashes... But it shouldent... FIXME
//	delete leakBands;
//	delete fsigmaco60; 	  
  }

PublishedElectronBackgroundRun8::PublishedElectronBackgroundRun8() 
  : ElectronBackgroundModel() {
//----- Initialize ---//
  gaussianExpected   = 896.5;		// According to note https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat3
  anomalousExpected  = 0.;		// to be computed
  anomalousSlope     = 0.;  		// the anomalous background was considered flat in S1
  normCo60	     = 1./43.3; 	//Factor to convert from the Co60 exposure to the DM exposure
  leak_RangeFraction = 26./32.;		//Factor to correct from the 3-35pe range where the leakage is computed to the actual 4-30pe range
  npco60	     = 17595;    	//Number of Co60 points passing all cuts between 3. and 35. pe in S1
  meanCo60	     = 0.022242;        // Constants for gaussian are taken from run8 note. 
  sigmaCo60	     = 0.135846;	// Constants for gaussian are taken from run8 note. 


  fsigmaco60 = new TF1("fsigmaco60", "(-0.438131/3.09+0.000597874/3.09*x)", 0, 50);  //Sigma of the Co60 band as a function of S1
   										     //3.09 is the number of signals in the 99.8 area. dividing by 3.09 gives the 1 sigma level.
  gausBands = NULL;      // bands for gaussian component of ER bkg
  leakBands = NULL;    	 // bands for leakage component of ER bkg

}

PublishedElectronBackgroundRun8::PublishedElectronBackgroundRun8(XeRun* r) 
   :  ElectronBackgroundModel("Published for run 8",r,PL_ANALYSIS) {

/*  Just for sake of completeness, reported here the old outcome from Offer's bands.
    XEPHYR produce slightly different bands from data-driven method, while Offer uses fitted functions.
    The outcome at the end is a bit different but compatible.
    Note: XEPHYR applies an additional cut on  20. < s1/s2 < 500. while building the bands, this leads to
    a rejection of small fraction of bkg and signal. The "gaussianExpected" value should not be modified
    and will give the correct result within this cut.
	
      	gausleak    anomleak    		neutleak    	Tot
      	843.809		0   			0.0256585   =	843.835  (band 11)
      	23.9303 	0.834164   		0.0216748   = 	24.7862  (band 10)
      	11.466  	0.602256   		0.0226586   =	12.091   (band 9 )
      	6.64066 	0.215476   		0.0235288   =	6.87967  (band 8 )
 	4.14485  	0.387094   		0.0248152   =	4.55676  (band 7 )
	2.66158  	0.194321   		0.0265611   =	2.88246  (band 6 )
	1.70718  	0.0259269   		0.0274043   =	1.76051  (band 5 )
	1.06426  	0.0216907   		0.0287881   =	1.11474  (band 4 )
	0.621569  	0.0319206   		0.0295772   =	0.683066 (band 3 )
	0.317961  	0.158781   		0.0308258   =	0.507568 (band 2 )
	0.11999  	0.128686   		0.0272043   =	0.27588  (band 1 )
	0.0166444 	0.25089   		0.0168588   =	0.284393 (band 0 ) 


*/

//----- Initialize According to note https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat3  ---//
  gaussianExpected   = 896.5;		// expected events for gaussian ER 
  anomalousExpected  = 0.;		// to be computed
  anomalousSlope     = 0.;  		// the anomalous background was considered flat in S1
  normCo60	     = 1./43.3; 	// Factor to convert from the Co60 exposure to the DM exposure
  leak_RangeFraction = 26./32.;		// Factor to correct from the 3-35pe range where the leakage is computed to the actual 4-30pe range
  npco60	     = 17595;    	// Number of Co60 points passing all cuts between 3. and 35. pe in S1
  meanCo60	     = 0.022242;        // Constants for gaussian are taken from run8 note. 
  sigmaCo60	     = 0.135846;	// Constants for gaussian are taken from run8 note. 


  fsigmaco60 = new TF1("fsigmaco60", "(-0.438131/3.09+0.000597874/3.09*x)", 0, 50);  //Sigma of the Co60 band as a function of S1
   										     //3.09 is the number of signals in the 99.8 area. dividing by 3.09 gives the 1 sigma level.
										     //Don't know why but sigma is returned negative.
  gausBands = NULL;      // bands for gaussian component of ER bkg
  leakBands = NULL;    	 // bands for leakage component of ER bkg

  isGausComputed = false;
  isLeakComputed = false;

}

bool PublishedElectronBackgroundRun8::printIt(int level){
  ElectronBackgroundModel::printIt(level);
  return true;
}


S1S2Bands* PublishedElectronBackgroundRun8::getGausBands()
{
   if(isGausComputed) return gausBands;
   return computeBands();
}

S1S2Bands* PublishedElectronBackgroundRun8::getLeakBands()
{
   if(isLeakComputed) return leakBands;
   return computeAnomalousBands();
}
  
S1S2Bands* PublishedElectronBackgroundRun8::computeBands()
{
 //Gaussian Contribution  ////////////////////////////////////
 // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat3
 // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:code_implementation
 ////////////////////////////////////////////////////////////

  trace(ENTRY,"PublishedElectronBackgroundRun8::computeBands");

  //if(isGausComputed) return gausBands;  // this function should not return a pointer if you want this feature. FIXME 

  // load REFERENCE bands
  S1S2Bands* refBands=run->getReferenceBands();  

  if(refBands==NULL) {
    cout<<"PublishedElectronBackgroundRun8::computeBands - Error - Reference bands are empty, return NULL"<<endl
	<<"Run: " << run->getName() << " may not be initialize, make shure that you executed XeRun::fillDataBands().  Return NULL."<<endl;
    return NULL;
  }

  // bkg bands which will be filled
  gausBands = new S1S2Bands( run->getNameSpace()+S1S2Data::getTypeName(ER_BACKGROUND),refBands,false);
  
  int    nSlices       = refBands->getNSlices();
  int    nBands        = refBands->getNBands();
  double AnalysisS1Max = run->getAnalysisS1Max();
  double AnalysisS1Min = run->getAnalysisS1Min();

  double testCounter =0. ; // this counter should be == gaussianExpected
			   // Tourn out to be slightly different due to the cut on s1/s2
			   // in band making: s1/s2 < 500 and s1/s2 > 20. 

  // Loop over S1 slices
  for(int s=0;s<nSlices;s++)
  {
    double s1Mi  =  refBands->getS1LowerEdge(s);
    double s1Ma  =  refBands->getS1UpperEdge(s);
    double s1    =  .5*(s1Mi+s1Ma);
   
    // loop over bands
    for(int b=0; b<nBands; b++)
    {
      //get the edges of the band in the flatten space for S1 slice = s	
      double BandMin = run->flatten(s1, s1*refBands->getS2overS1LowerEdge(b,s)); 
      double BandMax = run->flatten(s1, s1*refBands->getS2overS1UpperEdge(b,s));  

      double gaussExp =  gaussianExpected *(s1Ma-s1Mi)/(AnalysisS1Max-AnalysisS1Min)* 0.5 *  //norm factors: the ER is considered flat in S1 in flatten space, 0.5 because of TMath::Erf()
		(TMath::Erf( (BandMax-meanCo60)/sqrt(2.)/sigmaCo60 ) - TMath::Erf( (BandMin-meanCo60)/sqrt(2.)/sigmaCo60) ); // integral of gauss on band
      // Note --- TMath::Erf(X)  = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x 

      testCounter += gaussExp; 

      gausBands->fill(s,b,gaussExp);
    }
  }

  if(debugLevel > 0)
	cout << "PublishedElectronBackgroundRun8::computeBands() " << "- Debug -  Bkg testCounter = " << testCounter << " Expected = " << gaussianExpected << endl; 

  isGausComputed = true; // setting the flag.  

  trace(EXIT,"PublishedElectronBackgroundRun8::computeBands");
  return gausBands;

}

S1S2Bands* PublishedElectronBackgroundRun8::computeAnomalousBands()
{
  //  Anomalous Leakage Contribution ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat2#uncertainty_in_the_er_rejection
  // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:erleakagerun08flat3
  // https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:bgprediction:code_implementation
  //  
  // * Comment: leakage computed per band and s1 bin, it is averaged between 3 < S1 < 35 per "flat" slice, which is conceptually wrong
  //   since the gauss sigma changes with S1. However this is the official method for run8.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  trace(ENTRY,"PublishedElectronBackgroundRun8::computeAnomalousBands");

  //if(isLeakComputed) return leakBands;  // this function should not return a pointer if you want this feature. FIXME 

  // load REFERENCE bands
  S1S2Bands* refBands      = run->getReferenceBands();  

  // Co60 calibration data points with all cuts except S1. 
  S1S2Data*  ERcalibration = run->getS1S2Data(E_GAMMA_DATA);  

  if(refBands == NULL || ERcalibration == NULL) {
    cout<<"PublishedElectronBackgroundRun8::computeAnomalousBands() - Error - Reference bands are not defined. "<< endl
 	<< "Run: " << run->getName() << " may not be initialize, make shure that you executed XeRun::fillDataBands().  Return NULL."<<endl;
    return NULL;
  }

  //the bands that are going to be filled
  leakBands= new S1S2Bands(run->getNameSpace()+S1S2Data::getTypeName(ER_LEAKAGE) ,refBands,false);

  int    nSlices       = refBands->getNSlices();
  int    nBands        = refBands->getNBands();
  double AnalysisS1Max = run->getAnalysisS1Max();
  double AnalysisS1Min = run->getAnalysisS1Min();

  // Loop over S1 slices (this is crazy inefficient! try to FIXME!)
  for(int s=0;s<nSlices;s++)
  {
    double s1Mi  =  refBands->getS1LowerEdge(s);
    double s1Ma  =  refBands->getS1UpperEdge(s);
    double s1    =  .5*(s1Mi+s1Ma);
   
    // loop over bands
    for(int b=0; b<nBands; b++)
    {
      //get the edges of the band in the flatten space for S1 slice = s	
      double BandMin = run->flatten(s1, s1*refBands->getS2overS1LowerEdge(b,s)); 
      double BandMax = run->flatten(s1, s1*refBands->getS2overS1UpperEdge(b,s));  
      
      // counter of events that fall in the flatten "stripe" related to the band.	
      int evnt_in_stripe = 0;
	
      //loop over the events in the dataset, counting how many in that band "stripe" for 3 < S1 < 35 
      for(int event_itr =0; event_itr < ERcalibration->getNEvents(); event_itr++)
      {
	  double evnt_s1         = ERcalibration->getS1(event_itr);
	  double evnt_flat_s2os1 = run->flatten(evnt_s1, ERcalibration->getS2(event_itr));
	  if( evnt_s1 > 3. && evnt_s1 < 35. && evnt_flat_s2os1 > BandMin && evnt_flat_s2os1 < BandMax ) evnt_in_stripe++; 
      }
	// integral of gauss on stripe 
      double gaussExp =  npco60 *  0.5 * (TMath::Erf( BandMax/sqrt(2.)/(-1.*fsigmaco60->Eval(s1)) ) - TMath::Erf( BandMin /sqrt(2.) /(-1.*fsigmaco60->Eval(s1)) ) ); 
      // Note --- TMath::Erf(X)  = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x
      // Note -- fsigmaco60 returns negative values (don't know why) 

      // leakage for this S1 bin in this band	
      double leak_per_bin =0.;
      if( gaussExp < evnt_in_stripe )                        //norm for s1 bin                           //norm 3-35 PE      // norm to data   	
	leak_per_bin = ((double)evnt_in_stripe - gaussExp) * (s1Ma-s1Mi)/(AnalysisS1Max-AnalysisS1Min) * leak_RangeFraction * normCo60; 

	leakBands->fill(s,b,leak_per_bin);
    }
   }

//  setExpectedAnomalous(vectorSum(anomalousBkgs));
  trace(EXIT,"PublishedElectronBackgroundRun8::computeAnomalousBands");

  isLeakComputed = true; // setting the flag.

  return leakBands;
}








bool PublishedElectronBackgroundRun10::isRunCompatible(int r) {return r==10;}
PublishedElectronBackgroundRun10::~PublishedElectronBackgroundRun10(){}


PublishedElectronBackgroundRun10::PublishedElectronBackgroundRun10() 
                            : ElectronBackgroundModel() {}
PublishedElectronBackgroundRun10::PublishedElectronBackgroundRun10(XeRun* r)
 :  ElectronBackgroundModel("Published for run 10",r,PL_ANALYSIS) {
 
  double PRL_ANOMALOUS_SLOPE=0.225;
  
  double bg[DEFAULT_N_BANDS_RUN_10]={
           0.193037 ,0.343531 ,0.481669 ,0.731189 ,1.01777 ,1.30644
          ,1.8196 ,2.71279 ,3.98278 ,6.50712 ,14.7515 ,375.961 };

  double abg[DEFAULT_N_BANDS_RUN_10]={
           0.422222 ,0.0950062 ,0.0580075 ,0.0478701 ,0.0411387 ,0.0369015
          ,0.035904 ,0.0357247 ,0.0371685 ,0.0436412 ,0.0592695 ,1.3163};

  bkgs.resize(DEFAULT_N_BANDS_RUN_10);
  setBackground(bg);
  anomalousBkgs.resize(DEFAULT_N_BANDS_RUN_10);
  setAnomalousBackground(abg);
  setAnomalousSlope(PRL_ANOMALOUS_SLOPE);
}

bool PublishedElectronBackgroundRun10::printIt(int level){
  ElectronBackgroundModel::printIt(level);
  return true;
}

void PublishedElectronBackgroundRun10::setBackground(double *backgrounds){
  for(int b=0;b<DEFAULT_N_BANDS_RUN_10;b++) bkgs[b]=backgrounds[b];
}

void PublishedElectronBackgroundRun10::setAnomalousBackground(double *backgrounds){
  for(int b=0;b<DEFAULT_N_BANDS_RUN_10;b++) anomalousBkgs[b]=backgrounds[b];
}


S1S2Bands* PublishedElectronBackgroundRun10::computeBands(){
  trace(ENTRY,"PublishedElectronBackgroundRun10::computeBands");
  S1S2Bands* refBands=run->getReferenceBands();
  if(refBands==NULL) {
    cout<<"Error: Can't compute electron gaussian bands!"<<endl;
    return NULL;
  }
  int nBands=refBands->getNBands();
  if(nBands!=(int)bkgs.size()) {
    cout<<"Mismatch in number of bands in "<<Name<<endl;
  }
  S1S2Bands* GaussERBands=new S1S2Bands(
       run->getNameSpace()+S1S2Data::getTypeName(ER_BACKGROUND),refBands,false);
  if(ERGaussianFlatInS1){
    GaussERBands->fillUniformly(bkgs);
  }
  else {
    TabulatedDist *dist=run->getS1S2Bands(E_GAMMA_CUT_DATA)->newSliceXeDist();
    GaussERBands->fillAccordingToS1Dist(bkgs,dist->getSpectrum());
    delete dist;
  }
  setExpectedGaussian(vectorSum(bkgs));
  trace(EXIT,"PublishedElectronBackgroundRun10::computeBands");
  return GaussERBands;
}

S1S2Bands* PublishedElectronBackgroundRun10::computeAnomalousBands(){
  trace(ENTRY,"PublishedElectronBackgroundRun10::computeAnomalousBands");
  S1S2Bands* refBands=run->getS1S2Bands(AM_BE_CUT_DATA);
  if(refBands==NULL) {
    cout<<"Error: Can't compute anomalous electron bands!"<<endl;
    return NULL;
  }
  int nBands=refBands->getNBands();
  if(nBands!=(int)anomalousBkgs.size()) {
    cout<<"Mismatch in number of bands in "<<Name<<endl;
  }
  string title=run->getNameSpace()+S1S2Data::getTypeName(ER_LEAKAGE);
  S1S2Bands* ERAnomalousBands=new S1S2Bands(title,refBands,false);
  if(anomalousSlope>0.) {
    ERAnomalousBands->fillExponential(anomalousBkgs,anomalousSlope);
  }
  else {
    ERAnomalousBands->fillUniformly(anomalousBkgs);
  }
  setExpectedAnomalous(vectorSum(anomalousBkgs));
  trace(EXIT,"PublishedElectronBackgroundRun10::computeAnomalousBands");
  return ERAnomalousBands;
}

//------------------- Profile likelihood based on S1S2 run ----------------

S1S2PL::~S1S2PL(){



}


S1S2PL::S1S2PL(XeRun *r) : ProfileLikelihood() , RunComponent(){
if(printLevel > 1) cout << "S1S2PL::S1S2PL() - DEBUG : ENTERING" << endl;

  setRun(r);

  setName("S1S2 Profile likelihood on "+run->getName());

  // this is for initialize the XeRun, in case it is not.
  run->checkIt();  // NOTE: it is important that this method is called before initialize()

  if(!run->isBandFromMC()) {cout << "S1S2PL::S1S2PL() - WARNING: bands are computed from AmBe!!! "<< endl;
  			run->tabulateSignalLeffQy(); // this would not be called otherwise
  }

  initialize();

  run->setPValue(this);

  run->markAsChanged();  // NOTE: it is important that this method is called before initialize()
  
  run->checkIt();

  checkRunCompatibility();

 
  //update();
 
if(printLevel > 1) cout << "S1S2PL::S1S2PL() - DEBUG : LEAVING" << endl;
}

void S1S2PL::initialize(){
if(printLevel > 1) cout << "S1S2PL::initialize() - DEBUG : ENTERING" << endl;

  nBands = run->getReferenceBands()->getNBands();

  // setting flag to default
  isLEffTValueFit	=	DEFAULT_PAR_LEFF_TVALUE_FIT;
  isQyTValueFit		=	DEFAULT_PAR_QY_TVALUE_FIT;
  S1Likelihood		=	DEFAULT_S1_LIKELIHOOD;
  isSystBkgTValueFit	=	DEFAULT_PAR_SYST_BKG_FIT;
  isStatBkgTValueFit	=	DEFAULT_PAR_STAT_BKG_FIT;
  isSigAccFit           =       false;

  // initializing post fit parameter values to UNDEFINED
  SigmaForExclusion	     =  UNDEFINED; 
  LEffTValueForExclusion     =  UNDEFINED; 	
  DLEffTValueForExclusion    =  UNDEFINED; 
  QyTValueForExclusion	     =  UNDEFINED; 
  DQyTValueForExclusion	     =  UNDEFINED; 
  SystBkgTValueForExclusion  =  UNDEFINED; 
  DSystBkgTValueForExclusion =  UNDEFINED; 

  setWithSimulatedata(false);

  // initializing bkg normalization sigma, default uncertainty as if one use the last 
  // band (bkg dominated) to normalize. NOTE: if N_bkg = N_data_CR / N_bkg_CR * N_bkg_tot
  // then sigma_N_bkg ~= sqrt(N_data_CR) * N_bkg_tot / N_bkg_CR
  // This should be instead given by whom compute the norm factor FIXME.
  double N_data_CR = run->getS1S2Bands(DM_CUT_DATA)->getBandContent(nBands -1);   // counting starts from zero that's why -1
  N_bkg_tot = run->getS1S2Bands(ALL_BACKGROUNDS)->getCumulatedBandContent(nBands-1);
  double N_co_CR  = run->getS1S2Bands(E_GAMMA_CUT_DATA)->getBandContent(nBands -1);


  double N_Co_total =  run->getS1S2Bands(E_GAMMA_CUT_DATA)->getCumulatedBandContent(nBands-1); 

  SystBkgSigma	= sqrt( 1./ N_data_CR + 1./ N_co_CR) * N_data_CR / N_co_CR * N_Co_total ; 

  if(printLevel > -1) {
	cout << "S1S2PL::initialize()  -  INFO : " <<getName() << " parameters values " << endl;
	cout << "                      SystBkgSigma = " << SystBkgSigma << endl;
	cout << "                      N_bkg_tot    = " << N_bkg_tot    << endl;
	cout << "                      N_Co_total    = " << N_Co_total    << endl;
  }
  
  //initializing nuissance parameter
  ls = new SigmaParameter();
  addParameter(ls);
  lkSystBkg = new TSystBkgParameter(run->getNumber());
//  lkSystBkg->setInitialValue(-0.5);
  addParameter(lkSystBkg);
  le = new TLEffParameter();
  le->setMaximum(run->getLEff()->getMaximumTval());
  le->setMinimum(run->getLEff()->getMinimumTval());
//  le->setStep(run->getLEff()->getStepSize()); 
  addParameter(le);
  lqy = new TQyParameter();
  lqy->setMaximum(run->getQy()->getMaximumTval());
  lqy->setMinimum(run->getQy()->getMinimumTval()); 
//  lqy->setStep(run->getQy()->getStepSize()); 
  addParameter(lqy);

  // signal acceptance nuissance
  lkacc = new LKParameter(PAR_SIGN_ACCEPTANCE_TVALUE, NUISANCE_PARAMETER, "signalAcc",0.,0.01,-2.,2.);
  addParameter(lkacc);  

  for(int b = 0; b < nBands; b++){
	//initializing post fit parameter values to UNDEFINED
	StatBkgTValueForExclusion.push_back(UNDEFINED);
	DStatBkgTValueForExclusion.push_back(UNDEFINED);
	
	//initializing stat parameters
	StatBkgSigmas.push_back(0.);  // CHECKME FIXME
        lkStatBkgs.push_back(new TStatBkgParameter(b, run->getNumber()));
        addParameter(lkStatBkgs[b]);
	
  } 

  DataBand 	 = run->getS1S2Bands(DM_CUT_DATA);

  data_type      = DM_CUT_DATA;			// default value 

  BackgroundBand = run->getS1S2Bands(ALL_BACKGROUNDS);


if(printLevel > 1) cout << "S1S2PL::initialize() - DEBUG : LEAVING" << endl;

}



void          S1S2PL::fitLEffTValue(bool fit)      {isLEffTValueFit=fit;}
void          S1S2PL::fitQyTValue(bool fit)        {isQyTValueFit=fit;}
void          S1S2PL::fitSystBkgTValue(bool fit)   {isSystBkgTValueFit=fit;}
void          S1S2PL::fitStatBkgTValue(bool fit)   {isStatBkgTValueFit=fit;}
void          S1S2PL::withS1Likelihood(bool s1)    {S1Likelihood=s1;}
void          S1S2PL::setSystBkgSigma(double s)    {SystBkgSigma=s;}
LKParameter*  S1S2PL::getLEffTValueParameter()     {return le;}
LKParameter*  S1S2PL::getQyTValueParameter()       {return lqy;}
LKParameter*  S1S2PL::getSigmaParameter()          {return ls;}
LKParameter*  S1S2PL::getSystBkgTValueParameter()  {return lkSystBkg;}
vector<LKParameter*> S1S2PL::getStatBkgTValueParameters() {return lkStatBkgs;}

double  S1S2PL::getSigmaForExclusion()         {return SigmaForExclusion;}
double  S1S2PL::getLEffTValueForExclusion()    {return LEffTValueForExclusion;}
double  S1S2PL::getDLEffTValueForExclusion()   {return DLEffTValueForExclusion;}
double  S1S2PL::getQyTValueForExclusion()    {return QyTValueForExclusion;}
double  S1S2PL::getDQyTValueForExclusion()   {return DQyTValueForExclusion;}
double  S1S2PL::getSystBkgSigma()              {return SystBkgSigma;}

double S1S2PL::getSystBkgTValueForExclusion()  {
  return SystBkgTValueForExclusion; 
}
double S1S2PL::getDSystBkgTValueForExclusion() {
  return DSystBkgTValueForExclusion;
}
vector <double> S1S2PL::getStatBkgTValueForExclusion()  {
  return StatBkgTValueForExclusion; 
}
vector<double> S1S2PL::getDStatBkgTValueForExclusion() {
  return DStatBkgTValueForExclusion;
}

void S1S2PL::exclusionComputed(){
if(printLevel > 1) cout << "S1S2PL::exclusionComputed() - DEBUG : ENTERING" << endl;

  cout << "\n\n---------------------------S1S2PL::exclusionComputed(ENTER)-------------------------------" << endl;
  SigmaForExclusion=ls->getCurrentValue();
  if(printLevel>0) {
    cout<<"Last value for sigma          :"<<SigmaForExclusion<<endl;  
  }
  if(isLEffTValueFit){
    LEffTValueForExclusion=le->getCurrentValue();
    DLEffTValueForExclusion=le->getSigma();
    if(printLevel>0) {
      cout<<"Last value for LEff t-value   :"<<LEffTValueForExclusion
          <<" +/- "<<DLEffTValueForExclusion<<endl;
    }
  }

  if(isSystBkgTValueFit){
    SystBkgTValueForExclusion=  lkSystBkg->getCurrentValue();
    DSystBkgTValueForExclusion= lkSystBkg->getSigma();
    if(printLevel>0) {
      cout<<"Last value for norm bkg t-valu:"<<SystBkgTValueForExclusion
          <<" +/- "<<DSystBkgTValueForExclusion<<endl;
    }
  }

  if(isQyTValueFit){
    QyTValueForExclusion=lqy->getCurrentValue();
    DQyTValueForExclusion=lqy->getSigma();
    if(printLevel>0) {
      cout<<"Last value for Qy t-value   :"<<QyTValueForExclusion
          <<" +/- "<<DQyTValueForExclusion<<endl;
    }
  }

  if(isStatBkgTValueFit){
     for(int b = 0; b < nBands; b++){
	//initializing post fit parameter values to UNDEFINED
	StatBkgTValueForExclusion[b] =  lkStatBkgs[b]->getCurrentValue();
	DStatBkgTValueForExclusion[b] = lkStatBkgs[b]->getSigma();
        if(printLevel>0) {
            cout<<"Last value for StatBkg t-value band "<<b << " : " << StatBkgTValueForExclusion[b] 
          <<" +/- "<<DStatBkgTValueForExclusion[b] <<endl;
         }
	
       } 
  }

if(printLevel > 1) cout << "S1S2PL::exclusionComputed() - DEBUG : LEAVING" << endl;
  cout << "---------------------------S1S2PL::exclusionComputed(END)-------------------------------\n\n" << endl;
}

bool S1S2PL::update() {
if(printLevel > 1) cout << "S1S2PL::update() - DEBUG : ENTERING" << endl;
  //NOTE: This function is called when the wimp mass is changed
  //variables are set to their default again.
/*
  if(printLevel>2) {
     cout<<"Updating "<<Name<<endl;
  } 
  le=getParameter(PAR_LEFF_TVALUE);
  if(le==NULL){
    cout<<"No LEFF t-value parameter in "<<getName()<<endl;
    return false; 
  }
  lqy=getParameter(PAR_QY_TVALUE);
  if(lqy==NULL){
    cout<<"No QY t-value parameter in "<<getName()<<endl;
    return false; 
  }
*/

  if(run==NULL){
    cout<<"No XeRun has been defined, can't set up "<<endl;
    return false;
  }
  else if(run->withoutData()) {
    cout<<run->getName()<<" has no data, can't set up "<<getName()<<endl;
    return false;
  }
  ls=getParameter(PAR_SIGMA);
  if(ls==NULL){
    cout<<"No sigma parameter in "<<getName()<<endl;
    return false; 
  }

  //here the bands are recomputed and Signal is retabulated 
  run->update(); 
  if(!run->isBandFromMC()) {cout << "S1S2PL::S1S2PL() - WARNING: bands are computed from AmBe!!! "<< endl;
    			//run->tabulateSignalLeffQy(); // this would not be called otherwise
  }

  //Pointer is restored to the recomputed bkg bands
  BackgroundBand = run->getS1S2Bands(ALL_BACKGROUNDS);
  DataBand = run->getS1S2Bands(data_type);


  if(printLevel>1) {
    cout<<"Updating "<<Name;
    if(!ls->inCombinedMode()) cout<<" sigmaUnit="<<sigmaUnit;
    if(isSystBkgTValueFit) cout<<", syst bkg scaling="<<SystBkgSigma;
    cout<<endl;
  }


  if(!ls->inCombinedMode())  ls->setMinuitUnit(1.);  // Ale change sigma unit
  //if(!ls->inCombinedMode())  ls->setMinuitUnit(sigmaUnit);

  activateParameter(lkSystBkg,isSystBkgTValueFit);

  activateParameter(le,isLEffTValueFit);

  activateParameter(lqy,isQyTValueFit);

  activateParameter(lkacc,isSigAccFit);

  for(int b=0;b<nBands;b++) {
    //epsilons nuissance parameter (stat. uncertainties)
    lkStatBkgs[b]->setInitialValue(BackgroundBand->getBandContent(b) / N_bkg_tot);
//THIS IS WRONG!!! TEST_ALE
    ((TStatBkgParameter*)lkStatBkgs[b])->setStatError( sqrt(BackgroundBand->getBandContent(b)) /  N_bkg_tot );
    lkStatBkgs[b]->setStep(sqrt(BackgroundBand->getBandContent(b)) /  N_bkg_tot *0.01 );
    //((TStatBkgParameter*)lkStatBkgs[b])->setStatError( sqrt(BackgroundBand->getBandContent(b) * frac_Co_model) / ( N_bkg_tot  * frac_Co_model));
    //StatBkgSigmas[b] = sqrt(BackgroundBand->getBandContent(b) * frac_Co_model) / ( N_bkg_tot  * frac_Co_model);
    activateParameter(lkStatBkgs[b],isStatBkgTValueFit);

  }


  run->setSigmaNucleon(DEFAULT_SIGMA_NUCLEON);  //NOT sure of this CHECKME FIXME
//  run->computeSignal();  // Not sure needed              	CHECKME FIXME

  if(printLevel > 1) cout << "S1S2PL::update() - DEBUG : LEAVING" << endl;
  return run->checkAll();
}





void S1S2PL::setWimpMass(double mass) {
  if(printLevel>1) cout<<getName()<<" sets wimp mass to "<<mass<<endl;
  run->setWimpMass(mass);
}

double S1S2PL::getWimpMass() {
  return  run->getWimpMass();
}

bool S1S2PL::simulate(double sigma) { 
  run->simulate(sigma); 
  return true;
}





double S1S2PL::computeTheLogLikelihood(){

//NOTE: here likelihood is considered as the Log(likelihood) and it is maximized afterwards

//Retriving Parameter of Interest
  double sigma=getParameterValue(PAR_SIGMA);

//Retriving Nuissace parametercurrent fitting  t_value
  double t_NormBkg 	  = 	0.;
  double t_Leff		  = 	0.;
  double t_Qy 	 	  = 	0.;
  double t_epsilon_Bkg_b  = 	0.;
  double t_sig_acc  	  = 	0.;

  if(isLEffTValueFit)      t_Leff 	=   getParameterValue(PAR_LEFF_TVALUE);
  if(isQyTValueFit)	   t_Qy		=   getParameterValue(PAR_QY_TVALUE);
  if(isSystBkgTValueFit)   t_NormBkg   	=   getParameterValue(PAR_SYST_BKG_TVALUE);
  if(isSigAccFit)          t_sig_acc   	=   getParameterValue(PAR_SIGN_ACCEPTANCE_TVALUE);

//Check limit status
  if(std::isnan(t_Leff)){
    cout<<"Fit is unstable for NP Leff "<<endl;
    return VERY_SMALL;
  }
  if(std::isnan(t_Qy)){
    cout<<"Fit is unstable for NP Qy "<<endl;
    return VERY_SMALL;
  }
  if(std::isnan(t_NormBkg)){
    cout<<"Fit is unstable for NP NormBkg"<<endl;
    return VERY_SMALL;
  }

  if(std::isnan(t_sig_acc)){
    cout<<"Fit is unstable for NP signAcc"<<endl;
    return VERY_SMALL;
  }
  

  double LL=0;

  double N_Co_total =  run->getS1S2Bands(E_GAMMA_CUT_DATA)->getCumulatedBandContent(nBands-1) ; 
//  cout << "N_Co_total  " << N_Co_total << endl;
  // adding NP gaussian constraints, NOTE: we neglect additive constant factors
//  if(isSystBkgTValueFit)    LL 	= LL - t_NormBkg*t_NormBkg / 2.;  //just for test with hagar UNCOMMENT
  if(isQyTValueFit)         LL  = LL - t_Qy*t_Qy / 2.;
  if(isLEffTValueFit)       LL  = LL - t_Leff*t_Leff / 2.;
  if(isSigAccFit)           LL  = LL - t_sig_acc*t_sig_acc /2.;


 S1S2Bands *SignalBand_t = run->getSignalBand(t_Leff,t_Qy);
// S1S2Bands *Data = run->getS1S2Bands(ALL_BACKGROUNDS);  // ALE_TEST TEST_A
// S1S2Bands *Data = run->getS1S2Bands(DM_SIMULATED_DATA);  // ALE_TEST TEST_A
// S1S2Bands *Data = run->getS1S2Bands(DM_CUT_DATA);  // ALE_TEST TEST_A

 if(t_Leff !=0 || t_Qy !=0)  
     SignalBand_t = new S1S2Bands( "", run->interpolateLeffQy(t_Leff, t_Qy), true);

      if(debugLevel > 1) cout << "Leff "<< t_Leff << " Qy  " << t_Qy << " t_NormBkg " << t_NormBkg << " sigma " << sigma << endl;
 
  // loop over the bands 
  for(int b=0; b < nBands; b++){

    	  double LL_band = 0.;

      //number of event for signal in this band
      //double N_s       = sigma * 0.2332;  // test to get equal result as OLD style XEPHYR 
      double N_s       = sigma *  run->getSignalHandler()->getXsecMultiplier() * SignalBand_t->getBandContent(b);  // normalized to 10^-45 cm^2
      if(isSigAccFit) N_s = N_s * ( 1. + t_sig_acc * 0.12);

      // NOTE that t_epsilon_bkg_b is not a t_value, it is actually used as parameter
      // that goes between [0,1], this to avoid problems with negative values and when
      // the actual epsilon_bkg_b = 0
      double epsilon_bkg_b = BackgroundBand->getBandContent(b) / N_bkg_tot;
      t_epsilon_Bkg_b      = epsilon_bkg_b;
      if(isStatBkgTValueFit) t_epsilon_Bkg_b  = getParameterValue(PAR_STAT_BKG_TVALUE+b);

		          //Tot bkg events + norm Syst uncertainty  // fraction of bkg event per band b + Stat. uncetrainty
      double N_b       =   (N_bkg_tot + t_NormBkg*SystBkgSigma)  *  t_epsilon_Bkg_b ;

//      if(b == 0 || b ==1) N_b = N_b /2.;

      double N_obs     = DataBand->getBandContent(b); // data are specified with setData, default in DM_DATA
//      if( (b==0 || b==1) && data_type ==ASIMOV_DATA) N_obs = N_obs /2.;  

      double N_Co_b = run->getS1S2Bands(E_GAMMA_CUT_DATA)->getBandContent(b);

      if(debugLevel > 1 && data_type != ASIMOV_DATA) cout << "Band "<< b << " N_S  " << N_s << " sigma " << sigma << "  Histo content " << SignalBand_t->getBandContent(b) 
      			       << " N_b  " << N_b <<  " N_obs " << N_obs  << "  Nco_b  " << N_Co_b << "  t_epsilon_Bkg_b " << t_epsilon_Bkg_b<< endl;

     //Hagar style 
     if(N_s + N_b <= 0.) return VERY_SMALL;
     LL_band += N_obs*log(N_s + N_b) -1.*( N_s + N_b )  ;
     //LL_band += N_Co_b*epsilon_bkg_b * log(N_Co_total*t_epsilon_Bkg_b) - (N_Co_total*t_epsilon_Bkg_b);
     LL_band += N_Co_total*epsilon_bkg_b * log(N_Co_total*t_epsilon_Bkg_b) - (N_Co_total*t_epsilon_Bkg_b);
     
     if(S1Likelihood){
	
	// loop over the bins, check the content and compute
	for(int s1_bin=0; s1_bin < DataBand->getNSlices(); s1_bin++){

	   double obs_S1content = DataBand->getContent(s1_bin, b);

	   if (obs_S1content == 0) continue;  		//do not compute and also do no cansider cases ns_fs + nb_fb < 0

	   double ns_fs = SignalBand_t->getContent(s1_bin,b) * sigma * run->getSignalHandler()->getXsecMultiplier();
           if(isSigAccFit) ns_fs = ns_fs * ( 1. + t_sig_acc * 0.12); // acc uncertainty

	   double nb_fb = N_b * BackgroundBand->getContent(s1_bin,b)/BackgroundBand->getBandContent(b);
	   if(debugLevel > 1) cout << "Band  " << b << " Slice  " << s1_bin << "  ns_fs  " << ns_fs << "   nb_fb  " << nb_fb << "   obs_S1content  "<< obs_S1content << endl;
	   if(ns_fs + nb_fb > 0 ) LL_band += obs_S1content * ( log(ns_fs + nb_fb) - log(N_s + N_b) );
	   else return VERY_SMALL;
	}

     }
     
/*	




      //double N_obs     = run->getBandContent(ALL_BACKGROUNDS,b);
      double N_obs     = run->getBandContent(DM_CUT_DATA,b);
      if(debugLevel > 10) cout << "Band "<< b << " N_S  " << N_s << " sigma " << sigma << "  Histo content " << SignalBand_t->getBandContent(b) 
      			       << " N_b  " << N_b <<  " N_obs " << N_obs  << "  t_epsilon_Bkg_b " << t_epsilon_Bkg_b<< endl;
      if(debugLevel > 10) cout << "Band "<< b << " N_B  " << N_b << " N_bkg_tot " << N_bkg_tot << "  t_NormBkg " << t_NormBkg << "  SystBkgSigma  " << SystBkgSigma << endl; 
      
      // poisson on number of events in band	 
      //if(N_s + N_b > 0) LL_band += N_obs*log( N_s + N_b ) -1.*(N_s+N_b);
      if(N_s + N_b > 0) LL_band += PoissonDist::logPdf(N_obs , N_s + N_b);
      else return VERY_SMALL;
  
      // poisson constraint to stat. uncertainty NP (epsilon)
      //if(isStatBkgTValueFit) LL_band += log(TMath::Poisson( N_bkg_tot * epsilon_bkg_b, N_bkg_tot  * t_epsilon_Bkg_b ));
      //if(isStatBkgTValueFit) LL_band += PoissonDist::logPdf( N_bkg_tot * epsilon_bkg_b, N_bkg_tot  * t_epsilon_Bkg_b );
      if(isStatBkgTValueFit) LL_band += PoissonDist::logPdf( N_bkg_tot * epsilon_bkg_b * frac_Co_model, N_bkg_tot * frac_Co_model * t_epsilon_Bkg_b ); //TEST_ALE
      //if(isStatBkgTValueFit) LL_band +=  log(TMath::Poisson( N_bkg_tot  * t_epsilon_Bkg_b, N_bkg_tot * epsilon_bkg_b )); //TEST_ALE
     // if(isStatBkgTValueFit) LL_band +=  N_bkg_tot * epsilon_bkg_b * log(N_bkg_tot * t_epsilon_Bkg_b) -1.*(N_bkg_tot * t_epsilon_Bkg_b); //TEST_ALE

      //S1 Shape Likelihood (fs() and fb())
      if(S1Likelihood){
	//	cout << "NBands " << dataS1Bins.size() << " Nevent in band " << b << " = " << dataS1Bins[b].size() << endl;
*/
		/*
		//poisson style
		  for(int s =0; s < SignalBand_t->getNSlices(); s++){
		        double obs_s = Data->getContent(s,b);
        		double sig_s = SignalBand_t->getContent(s,b) * sigma * 1.E-6 ;                                    
        		double bkg_s = BackgroundBand->getContent(s,b);
        		if(sig_s + bkg_s > 0) LL_band += log( TMath::PoissonI(obs_s, sig_s+bkg_s));                      
       			 else LL_band+= VERY_SMALL;
     		}                                                                                                    
                          
        }  
 */
/*
	   //loop over data events for that band
	   for(unsigned int d=0; d < dataS1Bins[b].size(); d++) {
		int s1_bin = (dataS1Bins[b])[d];

//		cout << (dataS1Bins[b])[d] << endl;
		//NOTE: normalizzation is correct since SignalBand_t.and BackgroundBand are in number of events and not PDF normalized to 1.
		double ns_fs = SignalBand_t->getContent(s1_bin,b) * sigma *1.E-6;
                double nb_fb = N_b * BackgroundBand->getContent(s1_bin,b)/BackgroundBand->getBandContent(b);
	     	if(ns_fs + nb_fb > 0 ) LL_band += log( (ns_fs + nb_fb ) / (N_s + N_b) );
		else return VERY_SMALL;
		
//		cout << "Band " << b << " dataS1Bins " << s1_bin << " fs " << SignalBand_t->getContent(s1_bin,b) * sigma *1.E-6 << " fb " << BackgroundBand->getContent(s1_bin,b) <<
  //                      "  L  " << (SignalBand_t->getContent(s1_bin,b) * sigma *1.E-6 + BackgroundBand->getContent(s1_bin,b)) / (N_s + BackgroundBand->getBandContent(b)) <<
//			"  LL " <<  log( (SignalBand_t->getContent(s1_bin,b) * sigma *1.E-6 + BackgroundBand->getContent(s1_bin,b)) / (N_s + BackgroundBand->getBandContent(b))) << endl;
//		cout << "\t\t H. rescaled print : fs " << SignalBand_t->getContent(s1_bin,b) * sigma *1.E-6 / N_s << "   fb " << BackgroundBand->getContent(s1_bin,b) / N_b 
//			<< "   ( ) " << log(BackgroundBand->getContent(s1_bin,b) / N_bkg_tot *406.) <<endl;



  	   }
      }

      if(debugLevel > 1) cout << "S1S2PL::computeTheLogLikelihood() - Info : Likelihood value for band "<< LL_band << endl;
  */

    
      LL+= LL_band;	
  } 

 // cout << "S1S2PL::computeTheLogLikelihood() - DEBUG : t_Leff=" << t_Leff << "  t_Qy=" << t_Qy  <<"  t_NormBkg " << t_NormBkg << "  sigma " << sigma << endl;
 
  if(t_Leff !=0 || t_Qy !=0)  delete SignalBand_t;

//   std::cout << std::fixed;
   std::cout << std::setprecision(10);
  if(debugLevel > 1) cout << "S1S2PL::computeTheLogLikelihood() - Info : Likelihood value "<< LL << endl;

  return LL ; 
}


void S1S2PL::setData(int dataType) {

	data_type = dataType;
	//check
	if(run->getS1S2Bands(dataType) == NULL) {
		cout <<"S1S2PL::setData - FATAL ERROR - the chosen data bands are empty. Quit. " <<endl;
		exit(100);
	}

	DataBand = run->getS1S2Bands(dataType);
}

void S1S2PL::generateAsimov(double mu_prime){
	run->generateAsimovData(mu_prime);
}

void S1S2PL::generateToyDataset(double seed, double mu_prime){
	run->generateData(seed, mu_prime);
}


bool S1S2PL::printIt(int level){
/*  if(level<1) return true;
  printFlagsAndParameters();

  if(level<2) return true;
  cout<<endl<<"Band count"<<endl<<" Band ";
  for(int s=0;s<=N_SHAPES;s++) {
    cout<<rightJustify(S1S2Data::getTypeName(S1S2Data::mapS1Shape(s)),15);
  }
  cout<<endl;  
  for(int b=0;b<nBands;b++){
    cout<<setw(4)<<b;
    for(int s=0;s<N_SHAPES;s++) cout<<formatF(bandContents[s][b],15,4);
    cout<<endl;
  }

  if(level<3) return true;
  for(int b=0;b<nBands-1;b++) {
    printVector(*dataS1Bins[b],"Data S1 bins in Band "+formatI(b));
  } 

  if(level<4) return true;
  cout<<endl<<"Distributions of backgrounds"<<endl<<" Bin ";
  for(int s=0;s<N_SHAPES;s++) {
    cout<<rightJustify(S1S2Data::getTypeName(S1S2Data::mapS1Shape(s)),15);
  }
  cout<<endl;  
  for(int p=0;p<N_PE_POINTS;p++){
    cout<<setw(3)<<p;
    for(int s=0;s<N_SHAPES;s++) cout<<formatF(dists[s][p],15,4);
    cout<<endl;
  }
  cout<<endl;
  */
  return true;
}

void S1S2PL::printFlagsAndParameters() {
/*  

  activateParameter(lkacc,isSigAccFit);ProfileLikelihood::printFlagsAndParameters();
  displayFlag("Vary systematic bkg.",6,isSystBkgTValueFit
             ,DEFAULT_PAR_SYST_BKG_FIT);
  displayFlag("Vary statistical bkg.",6,isStatBkgTValueFit
             ,DEFAULT_PAR_STAT_BKG_FIT);
  displayFlag("Vary LEff",6,isLEffTValueFit,DEFAULT_PAR_LEFF_TVALUE_FIT);
  displayFlag("S1 likelihood",6,S1Likelihood,DEFAULT_S1_LIKELIHOOD);
*/
}

double S1S2PL::nSignalPerCm2()          {return run->nSignalPerCm2(0.);}
double S1S2PL::nSignalPerCm2(double lt) {return run->nSignalPerCm2(lt);}


//----------- ModelImporter -----------//

ModelImporter::ModelImporter(TFile *file, TString name) : histo_name(name) 
{
	
	f = file;

	if(f == NULL) {
		cout << "ModelImporter - WARNING :  file not specified for " + histo_name << endl; 
		model = NULL;
	}
 	else  import();
}

ModelImporter::~ModelImporter(){}

void ModelImporter::import()
{

	if(f == NULL) {
		cout << "ModelImporter - WARNING : file not specified for " + histo_name << endl; 
		model = NULL;
	}
	else if( f->FindKey(histo_name) == NULL) {
		if(histo_name == "") cout << "ModelImporter::import() - WARNING : name of histogram not specified." << endl;
		else cout << "ModelImporter - WARNING : name " + histo_name + " not found in file " << f->GetName()<< endl; 
		model = NULL;
	}
	else 	{ model = (TH2F*) f->Get(histo_name);  }//cout << "ModelImporter - " + histo_name + " Imported from file"<< endl;} 
}


void ModelImporter::setHistoName(TString name){
	histo_name = name;
	import();
}


void ModelImporter::setFile(TFile *file) { f = file;}


void ModelImporter::rebinX(int Rebin)  { //cout << "bins before " << model->GetNbinsX() << endl;
	 if(model == NULL) { cout <<"ModelImporter::rebinX  - WARNING : Cannot rebin, histo empty" << endl;}
	 else model->RebinX(Rebin);  //cout << "bins after " << model->GetNbinsX() << " Model Name " << model->GetName() << endl;
}

S1S2Bands * ModelImporter::fillBands(S1S2Bands *reference_band, bool isLog) 
{
	//-------------------------------   Important note	-------------------------------------------//
	// IMPORTANT::: keep the same TH2F binning for signal and bkg.
	// An approximation is made here: TH2F bins which falls only partially in a band will be 
	// assigned to the band depending on where is the bin center, this will slightly deform the bands, 
	// however this deformation will be the same for signal and bkg if the binning of histos are equal.
	//
	// Also, typical size of a "narrow" band in log(s2/s1) ~ 0.02  -> max TH2F bin width ~ 0.005
	//       typical size of a "narrow" band in  s2/s1 ~ 1  -> max TH2F bin width 0.25
	//
	// IMPORTANT:: If isLog is true then the histo is supposed to be in log(s2/s1) VS s1, otherwise in S2 VS S1
	//-----------------------------------------------------------------------------------------------------------//


	//Sanity check
	if(reference_band == NULL) {
		// if you are calling this function in this case there is something totally wrong in your code and needs to stop. 
		cout << "ModelImporter::fillBands() - Fatal Error, reference bands are not initialized. Quit." << endl;
		exit(100); 
	}

	// bands that will be filled, the structure is inherithed from "reference_band"
	// NOTE: copy flag = "false" means to initialize the band content to zero.
        S1S2Bands *Bands = new S1S2Bands( histo_name.Data(), reference_band,false);

	if(model == NULL) {
		cout << "ModelImporter::fillBands() - Error, TH2F not specified for " + histo_name <<  ".  Returning empty bands." << endl;

		return Bands;
	}
	
 	// Sanity check (assuming TH2F have all equal bins), Note: getMinimumBandSize(true) return minimum band size in log(s2/s1)
	if(isLog && Bands->getMinimumBandSize(isLog) < model->GetYaxis()->GetBinWidth(1)) {
		cout << "ModelImporter::fillBands() - ERROR : TH2F width exceding band size for " + histo_name << ". Returning empty bands" << endl;
		return Bands;
	}
	
	
	double counter_content=0.;

	// Loop over S1 slices
	for(int slice = 0; slice < Bands->getNSlices() ; slice++)
	{
		double s1 = Bands->getS1Center(slice);		// s1 value of the current slice
		int bin_x = model->GetXaxis()->FindBin(s1);	
		
		//sanity checks
		if( fabs(model->GetXaxis()->GetBinLowEdge(bin_x) - Bands->getS1LowerEdge(slice)) > 1.E-10
		    || fabs(model->GetXaxis()->GetBinUpEdge(bin_x) - Bands->getS1UpperEdge(slice)) > 1.E-10 ) 
		{
		  cout << "ModelImporter::fillBands() - ERROR : TH2F has not the standard s1 binning for " + histo_name << ". Returning empty bands" << endl;
		  cout << model->GetXaxis()->GetBinLowEdge(bin_x) << " ?= " << Bands->getS1LowerEdge(slice) << " | " << model->GetXaxis()->GetBinUpEdge(bin_x) << " ?= " <<  Bands->getS1UpperEdge(slice) << endl;
		  break;
		}	
		
		//Filling bands for the current s1 slice
		for( int bin_y =0; bin_y < model->GetNbinsY() ; bin_y++) 
		{
		  //get TH2F Y center 
		  double log_s2_over_s1 = model->GetYaxis()->GetBinCenter(bin_y);

		  //retrive the corresponding band bin number, note: bands are stored in s2/s1		  
		  int band = Bands->getSlice(slice)->getS2overS1Bins()->getBin(pow(10.,log_s2_over_s1));

		  if(isLog == false) band = Bands->getSlice(slice)->getS2overS1Bins()->getBin(log_s2_over_s1 / s1); // in this case the histo is assumed to be in s2 VS s1

		  if(band < 0.) {
		    if(reference_band->getDebugLevel() > 3) {
		       if (model->GetBinContent(bin_x, bin_y) > 0. && isLog  ) cout << "ModelImporter::fillBands() - WARNING : for histo " << histo_name << " bin containing " << model->GetBinContent(bin_x, bin_y) << " events is out of possible bands, S1= " << s1 <<  " and log(s2/s1)= "<< log_s2_over_s1 << ". Are you sure the Y axis is UNFLATTENED? "<< endl;
		       if (model->GetBinContent(bin_x, bin_y) > 0. && !isLog ) cout << "ModelImporter::fillBands() - WARNING : for histo " << histo_name << " bin containing " << model->GetBinContent(bin_x, bin_y) << " events is out of possible bands, S1= " << s1 <<  " and s2 = "<< log_s2_over_s1 <<  endl;
		     }


			 continue; // this bin is out of the possible bands, which will happen for some bins in almost all histos.
		  }
		  counter_content += model->GetBinContent(bin_x, bin_y);
		  Bands->fill(slice, band, model->GetBinContent(bin_x, bin_y));

		}

	}

                 if(reference_band->getDebugLevel() > 1)  cout << "ModelImporter::fillBands() - INFO :: The integral of bins included in bands of  " << histo_name << "  is  " << counter_content << endl; 
	return Bands;
}

TH2F ModelImporter::getModel(){ 
	if(model == NULL ) {
		cout << "ModelImporter::getModel()  - FATAL ERROR : histogram undefined. XEPHYR will probably crash..." << endl;  //error handling not defined FIXME
		return TH2F("EMPTY","",10,0,10,10,0,10); 
	}
	
	return *model;
}




//------------ IMPORTED ElectronBackground ----------------//
ImportedElectronBackground::~ImportedElectronBackground(){}

ImportedElectronBackground::ImportedElectronBackground() 
                            : ElectronBackgroundModel() {}

ImportedElectronBackground::ImportedElectronBackground(XeRun* r)
 :  ElectronBackgroundModel("Imported model for run " + r->getName(),r,PL_ANALYSIS),
    gaussBkg(r->getDataFile(), r->getName() + "_ERgauss"), 
    leakageBkg(r->getDataFile(),r->getName() + "_ERleak"){
 
 //Note: there are a couple of variables of ElectronBackgroundModel wich are set to 
 //default values by its constructor, these variable are unimportant for this class.


}

S1S2Bands* ImportedElectronBackground::computeBands(){

  // load REFERENCE bands
  S1S2Bands* r_bands = run->getReferenceBands();  

  // Sanity checks
  if(r_bands==NULL) {
    cout<<"ImportedElectronBackground::computeBands - Error - Reference bands are empty"<<endl
	<<"Run: " << run->getName() << " may not be initialize, make sure that you executed XeRun::fillDataBands(). "<<endl;
  }
  
  return gaussBkg.fillBands(r_bands, false);

}

S1S2Bands* ImportedElectronBackground::computeAnomalousBands(){

  // load REFERENCE bands
  S1S2Bands*  r_bands = run->getReferenceBands();  

  // Sanity checks
  if(r_bands==NULL) {
    cout<<"ImportedElectronBackground::computeBands - Error - Reference bands are empty"<<endl
	<<"Run: " << run->getName() << " may not be initialize, make sure that you executed XeRun::fillDataBands(). "<<endl;
  }
  
  return leakageBkg.fillBands(r_bands,false);

}

bool ImportedElectronBackground::printIt(int level){
  ElectronBackgroundModel::printIt(level);
  return true;
}

bool ImportedElectronBackground::isRunCompatible(int ) {return true;}


//----------- Imported Neutron model  -------------//

ImportedNeutronBackgroundModel::~ImportedNeutronBackgroundModel(){}

ImportedNeutronBackgroundModel::ImportedNeutronBackgroundModel() 
  : NeutronBackgroundModel() {}

ImportedNeutronBackgroundModel::ImportedNeutronBackgroundModel(XeRun* r)
 :  NeutronBackgroundModel("Imported for " + r->getName(),r,PL_ANALYSIS) , 
    neutronBkg(r->getDataFile(),r->getName() +  "_Neutron") { } 

bool ImportedNeutronBackgroundModel::isRunCompatible(int ) {return true;}


double ImportedNeutronBackgroundModel::expectedBelow(double pe){
  return run->getS1S2Bands(NR_BACKGROUND)->getContentBelowS1(pe);
}


S1S2Bands* ImportedNeutronBackgroundModel::computeBands(){

  // load REFERENCE bands
  S1S2Bands*  r_bands = run->getReferenceBands();  

  // Sanity checks
  if(r_bands==NULL) {
    cout<<"ImportedNeutronBackgroundModel::computeBands - Error - Reference bands are empty"<<endl
	<<"Run: " << run->getName() << " may not be initialize, make sure that you executed XeRun::fillDataBands(). "<<endl;
  }
  
  return neutronBkg.fillBands(r_bands, false);

}


//-------------- Signal Model -------------------//
SignalModel::SignalModel(XeRun *run) : XeObject()
{
  setName( "Signal Model for " + run->getName() + " and mass " + formatI(run->getWimpMass(), 4) + " GeV");

  currentRun 		=   run; 			 //Set the run
  mass 			=   run->getWimpMass();	 	 //current mass point 
  Default_XsectionNorm 	=   1.E-39;			 //cross section with which the signal events are computed by default.
  Xsecion_multiplier 	=   UNDEFINED;	         //Modifier of cross section for Minuit
  Xsection 		=   Default_XsectionNorm * Xsecion_multiplier; //reference cross section

  printlevel = 0;

  normOne = false;
  
} 

S1S2Bands SignalModel::computeSignalBands(double T_Leff, double T_Qy){

	// bands that will be filled with structure inherithed from "reference_band"
	// NOTE: copy flag = "false" means to initialize the band content to zero.
        S1S2Bands clone_bands(getName() + " Leff " + formatI(T_Leff,3) + " QY " + formatI(T_Qy,3) ,
				 currentRun->getReferenceBands(),false);

	return clone_bands; // To be modified. Ale. FIXME
}


double SignalModel::getXsecMultiplier(){
  return Xsecion_multiplier; 
}

void   SignalModel::setXsecMultiplier(double factor){

	if(normOne == true && Xsecion_multiplier != UNDEFINED) 
		cout << "SignalModel::setXsecMultiplier  - WARNING : you are re-setting a XsecMultiplier when was fixed to scale signal to 1. Ignore this warning if in combination mode." << endl;

	Xsecion_multiplier = factor;
}

bool SignalModel::setMass(double mass_point){

  setName( "Signal Model for " + currentRun->getName() + " and mass " + formatI(mass_point, 4) + " GeV");
  mass = mass_point;
  //Xsecion_modifier = getXsecModifier(mass_point); //FIXME this function would retrieve factor as a function of mass, but not ready yet
  Xsection = Default_XsectionNorm * Xsecion_multiplier;

  return true; // in this case is trivial and is always true, but for daughter may not be the case.

}

S1S2Bands SignalModel::computeReferenceBands() {

	if(printlevel > 0) cout << "SignalModel::computeReferenceBands()  - INFO : computing bands for " 				+ getName() << endl;

	// retrieving some variables.
	int nBands   		= currentRun->getNbands();
	int nSlices  		= currentRun->getNslices();
	double AnalysisS1Min	= currentRun->getAnalysisS1Min();
	double AnalysisS1Max	= currentRun->getAnalysisS1Max();

	//initializing bands	
    	S1S2Bands reference(currentRun->getNameSpace()+"Reference bands for mass " + formatI(mass, 4) ,
				currentRun ,nBands,nSlices,AnalysisS1Min,AnalysisS1Max);
	
	//safety check
	if(currentRun->getS1S2Data(AM_BE_CUT_DATA) == NULL) {
		cout<<"SignalModel::computeReferenceBands()  - ERROR : S1S2Data for AM_BE_CUT_DATA are not defined, returning empty reference bands!!!!! "<< endl;
		//HERE SHOULD THROW A FATAL ERROR, but error handling not yet implemented... FIXME
		return 	reference;
	}


	//computing reference bands based on AmBe CUT S1S2Data
	reference.compute(currentRun->getS1S2Data(AM_BE_CUT_DATA));	
	
	return reference;
}
	

SignalModel1D::SignalModel1D(XeRun *run) :SignalModel(run) {

    Defaultsignal = (TH1F*) run->getSignalFile()->Get("name");

    if(Defaultsignal==NULL) { cout<<"SignalModel1D:: no signal has been set" << endl;}

    PsigmaSignal = NULL;
    MsigmaSignal = NULL;

}


void SignalModel1D::normalizeToOne(){


	normOne = true;
	double integral = Defaultsignal->Integral();
	setXsecMultiplier(1./integral);


}

void SignalModel1D::setSignalName(TString nameHisto){

    Defaultsignal = (TH1F*) currentRun->getSignalFile()->Get(nameHisto);

    if(Defaultsignal==NULL) { cout<<"SignalModel1D:: no signal has been set, name "<< nameHisto << " NOT FOUND"  << endl; }
    else {cout << "SignalModel1D:: Signal has been set to model  -  " << nameHisto << endl; } 

}

void SignalModel1D::setLeffSysNames(TString nameHisto_psigma, TString nameHisto_msigma){

	PsigmaSignal = (TH1F*) currentRun->getSignalFile()->Get(nameHisto_psigma);

	MsigmaSignal = (TH1F*) currentRun->getSignalFile()->Get(nameHisto_msigma);

         
    if(PsigmaSignal==NULL || MsigmaSignal == NULL) { cout<<"SignalModel1D:: no Leff uncertainty  has been set, name "<< nameHisto_psigma << " NOT FOUND"  << endl; }
    else {cout << "SignalModel1D:: Leff Sys Signal has been set to model  -  " << nameHisto_psigma << endl; } 

}

double SignalModel1D::getLeffmultiplier(double t_val, int bin_number){
	
	double scale_multiplier = 1.;
	if(PsigmaSignal == NULL && MsigmaSignal == NULL) return scale_multiplier;
	
	double psigma = PsigmaSignal->GetBinContent(bin_number);
	double median = Defaultsignal->GetBinContent(bin_number);
        double msigma = MsigmaSignal->GetBinContent(bin_number);

	if(t_val > 0.) scale_multiplier = t_val * psigma / median  + 1. - t_val ;
	if(t_val < 0.) scale_multiplier =  1. + t_val -1.*(t_val * msigma / median) ;

	return scale_multiplier;

}
//------------------------ MCSignalModel----------------------------//

MCSignalModel::MCSignalModel(XeRun *run, int Rebin,TString suffix) : SignalModel(run) , 
			signalInput(run->getSignalFile(), computeName(0.,0.,1., suffix)) {
	rebin = Rebin; 
	signalInput.rebinX(rebin);

	suffisso = suffix;
}

bool MCSignalModel::setMass(double mass_input){
	
  	setName( "Signal Model for " + currentRun->getName() + " and mass " + formatI(mass_input, 4) + " GeV");

	//check consistency.
	if(currentRun->getWimpMass() != mass_input) cout << "MCSignalModel::setMass() - WARNING : mass set to value different from XeRun. " << endl;

	mass = mass_input;

 	signalInput.setHistoName(computeName(0.,0.,1.,suffisso)); 

	return true; // ModelImporter is supposed to return a bool after import... FIXME

}

TString MCSignalModel::computeName( double tLeff, double tQy, double LCE, TString suffix) {

  char bleff[100];
  // to cope with stupid naming scheme (lol)
  if(tLeff == 0.) tLeff = 0.0010; 
  if(tQy   == 0.) tQy   = 0.0010; 

  //Note: "mass" is a member of SignalModel	
  sprintf(bleff,"Mass%.0f_Leff%.4f_Qy%.4f_LCE%.4f",mass,tLeff,tQy,LCE);
  TString name("h2dc_");
  name.Append(bleff);
  name.Append(suffix);

  return name;

}


S1S2Bands MCSignalModel::computeSignalBands(double T_Leff, double T_Qy){
  

 
  // import the new histo with the specified Leff and Qy, the mass of the wimp has to
  // be defined previously, and is defined by default in the SignalModel costructor. 
  signalInput.setHistoName(computeName(T_Leff, T_Qy, 1., suffisso)); 

  //Rebinning the Model after loading a new histo
  signalInput.rebinX(rebin);

  cout << " MCSignalModel::computeSignalBands  Nevents signal " << signalInput.getModel().Integral() << endl;

  //Note parameters: isLog set to false
  S1S2Bands *returningBands = signalInput.fillBands(currentRun->getReferenceBands(), false);  //returning pointer to local variable should be avoided FIXME

  return *returningBands;



}


S1S2Bands MCSignalModel::computeReferenceBands(){

  
  // import the histo with default Leff and Qy, the mass of the wimp has to
  // be defined previously. 
  signalInput.setHistoName(computeName(0.,0.,1.,suffisso)); 

  
  //Rebinning the Model after loading a new histo
  signalInput.rebinX(rebin);

  // retrieving some variables.
  int nBands   		= currentRun->getNbands();
  int nSlices  		= currentRun->getNslices();
  double AnalysisS1Min	= currentRun->getAnalysisS1Min();
  double AnalysisS1Max	= currentRun->getAnalysisS1Max();

  //initializing bands	
  S1S2Bands reference(currentRun->getNameSpace()+"Reference bands for mass " + formatI(mass, 4) ,
				currentRun ,nBands,nSlices,AnalysisS1Min,AnalysisS1Max);

  //compute the bands from imported TH2F histo
  reference.compute(signalInput.getModel(), false);

  cout << "MCSignalModel::computeReferenceBands() INFO : " << reference.getName() << " N signal events rescaled to " << reference.getTotalContent() * getXsecMultiplier() << endl;

  return reference;


}

TH2F MCSignalModel::getHisto(){

  // import the histo with default Leff and Qy, the mass of the wimp has to
  // be defined previously. 
  signalInput.setHistoName(computeName(0.,0.,1.,suffisso)); 

  
  //Rebinning the Model after loading a new histo
  signalInput.rebinX(rebin);

  return signalInput.getModel();

}


