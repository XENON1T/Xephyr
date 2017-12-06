#include "XeStat.h"



int    XeStat::analysisMode = NO_ANALYSIS;

       XeStat::XeStat(TString nam) : errorHandler("Stat")    {
           setName(nam);
           setExperiment(1);
        }

void   XeStat::setAnalysisMode(int m){analysisMode=m;}
bool   XeStat::isCutsBased()         {return analysisMode==CUTS_ANALYSIS;}
bool   XeStat::isPL()                {return analysisMode==PL_ANALYSIS;}
int    XeStat::getAnalysisMode()     {return analysisMode;}
TString XeStat::getAnalysisModeName() {return getAnalysisModeName(analysisMode);}

TString XeStat::getSystematicModeName(int i){
  switch(i){
    case ONE_SIGMA_BELOW   : return "1 sigma below";
    case CENTRAL           : return "Central Value";
    case ONE_SIGMA_ABOVE   : return "1 sigma above";
  }
  return UNDEFINED_STRING;
}

bool XeStat::isAnalysisDefined(bool verbose){
  bool undef=analysisMode==NO_ANALYSIS;
  if(undef && verbose){
    cout<<endl<<" You forgot to set the Analysis Mode "<<endl;
  }
  return !undef;
}



TString  XeStat::getSigmaUnitName(int unit) {
  switch(unit) {
    case SIGMA_UNIT  : return "cross section";
    case EVENT_UNIT  : return "number of events";
  }
  return UNDEFINED_STRING;
}

TString  XeStat::getSigmaModeName(int mode) {
  switch(mode) {
    case ESTIMATED   : return "Estimated";
    case UPPER_LIMIT : return "Upper Limit";
  }
  return UNDEFINED_STRING;
}

TString XeStat::getSigmaLabel(int unit){
  switch(unit) {
    case SIGMA_UNIT  : return LABEL_SIGMA;
    case EVENT_UNIT  : return LABEL_EVT;
  }
  return UNDEFINED_STRING;
}

TString XeStat::getUpperSigmaLabel(int unit){
  switch(unit) {
    case SIGMA_UNIT  : return LABEL_SIGMA_MAX;
    case EVENT_UNIT  : return LABEL_EVT_MAX;
  }
  return UNDEFINED_STRING;
}

int XeStat::getSigmaLinLog(int unit) {return(unit==SIGMA_UNIT? LOG:LINEAR);}

TString  XeStat::getAnalysisModeName(int m) {
  switch(m) {
    case NO_ANALYSIS   : return "Undefined";
    case PL_ANALYSIS   : return "Profile likelihood";
    case CUTS_ANALYSIS : return "Cuts based";
  }
  return UNDEFINED_STRING;
}

bool XeStat::checkAnalysisMode(TString name, int requestedAnalysisMode){
  if(   requestedAnalysisMode==analysisMode 
     || requestedAnalysisMode==NO_ANALYSIS) return true;

  if(analysisMode==NO_ANALYSIS) {
    cout<<"On request of "<<name<<", setting analysis mode to "
        <<getAnalysisModeName(requestedAnalysisMode)<<endl;
    setAnalysisMode(requestedAnalysisMode);
    return true;
  }

  cout<<name<<" is not comptabile with the current analysis mode :"
       <<getAnalysisModeName()<<endl;
  return false;
}




//--------------- Parameters of a likelihood computations ------------

LKParameter::~LKParameter(){}

LKParameter::LKParameter(TString nam) : XeStat(nam) { t0=0.;}
LKParameter::LKParameter(int i, int typ,TString nam,double initV,double st
                    ,double mi,double ma) : XeStat(nam) {
  if(typ<0 || typ>=N_PARAMETER_TYPES){
    cout<<"Unknown LK parameter type "<<typ<<endl;
    return ;
  }

  Debug("LKParameter::LKParameter()  " + nam , " Enter");

  t0 = 0.;
  id=i;
  setType(typ);
  setInitialValue(initV);
  setStep(st);
  setMinuitUnit();
  setMinimum(mi);
  setMaximum(ma);
  initialize();
  setCommon(false);
  setCombinedMode(false);

  initialSigma = 1.; // tipicly is 1. for a standard NP, but of course not for a Stat NP.

  Debug("LKParameter::LKParameter()  " + nam , " Exit" );
}

void   LKParameter::setId(int i)              {id=i;}
void   LKParameter::setCommon(bool c)         {common=c;}
void   LKParameter::setType(int typ)          {type=typ;}
void   LKParameter::setInitialValue(double v) {initialValue=v;}
void   LKParameter::setCurrentValue(double v) {currentValue=v;}
void   LKParameter::setStep(double st)        {step=st;}
void   LKParameter::setSigma(double si)       {sigma=si;}
void   LKParameter::setMinimum(double mi)     {minimum=mi;}
void   LKParameter::setMaximum(double ma)     {maximum=ma;}
void   LKParameter::setCombinedMode(bool m)   {combinedMode=m;}

int    LKParameter::getType()                 {return type;}
int    LKParameter::getId()                   {return id;}
bool   LKParameter::isCommon()     {return common && experiment!=UNDEFINED_INT;}
bool   LKParameter::inCombinedMode()          {return combinedMode;}
double LKParameter::getInitialValue()         {return initialValue;}
double LKParameter::getCurrentValue()         {return currentValue;}
double LKParameter::getStep()                 {return step;}
double LKParameter::getMinimum()              {return minimum;}
double LKParameter::getMaximum()              {return maximum;}
double LKParameter::getSigma()                {return sigma;}

void   LKParameter::setCurrentValueInMinuitUnits(double v) {
  currentValue=v*MinuitUnit;
}

void   LKParameter::setSigmaInMinuitUnits(double si) {
  sigma=si*MinuitUnit;
}

void   LKParameter::setMinuitUnit(double s){
  Debug("Minuit Scaling of " + getName(), TString::Itoa(s,10));
  MinuitUnit=s;
}
double LKParameter::getCurrentValueInMinuitUnits() {
  return currentValue/MinuitUnit;
}

double LKParameter::getInitialValueInMinuitUnits() {
  return initialValue/MinuitUnit;
}

double LKParameter::getMinimumInMinuitUnits()     {
  return minimum/MinuitUnit;
}

double LKParameter::getStepInMinuitUnits()        {
  return step/MinuitUnit;
}

double LKParameter::getMaximumInMinuitUnits()     {
  return maximum/MinuitUnit;
}

void LKParameter::initialize() {
  setCurrentValue(initialValue); 
  setSigma(UNDEFINED);
}

double LKParameter::getLLGausConstraint(){
	//no constraint if free
	if(getType() == FREE_PARAMETER || getType() == FIXED_PARAMETER) return 0.;

	double t = getCurrentValue();
	return ( -1 * (t - t0)  * (t - t0)  /2. );  // t0 is the current measure 
}


bool LKParameter::compares(LKParameter* par,bool prnt){
  if(id!=par->getId()) {
    if(prnt) cout<<"Mismatch in id "<<id<<" vs "<<par->getId()
                  <<" for "<<getName()<<endl;
    return false;
  }
  if(type!=par->getType()) {
    if(prnt) cout<<"Mismatch in type "<<type<<" vs "<<par->getType()
                  <<" for "<<getName()<<endl;
    return false;
  }
  if(initialValue!=par->getInitialValue()) {
    if(prnt) cout<<"Mismatch in initial value "<<initialValue<<" vs "
                  <<par->getInitialValue()<<" for "<<getName()<<endl;
    return false;
  }
  if(step!=par->getStep()) {
    if(prnt) cout<<"Mismatch in step "<<step<<" vs "<<par->getStep()
                  <<" for "<<getName()<<endl;
    return false;
  }
  if(minimum!=par->getMinimum()) {
    if(prnt) cout<<"Mismatch in minimum "<<minimum<<" vs "
                  <<par->getMinimum()<<" for "<<getName()<<endl;
    return false;
  }
  if(maximum!=par->getMaximum()) {
    if(prnt) cout<<"Mismatch in maximum "<<maximum<<" vs "
                  <<par->getMaximum()<<" for "<<getName()<<endl;
    return false;
  }
  return true;
}

bool LKParameter::isActive()       {
  return    type==PARAMETER_OF_INTEREST
         || type==NUISANCE_PARAMETER
         || type==FROZEN_PARAMETER
         || type==FREE_PARAMETER;
}

bool LKParameter::isOfInterest()    {
  return    type==PARAMETER_OF_INTEREST
         || type==FROZEN_PARAMETER;
}

TString LKParameter::getTypeName()        {return getTypeName(type);}
TString LKParameter::getTypeName(int type){
  switch(type) {
    case FROZEN_PARAMETER      : return "Frozen     ";
    case PARAMETER_OF_INTEREST : return "Of interest";
    case NUISANCE_PARAMETER    : return "Nuisance   ";
    case FIXED_PARAMETER       : return "Fixed      ";
    case FREE_PARAMETER        : return "Free       ";
  } 
  return UNDEFINED_STRING;
}

void LKParameter::freeze(bool doFreeze){
  if(doFreeze) {
    if(type==PARAMETER_OF_INTEREST) type=FROZEN_PARAMETER;
  }
  else {
    if(type==FROZEN_PARAMETER) type=PARAMETER_OF_INTEREST;
  }
}

void LKParameter::printHeader(){
  cout<<formatI(id,3).Data()
      <<" "<<leftJustify(getName(),20).Data()
      <<" "<<getTypeName().Data() ; //<< endl;
}

void LKParameter::printInitial(){
  printHeader();
  cout<<formatG(initialValue,8,2).Data()
      <<formatG(step,8,2).Data()
      <<formatG(minimum,8,2).Data()
      <<formatG(maximum,8,2).Data()
      <<endl;
}

void LKParameter::printCurrent(bool withError){
  printHeader();
  cout<<" "<<formatRG(currentValue,8,2).Data();
  if(withError && (type==PARAMETER_OF_INTEREST || type==NUISANCE_PARAMETER || type==FREE_PARAMETER)) {
    if(sigma>0.) cout<<" +- "<<formatLG(sigma,8,2).Data();
    else         cout<<" +-     (undef)";
  }
  else           cout<<"                ";
  //cout<<" "<<formatRG(getCurrentValueInMinuitUnits(),8,2).Data()<<endl;
  cout<<" "<<endl;
}

SigmaParameter::~SigmaParameter(){}
SigmaParameter::SigmaParameter() : LKParameter(PAR_SIGMA,PARAMETER_OF_INTEREST
                ,"Sigma",0.,0.001,-50.,50.) {setCommon();}

TSystBkgParameter::~TSystBkgParameter(){}
TSystBkgParameter::TSystBkgParameter(int run) : LKParameter(PAR_SYST_BKG_TVALUE
        ,NUISANCE_PARAMETER ,getTheName(run),0.,0.01,-5.,5) {}

TString TSystBkgParameter::getTheName(int run){

  return "Syst. bkgd t-Value run"+formatI(run,2) ;
}


//----------TGauss--------///
TGaussParameter::~TGaussParameter(){}

// this is a t_value
TGaussParameter::TGaussParameter(int b,int run) : LKParameter(PAR_GAUSS_TVALUE+b
        ,NUISANCE_PARAMETER ,getTheName(b, run),0.,0.01,-5.,5.) { uncertainty = 1. ;}

TString TGaussParameter::getTheName(int b, int run){
  return "Gauss NP "+formatI(b,3) + " run"+formatI(run,2) ;
}


//------------------------//

TStatBkgParameter::~TStatBkgParameter(){}

TStatBkgParameter::TStatBkgParameter() : LKParameter("FuckYouGiveMeAName"){ stat_error = 0.;}

// this is not a t_value
TStatBkgParameter::TStatBkgParameter(int b,int run) : LKParameter(PAR_STAT_BKG_TVALUE+b
        ,NUISANCE_PARAMETER ,getTheName(b, run),0.1,0.01,0.,1.) {band=b; stat_error =0.;}

TString TStatBkgParameter::getTheName(int b, int run){
  return "Stat .bkgd for band"+formatI(b,3) + " run"+formatI(run,2) ;
}
 

void TStatBkgParameter::setStatError(double err){
    stat_error = err;
    setSigma(err);  // so I can use get sigma

    initialSigma = err;
}
void TStatBkgParameter::printCurrent(bool withError){
  printHeader();
  cout<<" "<<formatRG(currentValue,8,2);
  if(withError && (type==PARAMETER_OF_INTEREST || type==NUISANCE_PARAMETER)) {
    if(sigma>0.) cout<<" +- "<<formatLG(sigma,8,2);
    else         cout<<" +-     (undef)";
  }
  else           cout<<"                ";
  if(stat_error > 0.) cout<<" "<<formatRG(getCurrentValueInMinuitUnits(),8,2)<< "    Pull "<< (currentValue - initialValue) / stat_error << endl;
  else cout<<" "<<formatRG(getCurrentValueInMinuitUnits(),8,2)<< "    Pull " << endl;
}


TLEffParameter::~TLEffParameter(){}
TLEffParameter::TLEffParameter() : LKParameter(PAR_LEFF_TVALUE
              ,NUISANCE_PARAMETER,"LEff t-Value",0.,0.01,-1.,1.) {setCommon();}

TQyParameter::~TQyParameter(){}
TQyParameter::TQyParameter() : LKParameter(PAR_QY_TVALUE
              ,NUISANCE_PARAMETER,"Qy t-Value",0.,0.01,-1.,1) {setCommon();} //step for minimizzation set to 1 CHECKME FIXME



TEfficiencyParameter::~TEfficiencyParameter(){}
TEfficiencyParameter::TEfficiencyParameter() : LKParameter(PAR_EFFICIENCY_TVALUE
              , NUISANCE_PARAMETER,"Efficiency t-Value",0.,1.,-2.5,2.5) {}  //step for minimizzation set to 1 CHECKME FIXME


//--------------- virtual class for Likelihood computation --------

Likelihood::~Likelihood()                     {clear();}
Likelihood::Likelihood(TString nam) : XeStat(nam) {
  setup();
}

double Likelihood::getSigmaHat()  {return sigmaHat;}

double Likelihood::getLogD()      {return LogD;}


void Likelihood::setup(){
  currentId=0;
  combinedMode=false;
  forcedNParametersOfInterest=UNDEFINED_INT;
  seed =0;
  sigmaHat           = UNDEFINED;
  LogD               = UNDEFINED;
}

void Likelihood::clear(){
  TRAVERSE_PARAMETERS(it) {
    LKParameter *p=it->second; 
    deleteWithPointer(p);
  }
  clearTheParameters();
}

void Likelihood::clearTheParameters(){
  parameters.clear();
  MinuitParameters.clear();
}

void Likelihood::printInitialHeader(){
  cout<<" Id Name                  Type         Exp    Initial  Step    Min     Max"
      <<endl;
}

void Likelihood::printCurrentHeader(){
  cout<<" Id Name                 Type         Current Value               Uncertainty"
      <<endl;
}

void Likelihood::removeParameter(LKParameter *par, bool tolerant){
  removeParameter(par->getId(),tolerant);
}

void Likelihood::removeParameter(int id, bool tolerant){
  ParameterIterator it=parameters.find(id);
  if(it==parameters.end()){
    if(!tolerant) cout<<"Can't remove non exisiting parameter "<<id<<endl;
    return;
  }
  parameters.erase(it);
}

void Likelihood::addParameterTolerant(LKParameter* param){
  Info("Likelihood::Adding parameter ", param->getName());
  param->setExperiment(experiment);
  int id=param->getId();
  parameters[id]=param;
  currentId=max(id,currentId);
}

void Likelihood::activateParameter(LKParameter* param, bool active){
  if(active) addParameterTolerant(param);
  else       removeParameter(param,true);
}

void Likelihood::addParameter(LKParameter* param,int id){
  //adding parameter in case it possibly exist SAME, if exits wont add it
  // this is usefull for combination and parameter with known ID like sigma
  if(id==SAME) id=param->getId();

  //AUTO: case you don't know if parameter exist, if doesn't exist it assign a new ID
  else if(id==AUTO ) {
	  if( param->getId() == PAR_NOT_ASSIGNED) id  = ++currentId; 
  	  else                                    id  = param->getId();
  }
  cout<<"\tLikelihood::addParameter - Info : Adding parameter "<< 
	  param->getName() <<"  with ID " << id << "  to PL " << 
	  getName() << endl;

  //check if ID exist
  if(checkParameter(id,false)) parameters[id]=param;

  currentId=max(id,currentId);
  param->setId(id);
  param->setExperiment(experiment);

}

void Likelihood::replaceParameter(LKParameter* param){
  int id=param->getId();
  if(parameters[id]==NULL) {
    cout<<"  ==> Error, can't replace inexisting parameter"<<endl;
    return;
  }
  if(parameters[id]!=param) {
    delete parameters[id];
    parameters[id]=param;
     Info("Replacing parameter ", TString::Itoa(id,10) + " with existing common one");
  }
  else {
    Debug("replaceParameter" , "Keeping same parameter " + TString::Itoa(id,10));
  }
}

void Likelihood::addParameter(int id, int type,TString nam,double initialVal
                             ,double step,double mi, double ma) {
  if(id==AUTO) id=++currentId; 
  addParameter(new LKParameter(id,type,nam,initialVal,step,mi,ma));
}

bool Likelihood::parameterExists(int p) {
 return parameters.find(p)!=parameters.end();
}

bool Likelihood::checkParameter(int p,bool shouldExist) {
  bool exist=parameterExists(p);
  if(shouldExist){
    if(!exist){
      cout<<"Parameter "<<p<<" does not exist"<<endl;
      return false;
    }
  }
  else if(exist){
    cout<<"\t ==> Likelihood::checkParameter - Parameter with ID "<<p<<" already exists! Set to be corelated." << endl;
    return false;
  }
  return true;
}

void Likelihood::listParameters(){
  cout<<parameters.size()<<" parameters :";
  TRAVERSE_PARAMETERS(it) {
    LKParameter *p=it->second; 
    cout<<" "<<it->first<<"("<<p->getName()<<")";
  }
  cout<<endl;
}

void Likelihood::resetParameters(){
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    p->initialize(); 
  }
}

map<int,LKParameter*>* Likelihood::getParameters(){return & parameters;}

LKParameter* Likelihood::getParameter(int id){
  ParameterIterator it=parameters.find(id);
  if(it==parameters.end()) {
    cout<<"Can't find parameter "<<id<<endl;
    return NULL;
  }
  LKParameter *p=it->second; 
  return p;
}

double Likelihood::getParameterValue(int id){
  LKParameter* par=getParameter(id);
  if (par == NULL ) cout << "Parameter ID " << id << " Not defined --- returning " << UNDEFINED << endl;
  return par==NULL? UNDEFINED : par->getCurrentValue();
}

void Likelihood::setParameterValue(int id,double v){
  LKParameter* par=getParameter(id);
  if(par!=NULL) par->setCurrentValue(v);
}

void Likelihood::setParameterType(int id,int type){
  LKParameter* par=getParameter(id);
  if(par!=NULL) par->setType(type);
}

void Likelihood::ignoreParameter(int id){setParameterType(id,FIXED_PARAMETER);}

void Likelihood::setInitialValue(int id,double v){
  LKParameter* par=getParameter(id);
  if(par!=NULL) par->setInitialValue(v);
}

void Likelihood::printInitialParameters(){
  cout<<endl<<getName()<<" has "<<parameters.size()<<" parameters:"<<endl;
  printInitialHeader(); 
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    p->printInitial(); 
  }
  cout<<endl;
}

void Likelihood::printResultParameters(){
  printCurrentHeader(); 
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    p->printCurrent(true); 
  }
}

void Likelihood::printCurrentParameters(){
  printCurrentHeader(); 
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second;
    p->printCurrent(false); 
  }
}


TH1F Likelihood::getPullsHisto(){

   TString titleHisto = "Pulls for "+getName();
   TH1F pullsHisto("pullsHisto",titleHisto,getNTotalParameters(), 1, getNTotalParameters());

   int bin = 1;

   TRAVERSE_PARAMETERS(it) {
	LKParameter *p=it->second;

	double pull =  (p->getCurrentValue() - p->getInitialValue() ) / p->getInitialSigma();
	if(p->getType() == PARAMETER_OF_INTEREST) pullsHisto.SetBinContent(bin, p->getCurrentValue());
	else pullsHisto.SetBinContent(bin, pull);
	pullsHisto.SetBinError(bin, p->getSigma());
	TString paraName = p->getName();
	pullsHisto.GetXaxis()->SetBinLabel(bin, paraName.Data());
	
	bin++;	
   }	
   
   return pullsHisto;
    	
}

bool Likelihood::inCombinedMode()         {return combinedMode;}
int  Likelihood::getNMinuitParameters()   {return (int)MinuitParameters.size();}
int  Likelihood::getNTotalParameters()    {return (int)parameters.size();}
int  Likelihood::getNNuisanceParameters() {
  return getNParameters(NUISANCE_PARAMETER);
}

int  Likelihood::getNParameters(int t)  {
  int n=0;
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    if(p->getType()==t) n++;
  }
  return n;
}

void Likelihood::forceNParametersOfInterest(int nF){
  forcedNParametersOfInterest=nF;
}

int Likelihood::getNParametersForChi2(){
  if(forcedNParametersOfInterest>0) return forcedNParametersOfInterest;
  return nParametersOfInterest;
}

void Likelihood::setCombinedMode(bool m)  {
  combinedMode=m;
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    p->setCombinedMode(combinedMode);
  }
}

int  Likelihood::getNActiveParameters()  {
  int n=0;
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    if(p->isActive()) n++;
  }
  return n;
}

int  Likelihood::getNParametersOfInterest()  {
  int n=0;
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    if(p->isOfInterest()) n++;
  }
  return n;
}

bool Likelihood::checkConsistency(){
  nParametersOfInterest=getNParametersOfInterest();
  nActiveParameters=getNActiveParameters();
  nNuisanceParameters=getNNuisanceParameters();
  if(nParametersOfInterest==0){
    cout<<"No parameter of interest in "<<getName()<<endl;
    return false;
  }
  else if(nParametersOfInterest>5){
    cout<<"Too many parameters of interest in "<<getName()<<endl;
    return false;
  }
  return true;
}

int Likelihood::mapMinuitParameters(bool freeze){
  MinuitParameters.clear();
  TRAVERSE_PARAMETERS(it) { 
    LKParameter *p=it->second; 
    p->freeze(freeze); 
    int t=p->getType();
    if( t==NUISANCE_PARAMETER || t==PARAMETER_OF_INTEREST || t==FREE_PARAMETER) {
      MinuitParameters.push_back(p);
    }
  }
  int np=MinuitParameters.size();
  return np;
}

void Likelihood::setCurrentValues( const double *values,const double *errs) {
  int n=getNMinuitParameters();
  for(int p=0;p<n;p++) {
    MinuitParameters[p]->setCurrentValue(values[p]);
    if(errs!=NULL) MinuitParameters[p]->setSigma(errs[p]);
  }
}

void Likelihood::setCurrentValuesInMinuitUnits( const double *v,const double *e) {
  int n=getNMinuitParameters();
  for(int p=0;p<n;p++) {
    MinuitParameters[p]->setCurrentValueInMinuitUnits(v[p]);
    if(e!=NULL) MinuitParameters[p]->setSigmaInMinuitUnits(e[p]);
  }
}

//-- Only external static function could work!!! <-- ALE, THIS IS NOT TRUE! 
//look here: https://root.cern.ch/how-implement-mathematical-function-inside-framework
//Functor can call also methods of a class, clean this mess... FIXME.
static Likelihood*  currentLikelihood;

double LikelihoodEvaluate (const double * values) {
  currentLikelihood->setCurrentValuesInMinuitUnits(values); // write the current values to LKParameters


  if( errorHandler::globalPrintLevel < 1) { //Debug print
    cout<<"Evaluating likelihood "<<endl;
    currentLikelihood->printCurrentParameters();
  }
  double e= -1. * currentLikelihood->computeTheLogLikelihood();
  if(errorHandler::globalPrintLevel < 1) {
    cout<<"             .... result:"<<printTools::formatF(e,19,8)<<endl; 
  }
  return e;
}

double Likelihood::maximize(bool freezeParametersOfInterest){
  currentLikelihood=this;
  int np=mapMinuitParameters(freezeParametersOfInterest);

  nActiveParameters = getNActiveParameters();

  if(getPrintLevel() < 2) {
    cout<<"Finding maximum of "+getName()
        <<endl
        <<"Total of "<<nActiveParameters<<" active parameters, "
        <<doOrDont(freezeParametersOfInterest)
        <<" freeze parameters of interest, Minuit fits "<<np<<" param."
        <<endl;
    printCurrentParameters();
  }
  if(np==0){
    double e=computeTheLogLikelihood();
    if(getPrintLevel() < 2) {
      cout<<"Nothing to minimize, result:"<<formatF(e,19,8)<<endl; 
    }
    return e; 
  } 

 string m="Minuit2";
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(m,"Migrad"); // ALE_TEST  -- before was Migrad
  //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(m,"Simplex"); // ALE_TEST  -- before was Migrad
  if(min==NULL) {
    cout<<"Can't load "<<m<<"; is it part of your current ROOT installation?"
        <<endl;
    return UNDEFINED;
  }
   
  // set tolerance , etc...
  ROOT::Math::Functor f(&LikelihoodEvaluate,np);
  min->SetFunction(f);
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2       //TEST_ALE was 100000
  min->SetMaxIterations(100000);  // for GSL                     //was           10000
  min->SetTolerance(0.001); 					// was 0.01		
  min->SetPrintLevel(-1);   // quiet
  for(int i=0;i<np;i++){
    LKParameter* par=MinuitParameters[i];
   TString pnS= par->getName();
    string pn = pnS.Data(); 
    double v=par->getInitialValueInMinuitUnits();
    double s = par->getStepInMinuitUnits();
    //double s= par->getStepInMinuitUnits(); // don't ask why but the factor 100 is needed for better convergence
    double vmi=par->getMinimumInMinuitUnits();
    double vma=par->getMaximumInMinuitUnits();
    if(getPrintLevel() < 1) {
      double v0=par->getInitialValue();
      double s0=par->getStep();
      double vmi0=par->getMinimum();
      double vma0=par->getMaximum();
      cout<<"Setting parameter "<<(i+1)<<": "<<pn<<endl
      <<"      Initial value  :"<<v0<<" min:"<<vmi0<<" max:"<<vma0<<" step:"<<s0
      <<endl
      <<"      In Minuit units:"<< v<<" min:"<<vmi <<" max:"<<vma <<" step:"<<s 
      <<endl;
    }
    min->SetLimitedVariable(i,pn.c_str(),v,s,vmi,vma);
  } 
  // do the minimization
  min->Minimize(); 
  
  // write the post fit value to LKparameter 
  setCurrentValuesInMinuitUnits(min->X(),min->Errors());
  double ml= -1. * min->MinValue();
  if(getPrintLevel() < 2) {
    cout<<"ML "<<ml<<" achieved for "<<endl;
    printResultParameters();
  }

  //set sigma_hat to current and also logD
  if(!freezeParametersOfInterest) {
	sigmaHat = parameters[PAR_SIGMA]->getCurrentValue();
	LogD = ml;
   }
 
  delete min;
  return ml;
}



double Likelihood::maximizeNumerically(int numberOfToys, bool freezeParametersOfInterest){
  // Actually this function minimize the -log(L) for consistency with what was written before.
  
  double ML = VERY_LARGE; // maximum of the likelihood (we are minimizing actually)

  currentLikelihood=this; // this is needed to feed the input parameter into the likelihood, it can be done better.

  //retrieving number of parameters
  int np=mapMinuitParameters(freezeParametersOfInterest);

  nActiveParameters = getNActiveParameters();

  // Some printing
  if(getPrintLevel() < 1) {
    cout<<"Finding maximum NUMERICALLY of "+getName()
        <<endl
        <<"Total of "<<nActiveParameters<<" active parameters, "
        <<doOrDont(freezeParametersOfInterest)
        <<" freeze parameters of interest, Numerical fits of "<<np<<" param."
        <<endl;
    printCurrentParameters();
  }

  if(np==0){
    double e=computeTheLogLikelihood();
    if(getPrintLevel() <2) {
      cout<<"Nothing to maximize, result:"<<formatF(e,19,8)<<endl; 
    }
    return e; 
  } 

  //function that evaluates log(L) for the set of parameter specified.
  //for usage see: https://root.cern.ch/how-implement-mathematical-function-inside-framework
  ROOT::Math::Functor computeLL(&LikelihoodEvaluate,np);

  // initializzation of the vector of values of NP to feed to Functor
  double bestGuess_val_np[np] ; 	//initial values
  double min_val_np[np] ; 		//minimum of the range for random sampling
  double max_val_np[np] ;		//maximum of the range for random sampling 
  double temp_val_np[np] ;  		//temporary value
  double postFit_val_np[np] ;  		//post fit value
  double aux_postFit_val_np[np] ;  		//post fit value

  cout << "Toy generation settings : " << endl;
  for(int i=0; i < np; i++){
    LKParameter*   par  = MinuitParameters[i];
    bestGuess_val_np[i] = par->getInitialValueInMinuitUnits();
  
    if(par->getType() == PARAMETER_OF_INTEREST) {  // just around zero for the Xsec.
	min_val_np[i]       = -1.5; 
    	max_val_np[i]       =  1.;
    }	
    else if(par->getInitialValueInMinuitUnits()  == 0) {  // case of normal T-valued parameter variated of 2 sigma
	min_val_np[i]       = -2.; 
    	max_val_np[i]       =  2.;
    }
    else {  // case stat. parameter
	min_val_np[i]	   =  bestGuess_val_np[i] - par->getStepInMinuitUnits() * 20. ;	// minus 2 sigma;
	if(min_val_np[i] < 0.) min_val_np[i] = 0.;

	max_val_np[i]	   =   bestGuess_val_np[i] + par->getStepInMinuitUnits() * 20. ;// plus 2 sigma;
	if(max_val_np[i] >1.) max_val_np[i] = 1.;
    }

    temp_val_np[i] 	= 0.;
    postFit_val_np[i]   = 0.;
    aux_postFit_val_np[i]   = 0.;

    //printing the values 
    cout << par->getName() << "   Min  " << min_val_np[i] << "  Max  " << max_val_np[i] << endl;
  }

  cout << "---------> Generating...." << endl;


  TRandom3 rand;
/* /// This is very inefficient!!!!!!!!!!!
  //Generate random toys and compute the log(L) searching for minima
  for(int itr = 0; itr < numberOfToys; itr++){

	if(itr % 1000000 == 0) {
		cout << "Generated  " << itr << " toys of " << numberOfToys << endl;
		rand.SetSeed(seed+itr); // not really necessary... 
	}
	//get values
	for(int k=0; k < np; k++){
	    temp_val_np[k] = rand.Uniform(min_val_np[k], max_val_np[k]);  // very inefficient way, to be improved! FIXME
	}	

	//compute the -log(L) for that values of NP
	double temp_LL = computeLL(temp_val_np);  
	
	if(ML > temp_LL) {
		ML = temp_LL;
		for(int k=0; k < np; k++) postFit_val_np[k] = temp_val_np[k];
	}
  }
//----------------------------------------------------//
*/


  //Generate random toys and compute the log(L) searching for minima
  // Assuming independence of most NP

  rand.SetSeed(seed);
  for(int itr = 0; itr < numberOfToys; itr++){

	//Start with the Sigma parameter  (sigma is np-1)
        double step_0 = (max_val_np[np-1] - min_val_np[np-1])/((double)numberOfToys);
	temp_val_np[np-1] = min_val_np[np-1] + ((double)itr)*step_0; 
	//temp_val_np[np-1] = rand.Uniform(min_val_np[np-1], max_val_np[np-1]);

	//generate BKG normalization
  	for(int itr1 = 0; itr1 < numberOfToys; itr1++){
           double step_1 = (max_val_np[np-2] - min_val_np[np-2])/((double)numberOfToys);
	   temp_val_np[np-2] = min_val_np[np-2] + ((double)itr1)*step_1; 
	   //temp_val_np[np-2] = rand.Uniform(min_val_np[np-2], max_val_np[np-2]);
	   
           //generate Stat NP considering them independent from each other
	   //so for each of them the value that	maximize the log(L) is independet from the value of the others
	   for(int k=0; k < np-2; k++) temp_val_np[k]  = bestGuess_val_np[k];  //set to best guess

	   //find the best fit for each
	   for(int k=0; k < np-2; k++){

		double minimum = VERY_LARGE;
	   	for(int itr2 = 0; itr2 < numberOfToys; itr2++){

           	         double step_2 = (max_val_np[k] - min_val_np[k])/((double)numberOfToys);
	                 temp_val_np[k] = min_val_np[k] + ((double)itr2)*step_2; 
	    		//temp_val_np[k] = rand.Uniform(min_val_np[k], max_val_np[k]);  // very inefficient way, to be improved! FIXME

			double temp_LL = computeLL(temp_val_np);
			if(minimum > temp_LL) { minimum = temp_LL; aux_postFit_val_np[k] = temp_val_np[k] ;}
		}
	    }

 	    aux_postFit_val_np[np-2] = temp_val_np[np-2]; // set to temporary BKG norm
	    aux_postFit_val_np[np-1] = temp_val_np[np-1]; // set to temp Sigma
	    
            double temp_LL = computeLL(aux_postFit_val_np);

	     // This is the overall  minimizzation!
	     if(ML > temp_LL) {
                ML = temp_LL;
                for(int k=0; k < np; k++) postFit_val_np[k] = aux_postFit_val_np[k];
	     }
        	
	}

	if(itr % 10 == 0) {
		cout << "Generated  " << itr << " toys of " << numberOfToys << " current param values: "<< endl;
                for(int k=0; k < np; k++) cout <<"\t\tPar " << k << " =  " << postFit_val_np[k] << endl;
		cout << " Current ML = " << ML << endl;
	}

  }

  //Setting post fit values
  for(int i=0; i < np; i++){
    LKParameter*   par  = MinuitParameters[i];
    //if(par->getType() == PARAMETER_OF_INTEREST) {
    //	par->setCurrentValue(0.);
    //	  par->setInitialValue(0.);
    //}
    //else{
	par->setCurrentValue(postFit_val_np[i]);
	//par->setInitialValue(postFit_val_np[i]);
    //}
  }

  cout <<"\n--------------------\nPOST FIT numerically evaluated parameters " << endl;
  printCurrentParameters();

  return (-1.*ML); // return -1* Minimum (which is then the maximum)
}






//-------------- profile likelihood ----------------

ProfileLikelihood::~ProfileLikelihood(){}

ProfileLikelihood::ProfileLikelihood(TString nam): Likelihood(nam){
  setup();
}


void ProfileLikelihood::printFlagsAndParameters(){
  cout<<endl<<getName()<<endl;
}

void ProfileLikelihood::setup(){
  sigPar             = NULL;
}



double ProfileLikelihood::getWimpMass(){
	return 0.;
}




void ProfileLikelihood::estimateCrossSection() {
  if( getPrintLevel() < 2){
     cout<<"Fitting "<<getName()<<" once, for null hypothesis"<<endl;
  }
  resetParameters();
  LogD=maximize(false);
  LKParameter *parsig=getParameter(PAR_SIGMA);
  if(parsig==NULL) {
    cout<<"Technical bug!"<<getName()<<" does not have a 'sigma' param"<<endl; 
    return;
  }
  if( getPrintLevel() < 2){
     cout<<"Estimated sigma="<<sigmaHat<<" log(d)="<<LogD<<endl;
  }
}




double ProfileLikelihood::returnLimitHagar( double cl, bool forceMuhatComp){

    //This computes the limit in mu, not in Xsection, you need to scale	
    // valid for qtilde and q	
    double q_val    = pow(ROOT::Math::normal_quantile(1-cl,1), 2.); 

    double mu_limit = getSigmaAtQval(q_val, forceMuhatComp);

    if(mu_limit < 0. ) {
	    cout << "ProfiledLikelihood::returnLimitHagar - ERROR IN COMPUTING LIMITS" << endl;
	    exit(100);
    }
    return mu_limit; 

}
     

     
     
double ProfileLikelihood::getSigmaAtQval(double qval, bool forceMuhatComp){

  //return the POI for that q_val computed with q_tilde
  //good for limit and good for sensitivity
  
  TGraph *g=getGraphAtQval(qval, forceMuhatComp);
  printf ("diff= got %d points \n",g->GetN());
  TGraph *go=new TGraph(g->GetN(),g->GetY(), g->GetX());
  double bestGuess=go->Eval(qval);
  if (bestGuess<g->GetX()[0]) bestGuess=(-1);
  if (g->GetN()==20) return bestGuess=(-2);
  delete g;
  delete go;
  return bestGuess;
}



TGraph * ProfileLikelihood::getGraphAtQval(double qval, bool forceMuhatComp){

// compute the TGraph of qtilde for mu > mu_hat

  double llmin  = UNDEFINED; 
  
  // Do not maximize twice if somehow in some previous code you already did it.
  if(getLogD() !=UNDEFINED && getSigmaHat() != UNDEFINED && forceMuhatComp == false)
  	llmin = getLogD();
  else{ 
	cout << "recomputing.... mu_hat" << endl;
	llmin = maximize(false);

  }
  double sigmin= getSigmaHat();
  double sigma0,ll0;
  TGraph *g=new TGraph();

    cout << "sigma min = " << sigmin << endl; 
    cout << "sigma min = " << getParameterValue(PAR_SIGMA) << endl; 
  //--- Qtilde set mu_hat to zero if mu_hat <0 
  if (sigmin<0) {
    resetParameters();
    setParameterValue(PAR_SIGMA, 0 );
    llmin       =  maximize(true);
    sigmin=0;
	cout << "Resetting simgma min to zero, sigma_min=" << sigmin << "   due to qtilde" << endl;
  }

  double dstep=1;
  sigma0=sigmin;   // Qtilde --> sigma0 > sigmamin
  ll0=llmin;
  double pdiff= -2.*(ll0-llmin)-qval;
  int counter=0;
  g->SetPoint(g->GetN(),sigma0, -2.*(ll0-ll0));
  while (counter<20 && fabs(pdiff)>0.01) {
    counter++;
    sigma0=sigma0+dstep;
    resetParameters();
	  setParameterValue(PAR_SIGMA, sigma0 );

	  ll0       =  maximize(true);
	  printCurrentParameters();

	  g->SetPoint(g->GetN(),sigma0,-2.*(ll0-llmin));
	  double diff= -2.*(ll0-llmin)-qval;
	  if (pdiff*diff<0) dstep= -1.*dstep*0.5;
	  pdiff=diff;
	  printf (" \t  at: %f (min %f)  q_lim=%f  ll0=%f  qval=%f diff=%f next-step=%f \n",sigma0,sigmin, qval, ll0,-2*(ll0-llmin), pdiff, dstep);
	}


  return g;


}
TGraph* ProfileLikelihood::getGraphOfLogLikelihood( int n_points){
  // warning ! needs first to exclusion->prepareForLimit
  int n = n_points;

  TGraph* gr= new TGraph(n);

  double reference=0.;
  LKParameter* par=getParameter(PAR_SIGMA);
  double min = par->getMinimum();
  double max = par->getMaximum();

//  par->setMinimum(-5.);
  min=0;
  double step = (max - min) / ((double) n);

  double ll_Denominator = maximize(false);  // unconditional fit!!!
//  double ll_Denominator = maximizeNumerically(200,false);  // unconditional fit!!!
  double post_fit_sigma = par->getCurrentValue();

  cout << "min:  " << min << "   Max:  " << max << "   PostFit:   " << post_fit_sigma  << "  LL max " << ll_Denominator<< endl; 
//  ll_Denominator = maximize(false);  // unconditional fit!!!

  for(int i=0;i<n;i++){
  cout << "min:  " << min << "   Max:  " << max << "   PostFit:   " << post_fit_sigma  << "  LL max " << ll_Denominator<< endl; 
    resetParameters();

    double value = min + ((double)i)*step;

    setParameterValue(PAR_SIGMA, value);

    double ll_Numerator = maximize(true); // conditional fit!!!
    printCurrentParameters();
    double test_stat_q = -2. * (  ll_Numerator - ll_Denominator);
    gr->SetPoint(i,value , test_stat_q);

//   cout << " \t\t\t\t\t\t test val " << test_stat_q  << "  value  "  << value << endl;
  }


  cout << "min:  " << min << "   Max:  " << max << "   PostFit:   " << post_fit_sigma  << "  LL max " << ll_Denominator<< endl; 

  gr->SetTitle("Log Likelihood Scan on Cross Section");
  gr->GetXaxis()->SetTitle("#sigma #times 10^{-45} cm^{2}");
  gr->GetYaxis()->SetTitle("-2Log(L(#sigma)/L(#hat{#sigma})");

  return gr;

}

TGraph* ProfileLikelihood::getLikelihoodScanOfParameter( int n_points, LKParameter * par, double mu){
  // warning ! needs first to exclusion->prepareForLimit
  int n = n_points;

  TGraph* gr= new TGraph(n);

  double reference=0.;
  LKParameter* sig=getParameter(PAR_SIGMA);

  sig->setCurrentValue(0); // conditional fit set to zero for the denominator

  double min= par->getMinimum();
  double max = par->getMaximum();

  double step = (max - min) / ((double) n);

  double ll_Denominator = maximize(true); // conditional best fit

  int save_para_type = par->getType();
  par->setType(FIXED_PARAMETER);

  sig->setCurrentValue(mu); // conditional fit set to mu... Default is zero signal!

  for(int i=0;i<n;i++){
    resetParameters();

    double value = min + ((double)i)*step;

    sig->setCurrentValue(mu); // conditional fit set to mu... Default is zero signal!
    
    par->setCurrentValue( value);

    double ll_Numerator = maximize(true); // conditional fit!!!
    printCurrentParameters();
    //double test_stat_q = -1.* ll_Numerator;
     double test_stat_q = -2. * (  ll_Numerator - ll_Denominator);
    gr->SetPoint(i,value , test_stat_q);

//   cout << " \t\t\t\t\t\t test val " << test_stat_q  << "  value  "  << value << endl;
  }


  gr->SetTitle("Log Likelihood Scan on "+TString(par->getName()));
  gr->GetXaxis()->SetTitle(TString(par->getName()));
  gr->GetYaxis()->SetTitle("-2Log(L(#sigma)/L(#hat{#sigma})");

  par->setType(save_para_type);
  return gr;

}

TGraph* ProfileLikelihood::getGraphOfParameter( int n_points, int param_index){
  // warning ! needs first to exclusion->prepareForLimit
  int n = n_points;

  TGraph* gr= new TGraph(n);
  gr->GetXaxis()->SetTitle("#mu");
  gr->GetYaxis()->SetTitle("param tval");

  double reference=0.;
  LKParameter* par=getParameter(PAR_SIGMA);
  double min = par->getMinimum();
  double max = par->getMaximum();

  double step = (max - min) / ((double) n);

  for(int i=0;i<n;i++){
    resetParameters();

    double value = min + ((double)i)*step;
    setParameterValue(PAR_SIGMA, value);
    double ll=maximize(true);
     
   LKParameter *par = getParameter(param_index);
    double post_fit_val = par->getCurrentValue();
    gr->SetPoint(i,value , post_fit_val);
  }

  return gr;
}


//--------- Combining several profile likelihoods --------------------

CombinedProfileLikelihood::CombinedProfileLikelihood(TString n)
                   : ProfileLikelihood(n) {
  combinedMode=false;
  nCommon=0;
  setExperiment(ALL);
}

CombinedProfileLikelihood::~CombinedProfileLikelihood(){
  clear();
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood *pl=it->second;
    pl->clearTheParameters();
  }
}

ProfileLikelihood* CombinedProfileLikelihood::getProfile(int ex){
  if(exps.find(ex)==exps.end()) {
    cout<<"Couldn't find experiment "<<ex<<endl;
    return NULL;
  } 
  return exps[ex];
}

void CombinedProfileLikelihood::printFlagsAndParameters(){
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood *pl=it->second;
    cout<<endl<<"Profile Likelihood for experiment "<<pl->getExperiment()<<endl;
    pl->printFlagsAndParameters();
  }
}


void   CombinedProfileLikelihood::combine(ProfileLikelihood* pl) {
  int exp=pl->getExperiment();
  if(exps[exp]!=NULL) {
    cout<<"Can't combine experiment "<<exp<<": already exists!"<<endl;
    return;
  }
  exps[exp]=pl;
  combinedMode=true;
  pl->setCombinedMode();
}

bool CombinedProfileLikelihood::checkPValue(){
 Info("CombinedProfileLikelihood", "Building up combined PL " + getName() );

  // first common parameters
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood *pl=it->second;
    if(!pl->initialize()) return false;
    map<int,LKParameter*> *params= pl->getParameters();
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
      LKParameter *param=ip->second;
      if(param->isCommon()){
        int id=param->getId();
        if(parameters.find(id)==parameters.end()){
          if(getPrintLevel() < 2) {
             cout<<" adding common parameter "<<id<<" ("<<param->getName()
                 <<")"<<endl;
          }
          addParameter(param);
          param->setExperiment(ALL);
        }     
        else if(param->compares(parameters[id],true)) {
          if(getPrintLevel() < 2) {
             cout<<" skipping existing common parameter "<<id<<" ("
                 <<param->getName()<<")"<<endl;
          }
          pl->replaceParameter(parameters[id]);
        }
        else {
          cout<<"Inconsistency in common parameter "<<param->getName()<<endl;
          return false;
        }
      }
    }
  }
  TRAVERSE_PARAMETERS(it) {
    LKParameter *p=it->second; 
    if(p->isCommon()) {
      nCommon++; 
      p->setCombinedMode();
    }
  }
  if(!checkConsistency()) return false;
  if(nCommon==0){  
    cout<<"No common parameter in combined PL "<<getName()<<endl;
    return false;
  }
  sigPar=getParameter(PAR_SIGMA);
  if(sigPar==NULL) {
    cout<<"There should be a 'SIGMA' parameter in "<<getName()<<endl;
    return false;
  }

  // then non common parameters
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    int exp=it->first;
    if(!pl->initialize()) return false;
    map<int,LKParameter*> *params= pl->getParameters();
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
      LKParameter *param=ip->second;
      if(!param->isCommon()) {
        int id=ip->first;
        addParameter(param,AUTO);
        param->setExperiment(exp);
        if(getPrintLevel() < 2) {
           cout<<" adding specific parameter "<<id<<" ("<<param->getName()
               <<") for experiment "<<exp<<", new id="<<param->getId()<<endl;
        }
      } 
    }
  }
  if(getPrintLevel() < 2) printInitialParameters();
  return true;
}

double CombinedProfileLikelihood::computeTheLogLikelihood(){

  double ll=0;

  if(exps.size() > 0) {
     TRAVERSE_EXPERIMENTS(it) {
        ProfileLikelihood* pl=it->second;
        ll +=pl->computeTheLogLikelihood();
      }
  }
  else { cout << "CombinedProfileLikelihood::computeTheLogLikelihood()::  FATAL ERROR - no likelihood defiend ----> QUIT!"<< endl;
	  exit(100);
  }

  return ll;
}
/*
void CombinedProfileLikelihood::setWimpMass(double mass){
  if(printLevel>0) {
    cout<<"Setting mass for "<<getName()<<": "<<mass<<endl;
  }
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    pl->setWimpMass(mass);
  }
}

*/
void CombinedProfileLikelihood::setData(int dataType){

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    pl->setData(dataType);
  }
}

void CombinedProfileLikelihood::generateAsimov(double mu_prime){

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    pl->generateAsimov(mu_prime);
  }
}

void CombinedProfileLikelihood::generateToyDataset(double seed, double mu_prime){

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    pl->generateToyDataset(seed, mu_prime);
  }
}


double CombinedProfileLikelihood::getWimpMass(){

  double mass = 0.;
  int counter =0;
  double tempMass =0.;

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;

    if(counter == 0) 
          mass = pl->getWimpMass();
    else{
	  tempMass = pl->getWimpMass();
	  if(tempMass != mass) {
		cout << "CombinedProfileLikelihood::getWimpMass - FATAL ERROR: masses for different experiment are different" << endl;
		exit(100);
	    }
	}
    counter++;
    
  }

  return mass;
}


double CombinedProfileLikelihood::getSignalDefaultNorm(){

  double norm = 0.;
  int counter =0;
  double tempNorm =0.;

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;

    if(counter == 0) 
          norm = pl->getSignalDefaultNorm();
    else{
	  tempNorm = pl->getSignalDefaultNorm();
	  if(tempNorm != norm) {
		cout << "CombinedProfileLikelihood::getSignalDefaultNorm - FATAL ERROR: Norm for different experiment are different" << endl;
		exit(100);
	    }
	}
    counter++;
    
  }

  return norm;
}

double CombinedProfileLikelihood::getSignalMultiplier(){

  double norm = 0.;
  int counter =0;
  double tempNorm =0.;

  bool setToSmallest = false;

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;

    if(counter == 0) 
          norm = pl->getSignalMultiplier();
    else{
	  tempNorm = pl->getSignalMultiplier();
	  if(tempNorm < norm) {
		cout << "CombinedProfileLikelihood::getSignalMultiplier - WARNING: Norm for different experiment are different, setting to smallest" << endl;
		norm = tempNorm;
		setToSmallest = true;
	    }
	}
    counter++;
    
  }

  if(setToSmallest){
     TRAVERSE_EXPERIMENTS(it) {
         ProfileLikelihood* pl=it->second;
    	 pl->setSignalMultiplier(norm);

     }
  }

  return norm;
}

void CombinedProfileLikelihood::setSignalMultiplier(double val){

  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    pl->setSignalMultiplier(val);

  }


}
/*
bool CombinedProfileLikelihood::update(){
  sigToEvents=0;
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    initializedOK= initializedOK && pl->update();
    sigToEvents += pl->nSignalPerCm2();
  }
  if(printLevel>0) {
    cout<<"Overall sigma to events for "<<getName()<<": "<<sigToEvents<<endl;
  }
  return initializedOK;
}
*/


