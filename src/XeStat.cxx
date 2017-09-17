#include <fstream>
#include "XeStat.h"
#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"


int    XeStat::analysisMode = NO_ANALYSIS;
int    XeStat::printLevel=0;

       XeStat::~XeStat()                          {}
       XeStat::XeStat()           : XeObject()    {}
       XeStat::XeStat(string nam) : XeObject()    {setName(nam);}
void   XeStat::setAnalysisMode(int m){analysisMode=m;}
void   XeStat::setPrintLevel(int l)  {printLevel=l;}
bool   XeStat::isCutsBased()         {return analysisMode==CUTS_ANALYSIS;}
bool   XeStat::isPL()                {return analysisMode==PL_ANALYSIS;}
int    XeStat::getAnalysisMode()     {return analysisMode;}
int    XeStat::getPrintLevel()       {return printLevel;}
string XeStat::getAnalysisModeName() {return getAnalysisModeName(analysisMode);}

string XeStat::getSystematicModeName(int i){
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

XeRange* XeStat::newSigmaRange(XeRange *range, int unit){
  if(range!=NULL) return range;
  if(unit==EVENT_UNIT) return new LinearRange("Number of events",50,0.,10.);
  cout<<"No default Sigma range in SIGMA_UNIT"<<endl;
  return NULL;
}

string  XeStat::getSigmaUnitName(int unit) {
  switch(unit) {
    case SIGMA_UNIT  : return "cross section";
    case EVENT_UNIT  : return "number of events";
  }
  return UNDEFINED_STRING;
}

string  XeStat::getSigmaModeName(int mode) {
  switch(mode) {
    case ESTIMATED   : return "Estimated";
    case UPPER_LIMIT : return "Upper Limit";
  }
  return UNDEFINED_STRING;
}

string XeStat::getSigmaLabel(int unit){
  switch(unit) {
    case SIGMA_UNIT  : return LABEL_SIGMA;
    case EVENT_UNIT  : return LABEL_EVT;
  }
  return UNDEFINED_STRING;
}

string XeStat::getUpperSigmaLabel(int unit){
  switch(unit) {
    case SIGMA_UNIT  : return LABEL_SIGMA_MAX;
    case EVENT_UNIT  : return LABEL_EVT_MAX;
  }
  return UNDEFINED_STRING;
}

int XeStat::getSigmaLinLog(int unit) {return(unit==SIGMA_UNIT? LOG:LINEAR);}

string  XeStat::getAnalysisModeName(int m) {
  switch(m) {
    case NO_ANALYSIS   : return "Undefined";
    case PL_ANALYSIS   : return "Profile likelihood";
    case CUTS_ANALYSIS : return "Cuts based";
  }
  return UNDEFINED_STRING;
}

bool XeStat::checkAnalysisMode(string name, int requestedAnalysisMode){
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


//---------------------------- Rejection ----------------------------------
double Rejection::rejections[N_REJECTION_MODES]={.995,.9975,.999,1.};
 
Rejection::Rejection() :XeStat("Rejection") {}
Rejection::~Rejection(){}

int Rejection::getMode(double e){
  for(int i=0;i<N_REJECTION_MODES;i++){
    if(fabs(e-rejections[i])<1.E-6||fabs(e-100.*rejections[i])<1.E-4) return i;
  }
  cout<<"Rejection mode for "<<e<<" not yet implemented"<<endl;
  return UNDEFINED_INT;
}

double Rejection::getRejection(int mode) {return rejections[mode];}

string Rejection::getModeName(int mode) {
  if(mode<=0 || mode<N_REJECTION_MODES) return getModeName(rejections[mode]);
  return "Invalid mode";
}
string Rejection::getModeName(double r) {
  char title[30];
  sprintf(title,"Rejection at %.3g%%",100.*r); 
  return string(title);
}

//------  wrapper for distribitions ---------------------------

TRandom3 XeDist::random;

XeDist::XeDist()        :  XeStat() {}
XeDist::XeDist(string n):  XeStat(){setName(getTheName(n));}
XeDist::~XeDist(){}

string XeDist::getTheName(string nam)   {return "Distribution: "+nam;}
double XeDist::cdf(double x1,double x2) {return cdf(x2)-cdf(x1);}
double XeDist::rndm()                   {return random.Rndm();}

double XeDist::findSigma(TH1* h, double sigma){
  return findQuantile(h,GaussianDist::below(sigma));
}

double XeDist::logPdf(double x) {
  double p=pdf(x);
  return p>0.? log(p) : UNDEFINED ;
} 

double XeDist::findQuantile(TH1* h, double q){
  int bin=1;
  int nBins=h->GetNbinsX();
  if(q<=0.) {}
  else if(q>=1.) bin=nBins;
  else {
    double under=h->GetBinContent(0);
    double over=h->GetBinContent(nBins+1);
    double th=(under+h->Integral()+over)*q;
    double sum=under;
    for(bin=1;bin<nBins;bin++){
      sum+= h->GetBinContent(bin);
      if(sum>th) break;
    }
  }
  return h->GetBinCenter(bin);
}

void XeDist::normalize(TH1F* h,double n) {h->Scale(n/h->GetSumOfWeights());}

void XeDist::normalize(TGraph* g, double norm){
  double *x=g->GetX();
  double *y=g->GetY();
  int n=g->GetN();
  double sumy=0.;
  for(int i=0;i<n;i++) sumy+=y[i];
  if(sumy==0.) return;
  double w=norm/sumy;
  for(int i=0;i<n;i++) g->SetPoint(i,x[i],y[i]*w);
}


vector<double> XeDist::generateVector(int n){
  vector<double> v(n);
  for(int i=0;i<n;i++) v[i]=generate();
  return v;
}



//----------------------handling set of values with multiplicity -------------

XeValues::~XeValues(){}
XeValues::XeValues() : map<double,int>(), XeObject("") {reset();}
XeValues::XeValues(string nam) : map<double,int>(), XeObject(nam) {reset();}

XeValues::XeValues(string h, double *d,int n) : map<double,int>(), XeObject(h){
  reset();
  add(d,n);
}

XeValues::XeValues(string nam, DataSet* dataSet,int c)   
   : map<double,int>(), XeObject(nam) {
  reset();
  add(dataSet->getVector(c),dataSet->getNEvents());
}

XeValues::XeValues(string h,vector<double>& vd) : map<double,int>(),XeObject(h){
  reset();
  add(&(vd[0]),(int)vd.size());
}

int  XeValues::getNValues()         {return (int)size();}
int  XeValues::getNEntries()        {return nEntries;}
void XeValues::reset()              {clear(); nEntries=0;}
void XeValues::add(double d)        {(*this)[d]++; nEntries++;}
void XeValues::add(double *d,int n) {
  for(int i=0;i<n;i++) (*this)[d[i]]++;
  nEntries += n;
}

double XeValues::largestGap(){
  map<double,int>::iterator it=begin();
  double vmin=it->first;
  double lg=0.;
  it++;
  for(; it!=end();it++){
    double t=it->first;
    lg=max(lg,t-vmin);
    vmin=t;
  }
  return lg;
}

double XeValues::quantile(double f){
  map<double,int>::iterator it=begin();
  if(f<=0.) return it->first;
  else if(f<1.){
    int threshold=f*nEntries;
    int current=0;
    for(;it!=end();it++) {
      current+= it->second;
      if(current>=threshold) return it->first;
    }
 }
 it=end() ; 
 it--;
 return it->first;
}

double XeValues::findSigma(double sigma){
  return quantile(GaussianDist::below(sigma));
}

bool XeValues::printIt(int level){
  if(level<1) return true;
  cout<<Name<<" has "<<nEntries<<" entries with "<<size()<<" different values"
      <<endl;
  if(level<2) return true;
  int n=0;
  for(map<double,int>::iterator it=begin();it!=end();it++){
    double d=it->first;
    int    w=it->second;
    cout<<formatG(d,10,4);
    if(w>1) cout<<" (*"<<w<<") ";
    else    cout<<"      ";
    n++;
    if(n%5==0) cout<<endl;
  }
  if(n%5!=0) cout<<endl;
  return n>0;
}


// --------------- Tabulated Spectrum and distribution -------------------

XeSpectrum::~XeSpectrum() { if(rangeIsCreated) deleteWithPointer(range);}

void XeSpectrum::initialize(){
  range=NULL;
  rangeIsCreated=false;
  nBins=0;
}
XeSpectrum::XeSpectrum() : XeStat(){initialize();}

XeSpectrum::XeSpectrum(string n) : XeStat(getTheName(n)){ initialize(); }

XeSpectrum::XeSpectrum(string n, int nB, double min,double max,double *spec
  , int mode, int nValues) : XeStat() {
  setName(getTheName(n));
  setTheRange(nB,min,max);
  if(mode==SPECTRUM)    importSpectrum(spec);
  else if(mode==VALUES) fillValues(spec,nValues);
  else cout<<"Unknonwn mode in XeSpectrum::XeSpectrum : "<<mode<<endl;
}

XeSpectrum::XeSpectrum(string n, int nB, double min,double max
            ,vector<double>& spec, int mode) :XeStat(getTheName(n)) {
  int size=spec.size();
  if(mode==SPECTRUM && nB!=size){
    cout<<"Inconsistency between given number of bins ("<<nB
        <<") and size of spectrum ("<<size<<")"<<endl;
    nBins=0;
    return;     
  }
  setTheRange(nB,min,max);
  if(mode==SPECTRUM)    importSpectrum(spec);
  else if(mode==VALUES) fillValues(&(spec[0]),spec.size());
  else cout<<"Unknonwn mode in XeSpectrum::XeSprectrum : "<<mode<<endl;
}


XeSpectrum::XeSpectrum(string nam,int nB,double min,double max,XeValues *values)
           : XeStat(getTheName(nam)) {
  int n=values->size();
  if(n<2) {
    cout<<"Error! Can't create a SpectrumFromValues with so few points: "
        <<n<<endl;
    return;
  }
  setTheRange(nB,min,max);
 
  for(map<double,int>::iterator it=values->begin();it!=values->end();it++){
    double val=it->first;
    double weight=it->second;   
    int bin=range->getIndex(val);
    if(bin>=0 && bin<nBins) spectrum[bin]+= weight;
  }
}

XeSpectrum XeSpectrum::operator+(XeSpectrum& other) {
  XeSpectrum spect(*this);
  if(range == other.range){
    double *s=spect.getSpectrum();
    double *o=other.getSpectrum();
    for(int b=0;b<nBins;b++) s[b]+=o[b];
  }
  else cout<<"Can't add spectra "<<getName()<<" with "<<other.getName();
  return spect;
}

XeSpectrum XeSpectrum::operator*(double a ) {
  XeSpectrum spect(*this);
  double *s=spect.getSpectrum();
  for(int b=0;b<nBins;b++) s[b]+=a*s[b];
  return spect;
}

void XeSpectrum::fillValues(double* values, int nValues){
  for(int i=0;i<nValues;i++) {
    int bin=range->getIndex(values[i]);
    if(bin>=0 && bin<nBins) spectrum[bin]+= 1.;
  }
}

void XeSpectrum::setTheRange(int nB, double min,double max){
  initialize();
  range=new LinearRange("Spectrum",nB+1,min,max);
  XeRange *rdef=range->findDefaultRange();
  if(rdef==NULL) rangeIsCreated=true;
  else {
    delete range;
    range=(LinearRange*) rdef;
  }
  setTheRange(range);
}

XeSpectrum::XeSpectrum(string n,LinearRange *lr, double *spec, int nValues
           , int mode) : XeStat(getTheName(n)){
  initialize();
  setTheRange(lr);
  if(mode==SPECTRUM)    importSpectrum(spec);
  else if(mode==VALUES) fillValues(spec,nValues);
  else cout<<"Unknonwn mode in XeSpectrum::XeSprectrum : "<<mode<<endl;
}

XeSpectrum::XeSpectrum(string n,LinearRange *lr ,vector<double>& spec, int mode)
  :XeStat(getTheName(n)) {
  initialize();
  setTheRange(lr);
  if(mode==SPECTRUM)    importSpectrum(spec);
  else if(mode==VALUES) fillValues(&(spec[0]),spec.size());
  else cout<<"Unknonwn mode in XeSpectrum::XeSprectrum : "<<mode<<endl;
}

string XeSpectrum::getTheName(string nam)    { return SPECTRUM_HEADER+nam;}

string XeSpectrum::getReducedName(){
  return removeHeader(getName(),SPECTRUM_HEADER);
} 

void XeSpectrum::setTheRange(LinearRange *lr){
  range=lr;
  nBins=range->getNPoints()-1;
  spectrum.resize(nBins);
  for(int b=0;b<nBins;b++) spectrum[b]=0.;
}

void XeSpectrum::importSpectrum(double *spec){
  if(spec!=NULL)  for(int i=0;i<nBins;i++) spectrum[i]=spec[i];
}

void XeSpectrum::importSpectrum(vector<double>& spec){
  reset();
  addSpectrum(spec);
}

void XeSpectrum::addSpectrum(vector<double>& spec){
  int vl=spec.size();
  if(vl!=nBins) {
    cout<<"Bins length incompatibility in XeSpectrum::addSpectrum : "<<vl
        <<" vs "<<nBins<<endl;
    return;
  }
  for(int i=0;i<nBins;i++) spectrum[i]+=spec[i];
}


void XeSpectrum::reset() {for(int i=0;i<nBins;i++) spectrum[i]=0.;}

void XeSpectrum::fillExponent(double x0, double norm,double xmin, double xmax){
  double xmi=min(xmin,range->getMin());
  double xma=max(xmax,range->getMax());
  pair<int,double> pi=range->getIndexAndFraction(xmi);
  int bi=pi.first;
  pair<int,double> pa=range->getIndexAndFraction(xma);
  int ba=pa.first;
  if(bi==ba) spectrum[bi]+= norm*(exp(xma/x0)-exp(xmi/x0))/x0;
  else {
    double fi=pi.second;
    double fa=pa.second;
    spectrum[bi]+= (1.-fi)*norm/x0*(exp(range->getNextValue(bi))-exp(xmi));
    spectrum[ba]+= fa*norm/x0*(exp(xma)-exp(range->getValue(ba)));
    if(bi<ba-1) {
      for(int b=bi+1;b<ba-1;b++) {
        spectrum[b]+= norm/x0
                    * (exp(range->getNextValue(b)/x0)-exp(range->getValue(b)/x0));
      }
    }
  }  
}

XeGraph* XeSpectrum::newGraph(string Xlab, string Ylab , XeStyle* sty
                              , string leg){
  XeGraph *g=new XeGraph(Name,nBins,Xlab,Ylab,sty,leg);
  for(int b=0;b<nBins;b++) g->setPoint(b,range->centerOfBin(b),spectrum[b]);
  return g;
}

void XeSpectrum::drawGraph(string options,string Xlab,string Ylab,XeStyle* sty){
   XeGraph *g=newGraph(Xlab,Ylab,sty);
   g->draw(options);
}

void XeSpectrum::drawGraphWithFrame(string options,string Xlab
                                   ,string Ylab,XeStyle* sty){
   XeGraph *g=newGraph(Xlab,Ylab,sty);
   g->drawWithFrame(options);
}

void XeSpectrum::add(XeSpectrum* other){
  if(! (*range == *other->getRange()) ) {
    cout<<"Can't add spectra of different specifications"<<endl;
    return;
  }
  double *original=getSpectrum();
  double *others=other->getSpectrum();
  for(int b=0;b<nBins;b++) original[b] += others[b];
}

TH1F*  XeSpectrum::newHistogram(string n,int plot){
  if(n=="") n=Name;
  TH1F *h=new TH1F(n.c_str(),n.c_str(),nBins,range->getMin(),range->getMax()); 
  for(int b=0;b<nBins;b++) h->Fill(range->centerOfBin(b),spectrum[b]);
  drawHist(h,plot);
  return h;
}

void XeSpectrum::drawHistogram(string options,string Xlab, string Ylab){
  TH1F *h=newHistogram();
  drawHist(h,Xlab,Ylab,options);
}

TH2F*  XeSpectrum::newHistogram(int nSpectra,int plot){
  TH2F* h=new TH2F(getNameChar(),getNameChar(),nBins
                  ,range->getMin(),range->getMax(),nSpectra,0.,nSpectra); 
  drawHist(h,plot);
  return h;
}

void XeSpectrum::fillHistogram(TH2F* hist, int spectr){
  for(int b=0;b<nBins;b++) {
    hist->Fill(range->centerOfBin(b),.5+spectr,spectrum[b]);
  }
}

double XeSpectrum::sum()        {return vectorSum(spectrum);}

void   XeSpectrum::setBin(int b, double content){
  if(b>=0 && b<nBins) spectrum[b]=content;
}

void XeSpectrum::fillDelta(double x,double w){
  pair<int,double> p=range->getIndexAndFraction(x);
  int b=p.first;
  if(b>=0 && b<nBins) spectrum[b]+=w;
}

pair<int,double> XeSpectrum::getBinAndFraction(double x){
  return range->getIndexAndFraction(x);
}

int             XeSpectrum::getBin(double x)   {return range->getIndex(x);}
int             XeSpectrum::getNBins()         {return range->getNBins();}
double          XeSpectrum::getMin()           {return range->getMin();}
double          XeSpectrum::getMax()           {return range->getMax();}
double*         XeSpectrum::getSpectrum()      {return &(spectrum[0]);}
LinearRange*    XeSpectrum::getRange()         {return range;}
vector<double>* XeSpectrum::getTheSpectrum()   {return &spectrum;}

double XeSpectrum::integral(double x1,double x2){

  if(x1==x2) return 0.;

  double sign=1.;
  double xmi=x1;
  double xma=x2;
  if(xmi>xma) {
    sign=-1.;
    xmi=x2;
    xma=x1;
  }

  pair<int,double> bfi=range->getTolerantIndexAndFraction(xmi);
  pair<int,double> bfa=range->getTolerantIndexAndFraction(xma);
  int    i1=bfi.first;
  int    i2=bfa.first;
  double f1=bfi.second;
  double f2=bfa.second;

  if(i1==i2) return sign*(f2-f1)*spectrum[i1];
  double e=(1.-f1)*spectrum[i1]+f2*spectrum[i2];
  for(int i=i1+1;i<i2;i++) e+=spectrum[i];
  return e*sign;
}


bool XeSpectrum::printIt(int level){
  bool ok=nBins>0;
  if(level<1) return ok;
  cout<<"    "<<nBins<<" bins of width "<<range->getStep()
      <<" from "<<range->getMin()
      <<" to "<<range->getMax()<<endl;
  double s=sum();
  cout<<"Sum: "<<s<<endl;
  if(s>0.){
    switch(level) {
      case 2:  printVector(spectrum); break;
      case 3:  printDetailed(true);   break;
      case 4:  printDetailed(false);  break;
    }
  }
  return ok;
}

void XeSpectrum::printDetailed(bool ignoreTail){
  int bStop=nBins;
  if(ignoreTail) {
    for(int b=nBins-1;b>=0;b--) if(spectrum[b]>0.) { bStop=b+1; break;}
  }
cout<<"bStop="<<bStop<<endl;
  double cumul=0.;
  for(int b=0;b<bStop;b++){
    cumul += spectrum[b];
    cout<<setw(4)<<b<<"  ["<<formatF(range->getValue(b),7,2)
                    <<","  <<formatF(range->getNextValue(b),7,2)
                    <<" ]    "<<formatG(spectrum[b],10,4)
                    <<"  " <<formatG(cumul,10,4)
                    <<endl;
  }
}

//------------------ Distributions --------------------------

TabulatedDist::~TabulatedDist(){}
TabulatedDist::TabulatedDist():XeSpectrum(), XeDist() {}
TabulatedDist::TabulatedDist(string nam): XeSpectrum(nam), XeDist(nam) {
  setName(XeDist::getTheName(nam));
}

TabulatedDist::TabulatedDist(string nam,LinearRange *lr, double *spect
             ,int nValues, int mode) : XeSpectrum(nam,lr), XeDist(nam) {
  setName(XeDist::getTheName(nam));
  if(mode==SPECTRUM) importSpectrum(spect);
  else if(mode==VALUES) fillValues(spect, nValues);
}
 
TabulatedDist::TabulatedDist(XeSpectrum *sp) : 
  XeSpectrum(sp->getName(), sp->getNBins(), sp->getMin(), sp->getMax() ) {
  setName(XeDist::getTheName(sp->getReducedName()));
  importSpectrum(sp->getSpectrum());
}

TabulatedDist::TabulatedDist(string nam, int nB, double min, double max
           , double *spect, int nValues, int mode) : XeSpectrum(nam,nB,min,max){
  setName(XeDist::getTheName(nam));
  if(mode==SPECTRUM) importSpectrum(spect);
  else if(mode==VALUES) fillValues(spect, nValues);
}


TabulatedDist::TabulatedDist(string nam, int nB, double min, double max
                         , XeValues *values) : XeSpectrum(nam,nB,min,max){
  setName(XeDist::getTheName(nam));
  int n=values->size();
  if(n<1) {
    cout<<"Error! Can't create a a TabulatedDist without any point"<<endl;
    return;
  }
 
  for(map<double,int>::iterator it=values->begin();it!=values->end();it++){
    double val=it->first;
    double weight=it->second;   
    int bin=range->getIndex(val);
    if(bin>=0 && bin<nBins) spectrum[bin]+= weight;
 }
  normalizeAndCumulate();
}


double* TabulatedDist::getCumulated() {return &(cdfVector[0]);} 

map<double,int>* TabulatedDist::getCdfMap()    {return &cdfMap;}
vector<double>*  TabulatedDist::getCdfVector() {return &cdfVector;};

void TabulatedDist::fillValues(double* spec,int nValues){
  XeSpectrum::fillValues(spec,nValues);
  normalizeAndCumulate();
}
void TabulatedDist::importSpectrum(double* spec){
  if(spec==NULL) return;
  XeSpectrum::importSpectrum(spec);
  normalizeAndCumulate();
}

void TabulatedDist::importSpectrum(XeSpectrum* spec){
  if(spec==NULL) return;
  if(getRange()!=spec->getRange()) {
    cout<<"Spectrum incompabitility bewteen TabulatedDist "<<getName()
        <<" and spectrum "<<spec->getName()<<endl;
  }
  importSpectrum(spec->getSpectrum());
}

void TabulatedDist::importSpectrum(vector<double>& spec){
  XeSpectrum::importSpectrum(spec);
  normalizeAndCumulate();
}

void TabulatedDist::normalizeAndCumulate(){
  cdfVector.resize(nBins);
  if(!normalizeVector(spectrum)){
    if(printLevel>1) cout<<"Can't normalize empty tabulated dist."<<getName()<<endl;
    for(int b=0;b<nBins;b++) cdfVector[b]=0.;
    return;
  }
  double c=0.;
  if(printLevel>4) cout<<"Building "<<Name<<endl;
  for(int b=0;b<nBins;b++){
    double w=spectrum[b];
    if(w<0.) {
      cout<<"Negative pdf in bin "<<b<<" of distribution "<<getName()<<endl;
    }
    else{
      if(w>0.) {
        c+=spectrum[b];
        cdfMap[c]=b;
      }
      cdfVector[b]=c;
      if(printLevel>4){
         cout<<setw(3)<<b<<", x= "<<setw(6)<<range->getValue(b)
             <<", pdf:"<<setw(10)<<spectrum[b]
             <<", cumulated :"<<setw(10)<<c<<endl;
      }
    }
  }
  if(c<0.999999 || c>1.00001) {
    cout<<"Invalid cumulated prob in "<<Name<<" = "<<formatF(c,11,9)<<endl;
    printIt(3);
  }
  if(printLevel>4) {
    cout<<"Cumulated cdf of "<<Name<<endl;
    map<double, int>::iterator it;
    for(it=cdfMap.begin();it!=cdfMap.end();it++){
      int bin=it->second;
      c=it->first;
      cout<<setw(3)<<bin<<setw(6)<<range->getValue(bin)<<setw(10)<<c<<endl;
    }
  }
}

void TabulatedDist::importDistribution(XeDist* dist) {
  for(int b=0;b<nBins;b++){
    spectrum[b]=dist->cdf(range->getValue(b),range->getNextValue(b));
  }
  normalizeAndCumulate();
}

void TabulatedDist::addDelta(double x, double w){
  int b=range->getIndex(x);
  if(b<0||b>nBins) return;
  spectrum[b]+=w;
  normalizeAndCumulate();
}

bool TabulatedDist::printIt(int level){
  bool ok=true;
  if(level<1) return ok;
  cout<<"Tabulated distribution:"<<getName()<<endl;
  ok=ok && XeSpectrum::printIt(level);
  if(level>2) ok=ok && printCumulated(level-1);
  return ok;
}

bool  TabulatedDist::printCumulated(int level)  {
  if(level<1) return true;
  cout<<"Cumulated values of "<<getName()<<", "<<nBins<<" bins of width "
       <<range->getStep()<<" from "<<range->getMin() <<" to "<<range->getMax()
       <<", last  cumulation="<<cdfVector[nBins-1]<<endl;
  if(level<2) return true;
  return printVector(cdfVector);
}

double TabulatedDist::cdf(double x){
  pair<int,double> p=range->getIndexAndFraction(x);
  int bin=p.first;
  if(bin<0) return 0.;
  if(bin==nBins) return 1.;
  double f=p.second;
  return cdfVector[bin]*(1.-f)+cdfVector[bin+1]*f;
}

double TabulatedDist::pdf(double x){
  pair<int,double> p=range->getIndexAndFraction(x);
  int bin=p.first;
  if(bin<0) return 0.;
  if(bin==nBins) return 1.;
  double f=p.second;
  return spectrum[bin]*(1.-f)+spectrum[bin+1]*f;
}

double TabulatedDist::generate(){
  double r=rndm();
  return quantileX(r);
}

double TabulatedDist::quantileX(double r){
  if(r<0.||r>1.) return UNDEFINED;
  return (quantileBin(r)+rndm())*range->getStep();
}

int TabulatedDist::quantileBin(double r){
  if(r<0.||r>1.) return -1;
  map<double,int>::iterator it=cdfMap.lower_bound(r);
  int bin=it->second;
  if(it==cdfMap.end()) bin=nBins-1;
  else if(bin==nBins) bin--;
  return bin;
}

XeGraph* TabulatedDist::newCumulatedGraph(string Xlab, string Ylab 
                                         , XeStyle* sty , string leg){
  XeGraph *g=new XeGraph(Name,nBins,Xlab,Ylab,sty,leg);
  for(int b=0;b<nBins;b++) g->setPoint(b,range->centerOfBin(b),cdfVector[b]);
  return g;
}

void TabulatedDist::drawCumulatedGraph(string options,string Xlab,string Ylab
                                     ,XeStyle* sty){
   XeGraph *g=newGraph(Xlab,Ylab,sty);
   g->draw(options);
}

void TabulatedDist::drawCumulatedGraphWithFrame(string options,string Xlab
                                   ,string Ylab,XeStyle* sty){
   XeGraph *g=newGraph(Xlab,Ylab,sty);
   g->drawWithFrame(options);
}

TH1F*  TabulatedDist::newCumulatedHistogram(string n,int plot){
  if(n=="") n=Name;
  TH1F *h=new TH1F(n.c_str(),n.c_str(),nBins,range->getMin(),range->getMax()); 
  for(int b=0;b<nBins;b++) h->Fill(range->centerOfBin(b),cdfVector[b]);
  drawHist(h,plot);
  return h;
}

void TabulatedDist::drawCumulatedHist(string options,string Xlab, string Ylab){
  TH1F *h=newHistogram();
  drawHist(h,Xlab,Ylab,options);
}


//---------------------  Simple uniform dist -------------------

UniformDist::UniformDist()                 :XeDist("Uniform"){setLimits(0.,1.);}
UniformDist::UniformDist(double a,double b):XeDist("Uniform"){setLimits(a,b);}
UniformDist::~UniformDist(){}
 
double UniformDist::cdf(double x)    {
  if(x<x0) return 0.;
  else if(x>x1) return 1.;
  return (x-x0)*onew;
}
void   UniformDist::setLimits(double a, double b){x0=a;x1=b;onew=1./(x1-x0);}
double UniformDist::pdf(double x)                {return (x<x0||x>x1)? 0.:onew;}
double UniformDist::generate()                   {return generate(x0,x1);}
double UniformDist::generate(double a, double b) {return a+(b-a)*rndm();}

//---------------------  Simple exponential dist -------------------
  
ExponentialDist::ExponentialDist(double t): XeDist("Exponential") {setX0(t);}

       ExponentialDist::~ExponentialDist()   {}
void   ExponentialDist::setX0(double x)      {x0=x;}
double ExponentialDist::cdf(double x)        {return x<=0.? 0: 1.-exp(-x/x0);}
double ExponentialDist::pdf(double x)        {return exp(-x/x0)/x0;}
double ExponentialDist::logPdf(double x)     {return -x/x0-log(x0);}
double ExponentialDist::generate()           {return generate(x0);}
double ExponentialDist::generate(double x)   {return random.Exp(x);}


//---------------------  Simple Chi2 dist -------------------
  
Chi2Dist::Chi2Dist(): XeDist("Chi2") {setNdof(0);}
Chi2Dist::Chi2Dist(int dof): XeDist("Chi2 with "+formatI(dof)+" dof") {
  setNdof(dof);
}

       Chi2Dist::~Chi2Dist()          {}
void   Chi2Dist::setNdof(int n)       {ndof=n;}
double Chi2Dist::pdf(double x)        {return pdf(x,ndof);}
double Chi2Dist::cdf(double x)        {return below(x,ndof);}
double Chi2Dist::generate()           {return generate(ndof);}
double Chi2Dist::below(double x,int d){return 1.-above(x,d);}
double Chi2Dist::generate(int n)      {
  double c=0;
  for(int i=0;i<n;i++) {
    double x=GaussianDist::generate(0.,1.);
    c+= x*x;
  }
  return c;
}

double Chi2Dist::above(double x, int dof){
  if(x<=0.) return 1.;
  return ROOT::Math::chisquared_cdf_c(x,(double) dof);
}

double Chi2Dist::pdf(double x, int dof){
  if(x<=0.) return 0.;
  return ROOT::Math::chisquared_pdf(x,(double) dof);
}

double Chi2Dist::between(double x,double y,int dof){
  return above(x,dof)-above(y,dof);
}

      
//---------------------  Simple Gaussian dist ----------------

GaussianDist::GaussianDist(): XeDist("Normal") {setMu(0.);setSigma(1.);}

GaussianDist::GaussianDist(double m, double s): XeDist("Gaussian") {
  setMu(m);
  setSigma(s);
}

       GaussianDist::~GaussianDist() {}
void   GaussianDist::setMu(double m)    {mu=m;}
void   GaussianDist::setSigma(double s) {sigma=s;}
double GaussianDist::cdf(double x)      {return below((x-mu)/sigma);}
double GaussianDist::pdf(double x)      {return exp(logPdf(x));}
double GaussianDist::logPdf(double x)   {return logPdf(x,mu,sigma);}
double GaussianDist::generate()         {return generate(mu,sigma);}
double GaussianDist::above(double x)    {return 1.-below(x);}
double GaussianDist::below(double x)    {return ROOT::Math::normal_cdf(x);}

double GaussianDist::generate(double m, double s) {return random.Gaus(m,s);}
double GaussianDist::inside(double xmi, double xma, double sigma){
  return below(xma/sigma)-below(xmi/sigma);
}
double GaussianDist::logPdf(double x, double m, double s){
  return -(LOGSQR2PI+log(s)+(x-m)*(x-m)/2./s/s);
}

//---------------------  Simple Poisson dist -------------------
  
       PoissonDist::PoissonDist(double m): XeDist("Poisson") {setMu(m);}
       PoissonDist::~PoissonDist()              {}
void   PoissonDist::setMu(double m)             {mu=m;}
double PoissonDist::pdf(int n,double m)         {return pdf((double)n,m);}
double PoissonDist::pdf(double n,double m)      {return TMath::Poisson(n,m);}
double PoissonDist::logPdf(int n, double m)     {return logPdf((double)n,m);}
double PoissonDist::logPdf(double n, double m)  {return log(pdf(n,m));}
double PoissonDist::pdf(double x)               {return pdf(x,mu);}
double PoissonDist::logPdf(double x)            {return log(pdf(x,mu));}
double PoissonDist::generate()                  {return generate(mu);}
double PoissonDist::generate(double m)     {return m<=0.?0.:random.Poisson(m);}
double PoissonDist::aboveEqual(int n, double m) {return 1.-below(n,m);}
double PoissonDist::above(int n, double m)      {return 1.-belowEqual(n,m);}

double PoissonDist::cdf(double x)                { 
  return belowEqual((int)(x+1.e-9),mu);
}
double PoissonDist::belowEqual(int n, double m) {
  return ROOT::Math::poisson_cdf(n,m);
}
double PoissonDist::below(int n, double m){
  return n<=0? 0. : ROOT::Math::poisson_cdf(n-1,m);
}
double PoissonDist::between(int n1, int n2, double m){
  return 1.-below(n1,m)-above(n2,m);
}

void PoissonDist::fillTable(XeTable& table, double m, double epsilon){
  table.clear();
  double cum=0.;
  double low=epsilon/2.;
  double up=1.-low;
  for(int n=0;;n++){
    double p=pdf(n,m);
    cum+=p;
    if(cum>low) table.add(n,p);
    if(cum>up) break;
  }
}

//------------- (not so) abstract class for signal and background ----------

SignalAndBackground::~SignalAndBackground() {}
SignalAndBackground::SignalAndBackground(double sToE, double b, double db
                     ,double e, double de) :XeStat("Signal and background") {
  reset();
  setSigToEvents(sToE);
  setBackground(b);  
  setDBackground(db);  
  setEfficiency(e);
  setDEfficiency(de);
}

int    SignalAndBackground::getNObserved()           {return nObserved;}
double SignalAndBackground::getBackground()          {return background;}
double SignalAndBackground::getDBackground()         {return dBackground;}
double SignalAndBackground::getExposure()            {return exposure;}
double SignalAndBackground::getEfficiency()          {return efficiency;}
double SignalAndBackground::getDEfficiency()         {return dEfficiency;}
double SignalAndBackground::getSigToEvents()         {return sigToEvents;}
void   SignalAndBackground::setExposure(double e)    {exposure=e;}
void   SignalAndBackground::setEfficiency(double e ) {efficiency=e;}
void   SignalAndBackground::setDEfficiency(double dE){dEfficiency=dE;}
void   SignalAndBackground::setDBackground(double dB){dBackground=dB;}
void   SignalAndBackground::setSigToEvents(double s) {sigToEvents=s;}
void   SignalAndBackground::setNObserved(int ob)      {nObserved=ob;}
void   SignalAndBackground::setBackground(double b)  {background=b;}
void   SignalAndBackground::setNObservedAsBackground(){
  setNObserved((int)(nBackgroundEvents()+.5));
}

void SignalAndBackground::fillPoissonTable(XeTable& table,double e,double sig){
  PoissonDist::fillTable(table,nTotalEvents(sig),e);
}


void   SignalAndBackground::reset() {
  setExposure();
  setSigToEvents();
  setNObserved();
  setBackground();
  setEfficiency();
  setDEfficiency();
  setDBackground();
}

double SignalAndBackground::nSignalEvents(double s) {
  return exposure*s*sigToEvents*efficiency;
}
double SignalAndBackground::nBackgroundEvents()  {return exposure*background;}

double SignalAndBackground::nTotalEvents(double s) {
  return exposure*(s*sigToEvents*efficiency+background);
}

void   SignalAndBackground::printSBParameters()   {
  cout<<"   Cross section to evts: "<<sigToEvents<<endl
      <<"   Nominal efficiency   : "<<efficiency<<" +/- "<<dEfficiency<<endl
      <<"   Expected background  : "<<background<<" +/- "<<dBackground<<endl
      <<"   Exposure             : "<<exposure<<endl
      <<"   Observed events      : "<<nObserved<<endl;
}

double SignalAndBackground::probability(int obs,double sig){
  double o=(obs==AUTO)? nObserved:obs;
  double mu=nTotalEvents(sig);
  return PoissonDist::pdf(o,mu);
}

void  SignalAndBackground::combine(SignalAndBackground *sb1
                                  ,SignalAndBackground *sb2){
  vector<SignalAndBackground*> vsb;
  vsb.push_back(sb1);
  vsb.push_back(sb2);
  combine(vsb);
}

void  SignalAndBackground::combine(vector<SignalAndBackground*> &vsb){
  reset();
  unsigned int l=vsb.size();
  double sv=0.;
  double sb=0.;
  double sdb=0.;
  for(unsigned int i=0;i<l;i++){
    SignalAndBackground* s=vsb[i];
    double w=s->getExposure();
    if(s->getDEfficiency()!=0.){
      cout<<"Can't combine signalAndBackground if dEff!=0"<<endl;
      continue;
    } 
    sv+=   w*s->getEfficiency()*s->getSigToEvents();
    sb+=   w*s->getBackground();
    sdb+=  pow(w*s->getDBackground(),2);
  }
  setSigToEvents(sv);
  setBackground(sb);
  setDBackground(pow(sdb,.5));
}

//---------------------------- Confidence Intervals -------------------------
 
       CI::~CI()                  {}
bool   CI::withUpperLimit()       {return withUpper;}
bool   CI::withLowerLimit()       {return withLower;}
bool   CI::isInside(double x)     {return x>=LowerLimit && x<=UpperLimit;}
double CI::getUpperLimit()        {return UpperLimit;}
double CI::getLowerLimit()        {return LowerLimit;}

bool CI::withCLs(int m)  {
  switch(m){
     case CI_TWO_SIDED: case CI_UP : case CLS_UP : return false;
  }
  return true;
}

bool CI::withLowerLimit(int m){
  switch(m){
     case UNDEFINED_INT: case CI_UP : case CLS_UP : return false;
  }
  return true;
}

bool CI::withUpperLimit(int m){
  switch(m){
     case UNDEFINED_INT: case CI_LOW: case CLS_LOW: return false;
  }
  return true;
}

bool CI::isSingleLowerLimit(int m){
  switch(m){
     case CI_LOW : case CLS_LOW : return true;
  }
  return false;
}

bool CI::isSingleUpperLimit(int m){
  switch(m){
     case CI_UP  : case CLS_UP  : return true;
  }
  return false;
}

CI::CI(int m) :solvable(0,VERY_LARGE) , XeStat("CI") {
  mode=m;
  withLower=withLowerLimit(mode);
  withUpper=withUpperLimit(mode);
  applyCLs(withCLs(mode));
}

bool   CI::isCLsApplied()        {return doApplyCLs;}
void   CI::applyCLs(bool doIt)   {doApplyCLs=doIt;}

// ------------ Poisson confidence intervals ---------------

PoissonCI::PoissonCI(double sigToE, double bkg, int mod) : CI(mod){
  SB=new SignalAndBackground(sigToE, bkg);
  deleteSB=true;
  LowerLimit=0.;
  UpperLimit=VERY_LARGE;
}

PoissonCI::PoissonCI(SignalAndBackground *sb,int mod) : CI(mod){
  setSignalAndBackground(sb);
  LowerLimit=0.;
  UpperLimit=VERY_LARGE;
}

void PoissonCI::setSignalAndBackground(SignalAndBackground *sb) {
  SB=sb;
  deleteSB=false;
}

     PoissonCI::~PoissonCI()               {if(deleteSB) delete SB;}
void PoissonCI::setBackground(double b)    {SB->setBackground(b);}
void PoissonCI::setExposure(double e)      {SB->setExposure(e);}
void PoissonCI::setSigToEvents(double v)   {SB->setSigToEvents(v);}
void PoissonCI::setNObserved(int observed) {SB->setNObserved(observed);}
void PoissonCI::setNObservedAsBackground() {SB->setNObservedAsBackground();}
void PoissonCI::setEfficiency(double e)    {SB->setEfficiency(e);}
bool PoissonCI::expectsBackground()        {return expectsBackground(mode);}
SignalAndBackground* PoissonCI::getSignalAndBackground() {return SB;}

bool PoissonCI::expectsBackground(int m)  {
  switch(m){
     case CLS_LOW: case CLS_UP: case CLS_TWO_SIDED: return true;
  }
  return false;
}

void PoissonCI::computeLimits(int nObs, double cl) {
  setNObserved(nObs);
  computeLimits(cl);
}

void PoissonCI::computeLimits(double cl) {
  LowerLimit=0.;
  UpperLimit=VERY_LARGE;
  int nObserved=SB->getNObserved();
  switch(mode) {
     case CI_UP:
     case CLS_UP:
          UpperLimit=oneLimit(mode,cl);
          break;
     case CI_LOW:
     case CLS_LOW:
          LowerLimit=oneLimit(mode,cl);
          break;
     case CI_TWO_SIDED:
          UpperLimit=oneLimit(CI_UP,.5+.5*cl);
          LowerLimit=oneLimit(CI_LOW,.5+.5*cl);
          break;
     case CLS_TWO_SIDED:
          UpperLimit=oneLimit(CLS_UP,.5+.5*cl);
          LowerLimit=oneLimit(CLS_LOW,.5+.5*cl);
  } 
  if(printLevel>0) {
    cout<<"Limits for nOBs="<<nObserved<<", CL="<<cl<<" are ["<<LowerLimit
        <<","<<UpperLimit<<"]"<<endl;
  }
}

double PoissonCI::computeUpperSigma(double cl){
  int m=isCLsApplied()? CLS_UP : CI_UP;
  return oneLimit(m,cl);
}

double PoissonCI::oneLimit(int cMode,double cl){
  int nObserved=SB->getNObserved();
  if(nObserved==0 && (cMode==CLS_LOW || cMode==CI_LOW)) return 0.;
  currentMode=cMode;
  double eventsFor1cm2=SB->nSignalEvents(1.);
  double xmi=0.0001/eventsFor1cm2;
  double xma=max(20,4*nObserved)/eventsFor1cm2;
  setSolverLimits(xmi,xma);
  if(printLevel>1) {
     cout<<"Events for 1 cm^2"<<eventsFor1cm2<<", setting limits at "<<xmi
         <<" and "<<xma<<endl;
  }
  double limit=solve(1.-cl);
  if(isSingleLowerLimit(currentMode) && limit==xmi) limit=0.;
  return limit;
}


double PoissonCI::getValue(double sigma){
  double t=SB->nTotalEvents(sigma);
  double b=SB->nBackgroundEvents();
  int nObserved=SB->getNObserved();
  switch(currentMode) {
    case CLS_UP : 
      return PoissonDist::belowEqual(nObserved,t)
           / PoissonDist::belowEqual(nObserved,b);
    case CLS_LOW: 
      return PoissonDist::aboveEqual(nObserved,t)
           / PoissonDist::belowEqual(nObserved,b);
    case CI_UP  : 
      return PoissonDist::belowEqual(nObserved,t);
    case CI_LOW : 
      return PoissonDist::aboveEqual(nObserved,t);
  }
  cout<<"Invalid mode in PoissonCI::getValue() "<<currentMode<<endl;
  return UNDEFINED;
}

double PoissonCI::upperLimit(int m, int nObs, double bkg, double CL){
  PoissonCI ci(1.,bkg,m);
  if(ci.withUpperLimit()){
    ci.computeLimits(nObs,CL);
    return ci.getUpperLimit();
  }
  cout<<"Given PoissonCI mode doesn't compute upperLimits"<<endl;
  return UNDEFINED;
}

double PoissonCI::lowerLimit(int m, int nObs, double bkg, double CL){
  PoissonCI ci(1.,bkg,m);
  if(ci.withLowerLimit()){
    ci.computeLimits(nObs,CL);
    return ci.getLowerLimit();
  }
  cout<<"Given PoissonCI mode doesn't compute lowerLimits"<<endl;
  return UNDEFINED;
}

double PoissonCI::coverage(double sigma, double CL){
  int  lowestAcceptable=-1;
  int  upmostAcceptable=-1;
  double total=SB->nTotalEvents(sigma);
  int rmax=max(50,int(4.*total));
  for(int r=0;r<rmax;r++){
    computeLimits(r,CL); 
    if(isInside(sigma)){
      if(printLevel>1) cout<<"Accepted "<<r<<endl;
      if(lowestAcceptable<0) lowestAcceptable=r;
      upmostAcceptable=r;
    }
    else {
      if(printLevel>1) cout<<"Rejected "<<r<<endl;
      if(lowestAcceptable>=0 || upmostAcceptable>=0) {
        if(printLevel>1) cout<<"We can give up!"<<endl;
        break;
      }
    }
  } 
  double c=0.;
  switch(mode) {
    case CI_TWO_SIDED:
    case CLS_TWO_SIDED:
         c=PoissonDist::between(lowestAcceptable,upmostAcceptable,total);
         break;  
    case CI_UP:
    case CLS_UP:
         c=PoissonDist::aboveEqual(lowestAcceptable,total);
         break;  
    case CI_LOW:
    case CLS_LOW:
         c=PoissonDist::belowEqual(upmostAcceptable,total);
  }
  if(printLevel>1) {
    cout<<"Between "<<lowestAcceptable<<" and "<<upmostAcceptable
        <<": "<<c<<endl;
  }
  return c;
}

//-------------  Upperlimit ----------------------------------------

XeSigma::~XeSigma() {}

XeSigma::XeSigma() : XeObject("Cross section") {}

XeSigma::XeSigma(double ma,string eTitle)
  : XeObject(getTheName(ma,eTitle)) {mass=ma;}

XeSigma::XeSigma(double ma,double sig, double evt, string nam) : XeObject() {
  set(ma,sig,evt,nam);
}

void XeSigma::set(double ma, double sig,double evt,string nam){
  setMassAndName(ma,nam);
  setValues(sig,evt);
}

void XeSigma::setMassAndName(double ma, string nam){
  setName(getTheName(ma,nam));
  mass=ma;
}
void XeSigma::setValues(double sig,double evt){ sigma=sig; events=evt; }

string XeSigma::getTheName(double mass,string eTitle){
  return addASpace(eTitle)+"cross section for "+formatF(mass,6,1);
}

double XeSigma::getMass()          {return mass;}
double XeSigma::getSigma()         {return sigma;}
double XeSigma::getEvents()        {return events;}
double XeSigma::getValue(int unit) {
  switch(unit) {
    case SIGMA_UNIT :  return sigma;
    case EVENT_UNIT : return events;
  }
  cout<<"Unknown unit :"<<unit<<endl;
  return UNDEFINED;
}

bool XeSigma::printIt(int level){
  if(level<1) return true;
  if(level==1) {
    cout<<setprecision(6)<<mass<<" "<<sigma<<" "<<events
        <<RESET_PRECISION<<endl;
  }
  else{
    cout<<getName()<<" "<<setprecision(6)<<sigma<<" cm^2, equivalent to "
        <<events<<" events"<<RESET_PRECISION<<endl;
  }
  return true;
}

XeLimit::~XeLimit() {}
XeLimit::XeLimit() : XeStat("Upper Lmit"){}

XeLimit::XeLimit(double mass,string nam) : XeStat(getTheName(mass,nam)) {
 setMassAndName(mass,nam);
}

void XeLimit::setMassAndName(double mass,string nam){
  estimated.set(mass,0.,0.,nam);
  limit.set(mass,0.,0.,nam);
}

XeLimit::XeLimit(double mass,double estimatedSigma
     , double estimatedEvents,double limitSigma,double limitEvents,string nam)
   : XeStat(getTheName(mass,nam)) {
   setMassAndName(mass,nam);
   setLimits(estimatedSigma,estimatedEvents,limitSigma,limitEvents);
}

void XeLimit::setLimits( double estimatedSigma, double estimatedEvents
               , double limitSigma,     double limitEvents){
  estimated.setValues(estimatedSigma,estimatedEvents);
  limit.setValues(limitSigma,limitEvents);
}

string XeLimit::getTheName(double mass,string eTitle){
  return "Limit of "+addASpace(eTitle)+", mass= "+formatF(mass,6,1);
}

double XeLimit::getMass()              {return limit.getMass();}
double XeLimit::getUpperSigma()        {return limit.getSigma();}
double XeLimit::getUpperEvents()       {return limit.getEvents();}
double XeLimit::getUpperLimit(int unit){return limit.getValue(unit);}
double XeLimit::getEstimatedSigma()    {return estimated.getSigma();}
double XeLimit::getEstimatedEvents()   {return estimated.getEvents();}
double XeLimit::getEstimated(int unit) {return estimated.getValue(unit);}

double XeLimit::getValue(int mode,int unit){
  switch(mode) {
    case ESTIMATED       : return getEstimated(unit);
    case UPPER_LIMIT     : return getUpperLimit(unit);
  }
  cout<<"Invalid mode "<<mode<<endl;
  return UNDEFINED;
}

bool XeLimit::printIt(int level){
  if(level<1) return true;
  if(level==1) {
    cout<<formatF(getMass(),12,5)
        <<" "<<formatG(getEstimatedSigma(),15,6)
        <<" "<<formatG(getEstimatedEvents(),12,6)
        <<" "<<formatG(getUpperSigma(),15,6)
        <<" "<<formatG(getUpperEvents(),12,6)
        <<endl;
  }
  else{
    cout<<getName()<<" :"<<endl<<setprecision(6)
        <<"   Estimated: "<<getEstimatedSigma() <<" cm^2, equivalent to "
                          <<getEstimatedEvents()<<" events"<<endl
        <<"   Limits:    "<<getUpperSigma() <<" cm^2, equivalent to "
                          <<getUpperEvents()<<" events"<<endl<<RESET_PRECISION;
  }
  return true;
}


//------------  Table of cross sections -------------------
XeSigmas::XeSigmas() : map<double, XeSigma*>() , XeStat("Cross Sections") {}

XeSigmas::XeSigmas(string runName) : 
    map<double, XeSigma*>() , XeStat("Cross sections "+runName) {}

XeSigmas::~XeSigmas() {
  for (SigmaIterator it=begin();it!=end();it++) delete it->second;
}


void XeSigmas::read(string fName){
  clear();
  ifstream myfile((getResultsDirectory()+fName).c_str());
  if(myfile.is_open()) {
    string line;
    while (getline(myfile,line)) {
      vector<double> v=decodeVector(line);
      int n=v.size();
      if(n==0) {}
      else if(n==2) add(v[0],v[1],UNDEFINED);
      else if(n==3) add(v[0],v[1],v[2]);
      else cout<<"Wrong simple Sigma input line: "<<line<<endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file "<<fName; 
}

void XeSigmas::add(double mass, double sigma, double event){
  if(find(mass)==end()) {
    (*this)[mass]=new XeSigma(mass,sigma, event,getName());
    attach((*this)[mass]);
  }
  else cout<<"Can't add already existing sigma for mass "<<mass<<endl;
}

void XeSigmas::add(XeSigma *limit){
  double mass=limit->getMass();
  if(find(mass)==end()) {
    (*this)[mass]=limit;
    attach((*this)[mass]);
  }
  else cout<<"Can't add already existing sigma for mass "<<mass<<endl;
}


XeGraph* XeSigmas::newSigmaGraph(int plot)  {return newGraph(SIGMA_UNIT,plot);}
XeGraph* XeSigmas::newEventsGraph(int plot) {return newGraph(EVENT_UNIT,plot);}

XeGraph* XeSigmas::newGraph(int unit, int plot) {
  int n=size();
  vector<double> masses(n);
  vector<double> limits(n);
  int i=0;
  for(SigmaIterator it=begin();it!=end();it++) {
    masses[i]=it->first;
    XeSigma* upper=it->second;
    limits[i]=upper->getValue(unit);
    i++;
  }
  XeGraph* g=new XeGraph(getName(),n,&(masses[0]),&(limits[0]),LABEL_MASS
                        ,getSigmaLabel(unit)); 
  g->setDefaultXScale(LOG);
  g->setDefaultYScale(getSigmaLinLog(unit));
  g->setLineColor(kBlue);
  g->drawIt(plot);
  return g;
}

bool XeSigmas::printIt(int level){
  bool ok=true;
  if(level<1) return true;
  cout<<getName()<<endl;
  for(SigmaIterator it=begin();it!=end();it++) {
    cout<<"  ";
    XeSigma* sigma=it->second;
    ok=ok && sigma->printIt(level);
  }
  return ok;
}

//------------  Table of limit cross sections -------------------
XeLimits::XeLimits() : map<double, XeLimit*>() , XeStat("Upper Limits") {}

XeLimits::XeLimits(string runName) : 
    map<double, XeLimit*>() , XeStat("Upper Limits "+runName) {}

XeLimits::~XeLimits() {
  for (UpperLimitIterator it=begin();it!=end();it++) delete it->second;
}


void XeLimits::read(string fName){
  clear();
  ifstream myfile((getResultsDirectory()+fName).c_str());
  if(myfile.is_open()) {
    string line;
    while (getline(myfile,line)) {
      vector<double> v=decodeVector(line);
      int n=v.size();
      if(n==0) {}
      else if(n==2) add(v[0],UNDEFINED,UNDEFINED,v[1],UNDEFINED);
      else if(n==3) add(v[0],UNDEFINED,UNDEFINED,v[1],v[2]);
      else if(n==5) add(v[0],v[1],v[2],v[3],v[4]);
      else cout<<"Wrong  limits input line: "<<line<<endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file "<<fName; 
}

void XeLimits::add( double mass, double estimatedSigma, double estimatedEvents
                  , double sigmaLimit, double eventsLimit ){
  if(find(mass)==end()) {
    (*this)[mass]=new XeLimit(mass,estimatedSigma, estimatedEvents
                 ,sigmaLimit,eventsLimit ,getName());
    attach((*this)[mass]);
  }
  else cout<<"Can't add already existing limit for mass "<<mass<<endl;
}

void XeLimits::add(XeLimit *limit){
  double mass=limit->getMass();
  if(find(mass)==end()) {
    (*this)[mass]=limit;
    attach((*this)[mass]);
  }
  else cout<<"Can't add already existing limit for mass "<<mass<<endl;
}


XeGraph* XeLimits::newGraphOfEstimatedSigma(int plot)    {
  return newGraph(ESTIMATED, SIGMA_UNIT,plot);
}
XeGraph* XeLimits::newGraphOfEstimatedEvents(int plot)   {
  return newGraph(ESTIMATED, EVENT_UNIT,plot);
}

XeGraph* XeLimits::newGraphOfUpperSigma(int plot)    {
  return newGraph(UPPER_LIMIT, SIGMA_UNIT,plot);
}
XeGraph* XeLimits::newGraphOfEventLimit(int plot)   {
  return newGraph(UPPER_LIMIT, EVENT_UNIT,plot);
}


XeGraph* XeLimits::newGraph(int mode,int unit, int plot) {
  int n=size();
  vector<double> masses(n);
  vector<double> limits(n);
  int i=0;
  for(UpperLimitIterator it=begin();it!=end();it++) {
    masses[i]=it->first;
    XeLimit* upper=it->second;
    limits[i]=upper->getValue(mode,unit);
    i++;
  }
  string title=getSigmaModeName(mode)+" "+getName()
              +" ("+getSigmaUnitName(unit)+")";
  XeGraph* g=new XeGraph(title,n,&(masses[0]),&(limits[0]),LABEL_MASS
                ,getUpperSigmaLabel(unit)); 
  g->setDefaultXScale(LOG);
  g->setDefaultYScale(getSigmaLinLog(unit));
  g->setLineColor(kBlue);
  g->drawIt(plot);
  return g;
}

bool XeLimits::printIt(int level){
  bool ok=true;
  if(level<1) return true;
  cout<<getName()<<endl;
  for(UpperLimitIterator it=begin();it!=end();it++) {
    cout<<"  ";
    XeLimit* upper=it->second;
    ok=ok && upper->printIt(level);
  }
  return ok;
}



// ---------------- Sensitivity bands -------------------------

XeSensitivity::~XeSensitivity() {}

XeSensitivity::XeSensitivity() : XeStat("Sensitivity") {}
 
XeSensitivity::XeSensitivity(double ma, double e,double* lims,string eTitle)
  : XeStat(getTheName(ma,eTitle)) {
  setLimits(lims);
  signalPerCm2=e;
  mass=ma;
}

XeSensitivity::XeSensitivity(double ma,double e,string eTitle)
  : XeStat(getTheName(ma,eTitle)) {mass=ma; signalPerCm2=e;}

string XeSensitivity::getTheName(double mass,string eTitle){
  return addASpace(eTitle)+"sensitivity for "+formatF(mass,6,1);
}

double   XeSensitivity::getMass()                       {return mass;}
double   XeSensitivity::getSignalPerCm2()               {return signalPerCm2;}
double * XeSensitivity::getLimits()                     {return limits;}
void     XeSensitivity::setLimit(int mode,double limit) {limits[mode]=limit;}
double   XeSensitivity::getLimit(int l,int unit) {
  double x= limits[l];
  if(unit==EVENT_UNIT) x*=signalPerCm2;
  return x;
}
void     XeSensitivity::setLimits(double* lims) {
  for(int l=0;l<N_SENSITIVITY_MODES;l++) limits[l]=lims[l];
}
bool     XeSensitivity::printIt(int level) {
  if(level<1) return true;
  cout<<setprecision(4)<<mass<<" "<<signalPerCm2;
  for(int l=0;l<N_SENSITIVITY_MODES;l++) cout<<" "<<limits[l];
  cout<<RESET_PRECISION<<endl;
  return true;
}

SensitivityBands::SensitivityBands() : 
    map<double, XeSensitivity*>() , XeStylized("Sensitivity bands") {}

SensitivityBands::SensitivityBands(string runName) : 
    map<double, XeSensitivity*>(),XeStylized("Sensitivity bands of "+runName){}

SensitivityBands::~SensitivityBands() {
  for (SensitivityIterator it=begin();it!=end();it++) delete it->second;
}


void SensitivityBands::read(string fName){
  clear();
  ifstream myfile((getResultsDirectory()+fName).c_str());
  if(myfile.is_open()) {
    string line;
    while (getline(myfile,line)) {
      vector<double> v=decodeVector(line);
      int n=v.size();
      if(n==0) {}
      else if(n==2+N_SENSITIVITY_MODES) add(v[0],v[1],&(v[2]));
      else cout<<"Wrong Sensitivity input line: "<<line<<endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file "<<fName; 
}

void SensitivityBands::add(double mass, double ePerCm2, double* limits){
  if(find(mass)==end()) {
    (*this)[mass]=new XeSensitivity(mass,ePerCm2,limits,getName());
    attach((*this)[mass]);
  }
  else cout<<"Can't add already existing limit for mass "<<mass<<endl;
}

void SensitivityBands::add(XeSensitivity *sensitivity){
  double mass=sensitivity->getMass();
  if(find(mass)==end()) {
     (*this)[mass]=sensitivity;
     attach((*this)[mass]);
  }
  else cout<<"Can't add already existing limit for mass "<<mass<<endl;
}

string SensitivityBands::bandName(int which){
  switch(which) {
    case MINUS_TWO_SIGMAS : return "-2 \\sigma";
    case MINUS_ONE_SIGMA  : return "-1 \\sigma";
    case MEDIAN           : return "median";
    case PLUS_ONE_SIGMA   : return "+1 \\sigma";
    case PLUS_TWO_SIGMAS  : return "+2 \\sigma";
  }
  cout<<"Unknown band type :"<<which<<endl;
  return UNDEFINED_STRING;
}

int SensitivityBands::lineStyle(int which){
  int lines[] = { DOTTED, DASH_DOTTED, PLAIN, DASH_DOTTED, DOTTED };
  return lines[which]; 
}

double SensitivityBands::numberOfSigmas(int which){
  double sigmas[] = { -2., -1., 0., 1., 2. };
  return sigmas[which]; 
}


XeGraph* SensitivityBands::newGraph(int which,int unit){
  int n=size();
  vector<double> masses(n);
  vector<double> limits(n);
  int i=0;
  for(SensitivityIterator it=begin();it!=end();it++) {
    masses[i]=it->first;
    XeSensitivity* sens=it->second;
    limits[i]=sens->getLimit(which,unit);
    i++;
  }
  string leg=bandName(which);
  XeGraph* g=new XeGraph(getName(),n,&(masses[0]),&(limits[0]),LABEL_MASS
                ,LABEL_SIGMA_MAX,NULL,leg); 
  g->setLineStyle(lineStyle(which));
  g->setDefaultXScale(LOG);
  g->setDefaultYScale(getSigmaLinLog(unit));
  return g;
}

XeMultiGraph* SensitivityBands::newMultiGraph(int unit){
  XeMultiGraph *mg=new XeMultiGraph(getName());
  for(int m=0;m<N_SENSITIVITY_MODES;m++) mg->add(newGraph(m,unit));
  return mg;
}

void SensitivityBands::draw(string ,int unit ){
  STATIC_CONST  int nBands=3;
  int colors[nBands]={5,8,5};
  int low[nBands]={MINUS_TWO_SIGMAS,MINUS_ONE_SIGMA,PLUS_ONE_SIGMA}; 
  int high[nBands]={MINUS_ONE_SIGMA,PLUS_ONE_SIGMA,PLUS_TWO_SIGMAS};
  TPolyLine *polygons[nBands];
  for(int i=0;i<nBands;i++){
    polygons[i]=newPolygon(low[i],high[i],unit);
    polygons[i]->SetLineColor(0);
    polygons[i]->SetFillColor(colors[i]);
    polygons[i]->Draw("F");
  } 
}

TPolyLine* SensitivityBands::newPolygon(int low, int high,int unit){
  int n=size();
  vector<double> x(2*n+1); 
  vector<double> y(2*n+1);
  int i=0; 
  for(SensitivityIterator it=begin();it!=end();it++) {
    XeSensitivity* sens=it->second;
    double mass=it->first;
    x[i]=mass;
    x[2*n-1-i]=mass;
    y[i]=sens->getLimit(high,unit);
    y[2*n-1-i]=sens->getLimit(low,unit);
    i++;
  }
  x[2*n]=x[0];
  y[2*n]=y[0];
  return new TPolyLine(2*n+1,&x[0],&y[0]);
}



double  SensitivityBands::getMinX() {return begin()->first;}
double  SensitivityBands::getMaxX() { 
  SensitivityIterator it=end()--; 
  return it->first;
}
double  SensitivityBands::getMaxY(){
  double m=0.;
  for(SensitivityIterator it=begin();it!=end();it++) {
    XeSensitivity* sens=it->second;
    m=max(m,sens->getLimit(PLUS_TWO_SIGMAS));
  }
  return m;
}

double  SensitivityBands::getMinY(){
  double m=VERY_LARGE;
  for(SensitivityIterator it=begin();it!=end();it++) {
    XeSensitivity* sens=it->second;
    m=min(m,sens->getLimit(MINUS_TWO_SIGMAS));
  }
  return m;
}

double  SensitivityBands::getMinYNotZero(){
  double m=VERY_LARGE;
  for(SensitivityIterator it=begin();it!=end();it++) {
    XeSensitivity* sens=it->second;
    double l=sens->getLimit(MINUS_TWO_SIGMAS);
    if(l>0.) m=min(m,l);
  }
  return m;
}

bool SensitivityBands::printIt(int level){
  bool ok=true;
  if(level<1) return true;
  for(SensitivityIterator it=begin();it!=end();it++) {
    XeSensitivity* sens=it->second;
    ok=ok && sens->printIt();
  }
  return ok;
}

// ---------------- Data Sets, real and toys ------------------

DataSet::DataSet() : XeStat()  {}
DataSet::DataSet(int nCol)              : XeStat()     {setColumns(nCol);}
DataSet::DataSet(string name, int nCol) : XeStat(name) {setColumns(nCol);}

void DataSet::setColumns(int nCol){
  nColumns=nCol;
  entries.resize(nColumns);
  nEvents=0;
  for(int c=0;c<nColumns;c++) entries[c]=new vector<double>();
}

DataSet::~DataSet(){clear();}

void DataSet::clear(){
  for(int c=0;c<nColumns;c++)  entries[c]->clear();
  nEvents=0;
}

// first index is event, second is parameter

bool    DataSet::update()  {
  if(printLevel>0) cout<<"Default void udpdate for "<<getName()<<endl;
  return true;
}

int     DataSet::getNEvents()                 {return nEvents;}
int     DataSet::getNColumns()                {return nColumns;}
double  DataSet::getValue(int entry,int col)  {
  if(col<0 || col>=nColumns){
    cout<<"Wrong column in "<<getName()<<" : "<<col<<endl;
  }
  return (*entries[col])[entry];
}

vector<double>* DataSet::getColumn(int c)     {return entries[c];}
double*         DataSet::getVector(int c)     {return &((*entries[c])[0]);}

void DataSet::getEntry(int e,double *values)   {
  for(int c=0;c<nColumns;c++) values[c]=(*entries[c])[e];
}

void    DataSet::addEntry(double* x) {
  for(int c=0;c<nColumns;c++) entries[c]->push_back(x[c]);
  nEvents++;
}

double DataSet::getMax(int col){
  double m=-VERY_LARGE;
  int ne=getNEvents();
  for(int e=0;e<ne;e++) m=max(m,(*entries[col])[e]);
  return m;
}

double DataSet::getMin(int col){
  double m=VERY_LARGE;
  int ne=getNEvents();
  for(int e=0;e<ne;e++) m=min(m,(*entries[col])[e]);
  return m;
}

bool DataSet::printIt(int level){
  int ne=getNEvents();
  bool ok=ne>0;
  cout<<getName()<<" has "<<ne<<" entries, "<<nColumns<<" columns"<<endl;
  if(level<0) return ok;
  for(int e=0;e<ne;e++){
    cout<<setw(4)<<e<<" : ";
    for(int c=0;c<nColumns;c++) cout<<setw(10)<<getValue(e,c)<<" ";
    cout<<endl;
  }
  return ok>0;
}

//------------ Simulated data sets for package testing -------------------- 

SimulatedDataSet::~SimulatedDataSet(){}

SimulatedDataSet::SimulatedDataSet() : DataSet() {}
SimulatedDataSet::SimulatedDataSet(string name, int nCol) : DataSet(name,nCol) {
  transient.resize(nCol); 
}

void SimulatedDataSet::initialize(){}

void SimulatedDataSet::simulate(double muSigma, double muBackground){
  int nS=PoissonDist::generate(muSigma);
  int nB=PoissonDist::generate(muBackground);
  generate(nS,nB);
}

void SimulatedDataSet::generate(int nS, int nB){
  initialize();
  clear();  
  for(int i=0;i<nS;i++){
    generateSignal(&transient[0]);
    addEntry(&transient[0]);
  }
  for(int i=0;i<nB;i++){
    generateBackground(&transient[0]);
    addEntry(&transient[0]);
  }
  if(XeStat::getPrintLevel()>0){
    cout<<"Generating data for"<<getName()<<", signal="<<nS
        <<", background="<<nB<<endl;
  }
}

//--------------- One dimension Simulated data set ----------------------

OneDimSimulatedDataSet::~OneDimSimulatedDataSet() {}
OneDimSimulatedDataSet::OneDimSimulatedDataSet() : SimulatedDataSet("",1) {
  signal=NULL;
  background=NULL;
}


XeDist* OneDimSimulatedDataSet::getBackground()         {return background;}
XeDist* OneDimSimulatedDataSet::getSignal()             {return signal;}
void    OneDimSimulatedDataSet::setSignal(XeDist* s)    {signal=s;initialize();}
void    OneDimSimulatedDataSet::setBackground(XeDist* b){
  background=b;
  initialize();
}

void OneDimSimulatedDataSet::generateSignal(double *x){
  *x=signal->generate();
}

void OneDimSimulatedDataSet::generateBackground(double *x){
  *x=background->generate();
}

void OneDimSimulatedDataSet::initialize(){
  if(signal==NULL || background==NULL) return;
  setName("Simulated "+signal->getName()+" + "+background->getName());
}

//-----------  Simulated Exponential signal + flat background -------------

SimulatedExponentAndBackground::~SimulatedExponentAndBackground() {
  delete exp;
  delete uni;
}

SimulatedExponentAndBackground::SimulatedExponentAndBackground(double x
  , double b) : OneDimSimulatedDataSet() {
  exp=new ExponentialDist(x);
  signal=exp;
  uni=new UniformDist(0.,b);  
  background=uni;
}


void SimulatedExponentAndBackground::setX0(double x)     {exp->setX0(x);}
void SimulatedExponentAndBackground::setXmaxB(double x)  {uni->setLimits(0.,x);}


//---------- Simulated Gaussian + flat background ----------------------

SimulatedGaussianAndBackground::~SimulatedGaussianAndBackground() {
  delete gau;
  delete uni;
}

SimulatedGaussianAndBackground::SimulatedGaussianAndBackground(double mu
  ,double sigma, double b) : OneDimSimulatedDataSet() {
  gau=new GaussianDist(mu,sigma);
  signal=gau;
  uni=new UniformDist(0.,b);  
  background=uni;
}

void SimulatedGaussianAndBackground::setMu(double x)     {gau->setMu(x);}
void SimulatedGaussianAndBackground::setSigma(double x)  {gau->setSigma(x);}
void SimulatedGaussianAndBackground::setXmaxB(double x)  {uni->setLimits(0.,x);}

//-------------------- Virtual class for calculating pValues ----

PValue::~PValue() {}
PValue::PValue() : XeStat() {initializeIt();}
PValue::PValue(string nam,int analMode) :  XeStat(nam) {
  setName(nam);
  initializeIt(analMode);
}

void PValue::initializeIt(int analMode){
  requestedAnalysisMode=analMode;
  initializedOK=false;
  inError=false;
  sigmaUnit=0.;
  lowerLimit=0.;
  upperLimit=DEFAULT_EVENT_UPPER_LIMIT;
}

bool PValue::isAnalysisModeOK(){
  return checkAnalysisMode(getName(),requestedAnalysisMode);
}

bool PValue::initialize(){
  if(!initializedOK) initializedOK=checkPValue();
  return initializedOK;
}

void   PValue::setSigmaUnit(double s)        {sigmaUnit=s;updateSigmaUnit();}
void   PValue::exclusionComputed()           {}
void   PValue::estimateCrossSection()        {}
void   PValue::setEventsLowerLimit(double l) {lowerLimit=l;}
void   PValue::setEventsUpperLimit(double l) {upperLimit=l;}
int    PValue::forceCLs()                    {return DO_NOT_FORCE;}
bool   PValue::isInError()                   {return inError;}
bool   PValue::isInitializedOK()             {return initializedOK;}
double PValue::pValueB()                     {return pValueS(0.);}
double PValue::eventsLowerLimit()            {return lowerLimit;}
double PValue::eventsUpperLimit()            {return upperLimit;}

bool   PValue::simulate(double ) {
  cout<<"No simulation implemented for "<<getName()<<endl;
  return false;
} 

double PValue::getSigmaHut() {
  cout<<"No sigma hat method implemented for "<<getName()<<endl;
  return UNDEFINED;
}

void   PValue::updateSigmaUnit() {
  if(printLevel>0){
    cout<<"Default 'do nothing' updateSigmaUnit"<<endl;
  }
}

double PValue::qExclusion(double) {
  cout<<"qExclusion isn't defined for this pValue!"<<endl;
  return UNDEFINED;
}

void   PValue::setWimpMass(double) {
  if(printLevel>0){
    cout<<"Default 'do nothing' setWimpMass"<<endl;
  }
}

//--------------  A simple PValue : Sigma and Background ---------------


PVcountingSB::~PVcountingSB(){}
PVcountingSB::PVcountingSB(double sigToE, double bkg, double eff) :
  PValue("S+B",CUTS_ANALYSIS), SignalAndBackground(sigToE,bkg,0.,eff,0.) {}


void   PVcountingSB::printFlagsAndParameters()   {
  if(isExperimentSet()) cout<<"Experiment "<<experiment <<"; ";
  printSBParameters();
}

double  PVcountingSB::pValueS(double sigma){
  double t=nTotalEvents(sigma);
  return PoissonDist::belowEqual(nObserved,t);
}   

void   PVcountingSB::setEventsUpperLimit(double){
  cout<<"Can't setEventsUpperLimit for "<<getName()<<endl;  
}

double PVcountingSB::eventsUpperLimit()  {
  return Exclusion::typicalUpperLimit(nObserved);
}

double PVcountingSB::getSigmaHut() {
  return (nObserved/exposure-background)/efficiency/sigToEvents;
}

double PVcountingSB::nSignalPerCm2()    {return nSignalEvents(1.);}
bool   PVcountingSB::checkPValue()      {return true;}
bool   PVcountingSB::update(){
  if(printLevel>0){
    cout<<"Sigma to events : "<<sigToEvents*exposure<<endl;
  }
  return true;
}

// -----------------  Exclusion finder --------------------------

Exclusion::~Exclusion(){}
Exclusion::Exclusion() : XeStat(), solvable(0.,1.)   {
  reset();
}
Exclusion::Exclusion(PValue *pV ) : XeStat(), solvable(0.,1.)   {
  reset();
  addPValue(pV);
  setName("Exclusion with "+pV->getName());
}

Exclusion::Exclusion(PValue** pV,int l ) : XeObject(), solvable(0.,1.)   {
 reset();
 if(printLevel>2){
   cout<<"Exclusion instantiated with "<<l<<" p-values"<<endl;
   for(int i=0;i<l;i++){
     PValue* p=pV[i];
     cout<<i<<" "<<p<<" "<<endl;
     cout<<p->getName()<<endl;
   }
 }
 addPValues(pV,l);
}

void Exclusion::reset(){
  setName("Exclusion");
  nPValues=0;
  wimpMass=0.;
  initializedOK=false;
  fisherCorrection=false;
  pValues.clear();
  applyCLs(true);
}

void Exclusion::setSigmaUnitAndLimits(double unit,double sMin, double sMax){
  setSolverLimits(sMin*unit,sMax*unit);
  if(printLevel>1) {
    cout<<"Setting search limits to ["<<sMin<<","<<sMax
        <<"] events = ["<<xmin<<","<<xmax<<"] cm^2"<<endl;
  }
  for(int i=0;i<nPValues;i++) pValues[i]->setSigmaUnit(unit);
}

bool   Exclusion::isCLsApplied()                    {return doApplyCLs;}
void   Exclusion::applyCLs(bool doIt)               {doApplyCLs=doIt;}
void   Exclusion::applyCLsAndSave(bool doIt)        {saveCLs();applyCLs(doIt);}
void   Exclusion::saveCLs()                         {savedApplyCLs=doApplyCLs;}
void   Exclusion::restoreCLs()                      {doApplyCLs=savedApplyCLs;}
void   Exclusion::applyFisherCorrection(bool doIt)  {fisherCorrection=doIt;}

void   Exclusion::addPValue(PValue *pV) {
  if(printLevel>1) cout<<"Adding "<<pV->getName()<<endl;
  if(!pV->isAnalysisModeOK()) return;
  else if(pV==NULL) setName(string("Exclusion by undefined pValue"));
  else {             
    pV->initialize();
    if(pV->isInError()){
      cout<<"Can't add pValue "<<pV->getName()<<"; it has an error"<<endl;
      return;
    }
    addName(pV->getName());
    pValues.push_back(pV);
    nPValues=pValues.size();
  }
}

double Exclusion::typicalUpperLimit(int n)    {
  if(n<5) return 20;
  else if(n<10) return 30;
  return 3*n;
}

void Exclusion::addPValues(PValue** pV,int l){
  if(printLevel>1)cout<<"Adding "<<l<<" pValues "<<endl;
  for(int i=0;i<l;i++) {
    if(printLevel>2)cout<<"Adding pValue "<<i<<" :"<<pV[i]->getName()<<endl;
    addPValue(pV[i]);
  }
}

bool Exclusion::initialize(){
  if(initializedOK) return true;
  if(!isAnalysisDefined(true))return false;
  if(nPValues==0) {
    cout<<"No pValue calculator defined for "<<getName()<<endl;
    return false;
  }
  int of=DO_NOT_FORCE;
  for(int i=0;i<nPValues;i++) {
    if(!pValues[i]->initialize()) return false;
    int f=pValues[i]->forceCLs();
    if(f==FORCE_TRUE || f==FORCE_FALSE){
      if(of==DO_NOT_FORCE) of=f;
      else if(of!=f) cout<<"****** Incompatibility in CLs forcing mode"<<endl;
    }
  }
  if(of!=DO_NOT_FORCE){
    bool cls=(of==FORCE_TRUE);
    if(cls!=doApplyCLs){
      doApplyCLs=cls;
      cout<<"CLs mode of "<<getName()<<" forced to "
           <<yesOrNo(doApplyCLs)<<endl;
    }
  }
  if(printLevel>1) {
    cout<<"Initialized "<<getName();
    if(doApplyCLs) cout<<" with CLs"<<endl;
    else           cout<<" without CLs"<<endl;
  }
  if(!update()) return false;
  initializedOK=true;
  return true;
}

double Exclusion::getValue(double sigma){
  if(printLevel>2) cout<<"Getting pValue for sigma= "<<sigma<<endl;
  return correctedPValue(sigma,true);
}

bool Exclusion::update() {
  for(int i=0;i<nPValues;i++){
    if(!pValues[i]->update()){
      cout<<"---> Something fishy in update mechanism !"<<endl;
      return false;
    }
  }
  return true;
}
 
bool Exclusion::prepareForLimit(){
  if(printLevel>1) {
    cout<<"Preparing "<<getName()<<" with "<<nPValues<<" p-values"<<endl;
  }
  if(!initialize()) {
    cout<<"---> Something fishy in the set-up!"<<endl;
    return false;
  }
  signalPerCm2=0.;
  eventMax=0.;
  eventMin=-100.;
  for(int i=0;i<nPValues;i++) {
    signalPerCm2+=pValues[i]->nSignalPerCm2();
    double upl=pValues[i]->eventsUpperLimit();
    double low=pValues[i]->eventsLowerLimit();
    if(printLevel>1) {
      cout<<"Events limit for "
         <<pValues[i]->getName()<<" low="<<low<<", up="<<upl<<endl;
    }
    eventMax=max(eventMax,upl);
    eventMin=max(eventMin,low);
  }
  if(printLevel>0) {
    cout<<"Expected number of events "<<signalPerCm2<<" for 1 cm^2, range "
        <<eventMin<<" to "<<eventMax<<endl;
  }
  if(signalPerCm2<=0.) {
    cout<<"Can't compute limit for given mass"<<endl;
    return false;
  }
  setSigmaUnitAndLimits(1.,eventMin,eventMax);
  //setSigmaUnitAndLimits(1./signalPerCm2,eventMin,eventMax);
  limit.setMassAndName(wimpMass,getName());
  return update();
}


XeSigma* Exclusion::computeCrossSection(){
  if(printLevel>0){
    cout<<endl
        <<"========== Computing cross section for Wimp mass: "
        <<setprecision(4)<<wimpMass<<RESET_PRECISION<<endl;
  }
  if(!prepareForLimit()) return NULL;
  for(int i=0;i<nPValues;i++) pValues[i]->estimateCrossSection();
  double estimatedSigma=pValues[0]->getSigmaHut();
  double estimatedEvents=estimatedSigma*signalPerCm2;
  if(printLevel>0) {
    cout<<" Sigma computed, estimated events:"<<estimatedEvents
        <<", "<<estimatedSigma<<" cm^2"<<endl;
  }
  return new XeSigma(wimpMass,estimatedSigma,estimatedEvents);
}

double Exclusion::computeUpperSigma(double cl){
  computeLimit(cl);
  return getUpperSigma();
}

double   Exclusion::getUpperSigma()       {return limit.getUpperSigma(); }
double   Exclusion::getUpperEvents()      {return limit.getUpperEvents(); }
double   Exclusion::getEstimatedSigma()   {return limit.getEstimatedSigma(); }
double   Exclusion::getEstimatedEvents()  {return limit.getEstimatedEvents(); }
XeLimit* Exclusion::getLimit()            {return &limit;}

void Exclusion::computeLimit(double cl){
  stopper* stop=NULL;
  if(printLevel>0){
    cout<<endl
        <<"========== Upper sigma cl="<<cl<<" for Wimp mass: "
        <<setprecision(4)<<wimpMass<<RESET_PRECISION<<endl;
    stop=new stopper("Compute limit");
  }
  if(!prepareForLimit()) return;
  if(printLevel>1) cout<<"Estimate cross section"<<endl;
  for(int i=0;i<nPValues;i++) pValues[i]->estimateCrossSection();
  savedPBkg=computePB();
  double estimatedSigma=pValues[0]->getSigmaHut();
             double estimatedEvents=estimatedSigma*signalPerCm2;
  if(printLevel>1) {
    cout<<" Sigma fitted, estimated events:"<<estimatedEvents
        <<", saved pValue for background: "<<savedPBkg<<endl;
  }

  if(fisherCorrection) cl=cl*(nPValues+1)/2./nPValues;
  double sigmaUp=solve(1.-cl);
  double eventsUp=signalPerCm2*sigmaUp;
  for(int i=0;i<nPValues;i++) pValues[i]->exclusionComputed();
  limit.setLimits(estimatedSigma, estimatedEvents,sigmaUp,eventsUp);
  if(printLevel>0) {
   stop->print();
   delete stop;
   limit.print(2); 
  }
}

void Exclusion::simulate(double sigma){
  for(int p=0;p<nPValues;p++) pValues[p]->simulate(sigma);
}

double Exclusion::qExclusion(double sig){
  double q=1.;
  for(int p=0;p<nPValues;p++) q *= pValues[p]->qExclusion(sig);
  return q; 
}



void Exclusion::simulateNQValues(int n, XeValues *qs, XeValues *sHats
                  ,XeValues *eHats, double sigSim , double sigTest){
  applyCLsAndSave(false);
  if(qs!=NULL) qs->clear();
  if(sHats!=NULL) sHats->clear();
  if(eHats!=NULL) eHats->clear();
  for(int i=0;i<n;i++){
    simulate(sigSim);
    if(!prepareForLimit()) return;
    for(int p=0;p<nPValues;p++) pValues[p]->estimateCrossSection();
    double q=qExclusion(sigTest);
    if(qs!=NULL) qs->add(q);
    if(sHats!=NULL) sHats->add(pValues[0]->getSigmaHut());
    if(eHats!=NULL) eHats->add(pValues[0]->getSigmaHut()*signalPerCm2);
  }
  restoreCLs();
}

XeGraph* Exclusion::newGraphOfQExclusion(XeRange *crossSection, int unit){
  XeRange* sr=newSigmaRange(crossSection,unit);
  XeGraph* gr=new XeGraph("Q exclusion of "+getName());
  gr->setLineColor(kBlue);
  gr->setXYLabels(getSigmaLabel(unit),LABEL_Q_EXCLUSION);
  int n=sr->getNPoints();
  gr->set(n);
  for(int i=0;i<n;i++){
    double cross=sr->getValue(i);
    if(unit==EVENT_UNIT) cross=cross/signalPerCm2;
    double q=1.;
    for(int p=0;p<nPValues;p++) pValues[p]->estimateCrossSection();
    q=qExclusion(cross);
    gr->setPoint(i,cross,q);
  }
  return gr;
}

void Exclusion::simulateNUpperLimits(int n, XeValues* sigmas, XeValues* events
                                    , double cl){
  for(int i=0;i<n;i++){
    simulate(0.);
    computeLimit(cl);
    if(sigmas!=NULL) sigmas->add(getUpperSigma());
    if(events!=NULL) events->add(getUpperEvents());
  }
}

XeSensitivity* Exclusion::newSensitivity(double mass, int nSimul, double cl){
  cout<<"Computing sensitity with "<<nSimul<<" simulations, WIMP mass="
      <<getWimpMass()<<" ... be patient!"<<endl;  
  XeValues sigmas;
  setWimpMass(mass);
  simulateNUpperLimits(nSimul,&sigmas,NULL,cl);
  XeSensitivity *sens=new XeSensitivity(mass,signalPerCm2,getName());
  for(int w=0;w<N_SENSITIVITY_MODES;w++) {
    sens->setLimit(w,sigmas.findSigma(SensitivityBands::numberOfSigmas(w)));
  }
  return sens;
}

SensitivityBands* Exclusion::newSensitivityBands(int nSimul, XeRange* mr 
                                                , double cl){
  SensitivityBands* bands=new SensitivityBands(getName());
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  for(int i=0;i<n;i++){
    double mass=mr->getValue(i);
    setWimpMass(mass);
    update();
    XeSensitivity *sens=newSensitivity(mass,nSimul,cl);
    bands->add(sens);
  }

  return bands;
}

double Exclusion::getWimpMass() {return wimpMass;}

void Exclusion::setWimpMass(double mass){
  if(printLevel>1) {
    cout<<"Exclusion finder set wimp Mass to "<<mass<<endl;
  }
  for(int i=0;i<nPValues;i++) pValues[i]->setWimpMass(mass);
  wimpMass=mass;
}

double Exclusion::correctedPValue(double sigma, bool useSaved){
  if(printLevel>1) {
    cout<<"Computing pValue for real sigma="<<sigma<<endl;
  }
  if(!initialize()) return 1.;
  double Ps=1.;
  for(int i=0;i<nPValues;i++) {
    double p=pValues[i]->pValueS(sigma);
    if(printLevel>1)  cout<<"      --> Ps="<<formatF(p,19,8)<<endl;
    if(p<0. || p>1.) {
      cout<<pValues[i]->getName()<<" returned a wrong value for Ps: "
           <<Ps<<endl;
      return 1.;
    }
    Ps*= p;
  }
  double Pb=1.;
  if(isCLsApplied()){
    if(!useSaved) savedPBkg=computePB();
    Pb=savedPBkg;
    if(printLevel>1)  cout<<"   --> Pb="<<formatF(Pb,19,8)<<endl;
    if(printLevel>1)  cout<<"   --> Ps/Pb="<<formatF(Ps/Pb,19,8)<<endl;
  }
  double ratio=Ps/Pb; 
  if(nPValues==1) return ratio;
  else return Chi2Dist::above(-2*log(ratio),2*nPValues);
}

double Exclusion::computePB(){

  double Pb=1.;
  for(int i=0;i<nPValues;i++) {
    double p=pValues[i]->pValueB();
    if(printLevel>1)  cout<<"      --> Pb="<<formatF(p,19,8)<<endl;
    if(p<0. || p>1.) {
      cout<<pValues[i]->getName()<<" returned a wrong value for Pb: "
          <<p<<endl;
      return 1.;
     }
     Pb *= p;
  }
  if(printLevel>1)  cout<<"   --> Pb="<<formatF(Pb,19,8)<<endl;
  return 0.5;   //Ale  FIXME pValueB = pvalS(0.) which returns always 0.
  //return Pb;
}

XeGraph* Exclusion::newGraphOfPValue(XeRange *sr,int plot){
  update();
  XeGraph* gr=new XeGraph("P values of "+getName());
  gr->setLineColor(kBlue);
  gr->setXYLabels(LABEL_SIGMA,LABEL_PVALUE);
  int n=sr->getNPoints();
  gr->set(n);
  for(int i=0;i<n;i++){
    double sigma=sr->getValue(i);
    double cp=correctedPValue(sigma,false);
    gr->setPoint(i,sigma,cp);
  }
  gr->drawIt(plot);
  return gr;
}

XeGraph* Exclusion::newGraphOfUpperSigma(XeRange* mr,int plot, double cl){
  return newGraphOfUpperLimit(SIGMA_UNIT,mr,cl,plot);
}

XeGraph* Exclusion::newGraphOfUpperEvents(XeRange* mr,int plot, double cl){
  return newGraphOfUpperLimit(EVENT_UNIT,mr,cl,plot);
}

XeGraph* Exclusion::newGraphOfUpperLimit(int unit, XeRange* mr,int plot
 , double cl){
  XeLimits* limits=newUpperLimits(mr,cl);
  XeGraph *g=limits->newGraph(UPPER_LIMIT,unit);
  delete limits;
  g->drawIt(plot);
  return g;
}

XeLimits* Exclusion::newUpperLimits(XeRange* mr,double cl){
  if(mr==NULL) mr=XeRange::getDefaultMassRange();
  int n=mr->getNPoints();
  XeLimits* limits=new XeLimits(getName());
  for(int i=0;i<n;i++){
    double mass=mr->getValue(i);
    setWimpMass(mass);
    computeLimit(cl);
    limits->add(new XeLimit(limit));
  }
  return limits;
}


//--------------- Parameters of a likelihood computations ------------

LKParameter::~LKParameter(){}
LKParameter::LKParameter() : XeStat() { t0=0.;}
LKParameter::LKParameter(int i, int typ,string nam,double initV,double st
                    ,double mi,double ma) : XeStat() {
  if(typ<0 || typ>=N_PARAMETER_TYPES){
    cout<<"Unknown LK parameter type "<<typ<<endl;
    return ;
  }

  if(printLevel>1) cout <<"LKParameter::LKParameter()  " << nam << " Enter"<< endl;

  t0 = 0.;
  setName(nam);
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

  if(printLevel>1) cout <<"LKParameter::LKParameter()  " << nam << " Exit" << endl;
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
  if(printLevel>1) cout<<"Minuit Scaling of "<<Name<<": "<<s<<endl;
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
  if(Name!=getName()) {
    if(prnt) cout<<"Mismatch in name "<<Name<<" vs "
                  <<par->getName()<<endl;
    return false;
  }
  if(id!=par->getId()) {
    if(prnt) cout<<"Mismatch in id "<<id<<" vs "<<par->getId()
                  <<" for "<<Name<<endl;
    return false;
  }
  if(type!=par->getType()) {
    if(prnt) cout<<"Mismatch in type "<<type<<" vs "<<par->getType()
                  <<" for "<<Name<<endl;
    return false;
  }
  if(initialValue!=par->getInitialValue()) {
    if(prnt) cout<<"Mismatch in initial value "<<initialValue<<" vs "
                  <<par->getInitialValue()<<" for "<<Name<<endl;
    return false;
  }
  if(step!=par->getStep()) {
    if(prnt) cout<<"Mismatch in step "<<step<<" vs "<<par->getStep()
                  <<" for "<<Name<<endl;
    return false;
  }
  if(minimum!=par->getMinimum()) {
    if(prnt) cout<<"Mismatch in minimum "<<minimum<<" vs "
                  <<par->getMinimum()<<" for "<<Name<<endl;
    return false;
  }
  if(maximum!=par->getMaximum()) {
    if(prnt) cout<<"Mismatch in maximum "<<maximum<<" vs "
                  <<par->getMaximum()<<" for "<<Name<<endl;
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

string LKParameter::getTypeName()        {return getTypeName(type);}
string LKParameter::getTypeName(int type){
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
  cout<<formatI(id,3)
      <<" "<<leftJustify(getName(),20)
      <<" "<<getTypeName()
      <<" "<<positiveFlag(experiment);
}

void LKParameter::printInitial(){
  printHeader();
  cout<<formatG(initialValue,8,2)
      <<formatG(step,8,2)
      <<formatG(minimum,8,2)
      <<formatG(maximum,8,2)
      <<endl;
}

void LKParameter::printCurrent(bool withError){
  printHeader();
  cout<<" "<<formatRG(currentValue,8,2);
  if(withError && (type==PARAMETER_OF_INTEREST || type==NUISANCE_PARAMETER || type==FREE_PARAMETER)) {
    if(sigma>0.) cout<<" +- "<<formatLG(sigma,8,2);
    else         cout<<" +-     (undef)";
  }
  else           cout<<"                ";
  cout<<" "<<formatRG(getCurrentValueInMinuitUnits(),8,2)<<endl;
}

SigmaParameter::~SigmaParameter(){}
SigmaParameter::SigmaParameter() : LKParameter(PAR_SIGMA,PARAMETER_OF_INTEREST
                ,"Sigma",0.,0.001,-50.,50.) {setCommon();}

TSystBkgParameter::~TSystBkgParameter(){}
TSystBkgParameter::TSystBkgParameter(int run) : LKParameter(PAR_SYST_BKG_TVALUE
        ,NUISANCE_PARAMETER ,getTheName(run),0.,0.01,-5.,5) {}

string TSystBkgParameter::getTheName(int run){

  return "Syst. bkgd t-Value run"+formatI(run,2) ;
}


//----------TGauss--------///
TGaussParameter::~TGaussParameter(){}

// this is a t_value
TGaussParameter::TGaussParameter(int b,int run) : LKParameter(PAR_GAUSS_TVALUE+b
        ,NUISANCE_PARAMETER ,getTheName(b, run),0.,0.01,-5.,5.) { uncertainty = 1. ;}

string TGaussParameter::getTheName(int b, int run){
  return "Gauss NP "+formatI(b,3) + " run"+formatI(run,2) ;
}


//------------------------//

TStatBkgParameter::~TStatBkgParameter(){}

TStatBkgParameter::TStatBkgParameter() : LKParameter(){ stat_error = 0.;}

// this is not a t_value
TStatBkgParameter::TStatBkgParameter(int b,int run) : LKParameter(PAR_STAT_BKG_TVALUE+b
        ,NUISANCE_PARAMETER ,getTheName(b, run),0.1,0.01,0.,1.) {band=b; stat_error =0.;}

string TStatBkgParameter::getTheName(int b, int run){
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
Likelihood::Likelihood() : XeStat()           {setup();}
Likelihood::Likelihood(string nam) : XeStat() {
  setName(nam);
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
  cout<<" Id Name                 Type         Exp    Initial  Step    Min     Max"
      <<endl;
}

void Likelihood::printCurrentHeader(){
  cout<<" Id Name                 Type         Exp  Value               For Minuit"
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
  if(printLevel>1){
    cout<<"Adding parameter "<<param->getName()<<endl;
  }
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
    if(printLevel>0) {
      cout<<"Replacing parameter "<<id<<" with existing common one"<<endl;
    }
  }
  else {
    if(printLevel>0) cout<<"Keeping same parameter "<<id<<endl;
  }
}

void Likelihood::addParameter(int id, int type,string nam,double initialVal
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
  cout<<endl<<Name<<" has "<<parameters.size()<<" parameters:"<<endl;
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
  MethodCounter::count("Compute LogLikelihood");
  currentLikelihood->setCurrentValuesInMinuitUnits(values); // write the current values to LKParameters
  int printLevel=XeStat::getPrintLevel();
  if(printLevel>3) {
    cout<<"Evaluating likelihood "<<endl;
    currentLikelihood->printCurrentParameters();
  }
  double e= -1. * currentLikelihood->computeTheLogLikelihood();
  if(printLevel>3) {
    cout<<"             .... result:"<<XeCore::formatF(e,19,8)<<endl; 
  }
  return e;
}

double Likelihood::maximize(bool freezeParametersOfInterest){
  currentLikelihood=this;
  int np=mapMinuitParameters(freezeParametersOfInterest);

  nActiveParameters = getNActiveParameters();

  if(printLevel>1) {
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
    if(printLevel>1) {
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
//  min->SetPrintLevel(0);
  min->SetPrintLevel(max(0,printLevel-3));
  for(int i=0;i<np;i++){
    LKParameter* par=MinuitParameters[i];
    string pn=par->getName();
    double v=par->getInitialValueInMinuitUnits();
    double s = par->getStepInMinuitUnits();
    //double s= par->getStepInMinuitUnits(); // don't ask why but the factor 100 is needed for better convergence
    double vmi=par->getMinimumInMinuitUnits();
    double vma=par->getMaximumInMinuitUnits();
    if(printLevel>2) {
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
  MethodCounter::count("Calls to Minuit");
  min->Minimize(); 
  
  // write the post fit value to LKparameter 
  setCurrentValuesInMinuitUnits(min->X(),min->Errors());
  double ml= -1. * min->MinValue();
  if(printLevel>1) {
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
  if(printLevel>1) {
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
    if(printLevel>1) {
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




double Likelihood::shapeLikelihood(vector<double>* values
           ,vector<XeDist*>& dists,vector<double>& weight){
  int nd=dists.size();
  if(nd!=(int)weight.size()){
    cout<<"Inconsistency in parameters in shapeLikelihood"<<endl;
    return UNDEFINED;
  }
  vector<double> norm(nd);
  normalizeVector(weight,norm);
  int nv=values->size();
  double ll=0.;
  for(int iv=0;iv<nv;iv++){
    double v=(*values)[iv];
    double p=0;
    for(int d=0;d<nd;d++) p+= norm[d]*(dists[d]->pdf(v));
    if(p<=0.) return VERY_SMALL_LOG;
    ll += log(p);
  }
  return ll;
}

double Likelihood::shapeLikelihood(vector<double>* values, XeDist* dist){
  vector<double> one(1);
  vector<XeDist*> dists(1);
  one[0]=1.;
  dists[0]=dist;
  return shapeLikelihood(values,dists,one);
}

double Likelihood::shapeLikelihood(vector<int>* bins,int nDist, double** dists
              , double *norm){
  int nPoints=bins->size();

  double ll=0.;
  for(int p=0;p<nPoints;p++){
    int bin=(*bins)[p];
    double prob=0;
    for(int d=0;d<nDist;d++) prob+= norm[d] * (dists[d])[bin];
    if(prob<=0.) return VERY_SMALL_LOG;
    ll += log(prob);
  }
  return ll;
}


//----  Likelihood from a data set -----------------

string LikelihoodFromDataSet::getTheName(DataSet* data){
  string dn=data->getName();
  return "Maximum Likelihood on "+dn;
}

LikelihoodFromDataSet::~LikelihoodFromDataSet() {}

LikelihoodFromDataSet::LikelihoodFromDataSet(DataSet* da)
 : Likelihood(getTheName(da)) {
  data=da;
  values.resize(data->getNColumns());
}

double LikelihoodFromDataSet::computeTheLogLikelihood(){
  double l=0;
  int ne=data->getNEvents();
  for(int e=0;e<ne;e++) {
    data->getEntry(e,&values[0]);
    l+= computeIndividualLog(&values[0]);
  }
  return l;
}

DataSet* LikelihoodFromDataSet::getDataSet() {return data;}

//-------------- profile likelihood ----------------

ProfileLikelihood::~ProfileLikelihood(){}

ProfileLikelihood::ProfileLikelihood(string nam): Likelihood(nam), PValue(){
  setName(nam);
  initializeIt(PL_ANALYSIS);
  setup();
}

ProfileLikelihood::ProfileLikelihood(): Likelihood(), PValue(){
  initializeIt(PL_ANALYSIS);
  setup();
}

void ProfileLikelihood::printFlagsAndParameters(){
  cout<<endl<<getName()<<endl;
}

void ProfileLikelihood::setup(){
  sigPar             = NULL;
  lowerLimit         = 0.;
}


bool ProfileLikelihood::checkPValue(){
  update();
  if(!checkConsistency()) return false;
  sigPar=getParameter(PAR_SIGMA);
  if(sigPar==NULL) {
    cout<<"There should be a 'SIGMA' parameter in "<<Name<<endl;
    return false;
  }
  else if(!sigPar->isOfInterest()){
    cout<<"The 'SIGMA' parameter should be of interest in "<<Name<<endl;
    return false;
  }
  if(printLevel>-1) printFlagsAndParameters();
  int threshold=inCombinedMode()? 0:-1;
  if(printLevel>threshold) printInitialParameters();
  return true;
}

double ProfileLikelihood::getWimpMass(){
	return 0.;
}

void ProfileLikelihood::updateSigmaUnit(){
  if(printLevel>0){
     cout<<Name<<" has set sigmaUnit to : "<<sigmaUnit<<endl;
  }
  if(sigPar==NULL) cout<<"No cross section parameter in "<<Name<<endl;
  else {
    sigPar->setMinimum(-50.);   //TEST_A 
    sigPar->setMaximum(50.); 
    sigPar->setInitialValue(0.); 
    sigPar->setStep(0.1); 
    sigPar->setMinuitUnit(1.0);
//    sigPar->setMinuitUnit(sigmaUnit);
/*    sigPar->setMinimum(eventsLowerLimit()*sigmaUnit); 
    sigPar->setMaximum(eventsUpperLimit()*sigmaUnit); 
    sigPar->setInitialValue(2.*sigmaUnit); 
    sigPar->setStep(.2*sigmaUnit); 
    sigPar->setMinuitUnit(sigmaUnit);
*/	//Ale sigma change
    if(printLevel>1) {
      cout<<"Parameters after setting the scale of sigma"<<endl;
      printInitialParameters();
    }
  }
}



void ProfileLikelihood::estimateCrossSection() {
  if(printLevel>1){
     cout<<"Fitting "<<Name<<" once, for null hypothesis"<<endl;
  }
  resetParameters();
  LogD=maximize(false);
  LKParameter *parsig=getParameter(PAR_SIGMA);
  if(parsig==NULL) {
    cout<<"Technical bug!"<<getName()<<" does not have a 'sigma' param"<<endl; 
    return;
  }
  if(printLevel>1){
     cout<<"Estimated sigma="<<sigmaHat<<" log(d)="<<LogD<<endl;
  }
}

XeGraph* ProfileLikelihood::newGraphOfLogLikelihood(EventRange *er){
  // warning ! needs first to exclusion->prepareForLimit
  if(er==NULL) er=XeRange::getDefaultEventRange();
  XeGraph* gr=new XeGraph("LogLikeLihood"+getName());
  gr->setLineColor(kBlue);
  gr->setXYLabels(LABEL_EVT,LABEL_LOG_LIKELIHOOD);
  int n=er->getNPoints();
  gr->set(n);
  double reference=0.;
  for(int i=0;i<n;i++){
    resetParameters();
    double events=er->getValue(i);
    setParameterValue(PAR_SIGMA,events*sigmaUnit);
    double ll=maximize(true);
    if(i==0) reference=ll;
    gr->setPoint(i,events,reference-ll);
  }
  return gr;
}

/*
double ProfileLikelihood::getMaximum(){

  LKParameter* par=getParameter(PAR_SIGMA);
  double LL = maximize(false);
  double post_fit_sigma = par->getCurrentValue();

  return post_fit_sigma;
}

*/


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

TGraph* ProfileLikelihood::getLikelihoodScanOfParameter( int n_points, LKParameter * par){
  // warning ! needs first to exclusion->prepareForLimit
  int n = n_points;

  TGraph* gr= new TGraph(n);

  double reference=0.;
  LKParameter* sig=getParameter(PAR_SIGMA);

  sig->setCurrentValue(0.); // conditional to no signal!

  double min= par->getMinimum();
  double max = par->getMaximum();

  double step = (max - min) / ((double) n);

  double ll_Denominator = maximize(true); // conditional best fit

    
  for(int i=0;i<n;i++){
    resetParameters();

    double value = min + ((double)i)*step;

    par->setCurrentValue( value);

    double ll_Numerator = maximize(false); // conditional fit!!!
    printCurrentParameters();
    double test_stat_q = -1.* ll_Numerator;
//   double test_stat_q = -2. * (  ll_Numerator - ll_Denominator);
    gr->SetPoint(i,value , test_stat_q);

//   cout << " \t\t\t\t\t\t test val " << test_stat_q  << "  value  "  << value << endl;
  }


  gr->SetTitle("Log Likelihood Scan on "+TString(par->getName()));
  gr->GetXaxis()->SetTitle(TString(par->getName()));
  gr->GetYaxis()->SetTitle("-2Log(L(#sigma)/L(#hat{#sigma})");

  return gr;

}

TGraph* ProfileLikelihood::getGraphOfParameter( int n_points, int param_index){
  // warning ! needs first to exclusion->prepareForLimit
  int n = n_points;

  TGraph* gr= new TGraph(n);
  gr->GetXaxis()->SetTitle("#sigma #times 10^{-45} cm^{2}");
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

double ProfileLikelihood::qExclusion(double sigma){
  if(printLevel>1){
     cout<<"Calculating q for sigma="<<sigma<<endl;
  }
  if(!initialize()) {
    cout<<"Technical bug; see author"<<endl;
    return VERY_LARGE;
  }
  LKParameter *parsig=getParameter(PAR_SIGMA);
  if(parsig==NULL) {
    cout<<"Technical bug!"<<getName()<<" does not have a 'sigma' param"<<endl; 
    return VERY_LARGE;
  }

  if(sigmaHat>=sigma) {
    if(printLevel>1){
      cout<<"No point in calculating L for given sigma, return 0."<<endl;
    }
    return 0.;
  }
  resetParameters();
  setParameterValue(PAR_SIGMA,sigma);
  double LogN=maximize(true);
  double q=2*(LogD-LogN);
  if(printLevel>1){
     cout<<"For sigma="<<sigma<<", LogN="<<LogN<<", LogD="<<LogD
         <<" q="<<q<<endl;
  }
  return q;
}


double ProfileLikelihood::pValueS(double sigma){
  if(!initialize()){
    cout<<"Can't compute p-Value of ill formed "<<getName()<<endl;
    return UNDEFINED;
  }
  // Ale. not sure pValueS(0) = pValueB() as it is now.
  // where is the distinction between expected and observed pVal? 
  // CHECKME
  int    nc=getNParametersForChi2();
  double q=qExclusion(sigma);
  double p=Chi2Dist::above(q,nc)/2.; // Divided by two is necessary see paper asympt. approx. 
  if(printLevel>2){
    cout<<"q ="<<q<<", nParameters Of Interest="<<nParametersOfInterest
        <<", chi2 computed with "<<nc<<" d.o.f" 
        <<", pValue="<<p<<endl;
  }
  return p;
}

//--------- Combining several profile likelihoods --------------------

CombinedProfileLikelihood::CombinedProfileLikelihood(string n)
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

double CombinedProfileLikelihood::nSignalPerCm2()     {return sigToEvents;}

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
  if(printLevel>0) cout<<"Building up combined PL "<<Name<<endl;

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
          if(printLevel>0) {
             cout<<" adding common parameter "<<id<<" ("<<param->getName()
                 <<")"<<endl;
          }
          addParameter(param);
          param->setExperiment(ALL);
        }     
        else if(param->compares(parameters[id],true)) {
          if(printLevel>0) {
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
    cout<<"There should be a 'SIGMA' parameter in "<<Name<<endl;
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
        if(printLevel>0) {
           cout<<" adding specific parameter "<<id<<" ("<<param->getName()
               <<") for experiment "<<exp<<", new id="<<param->getId()<<endl;
        }
      } 
    }
  }
  if(printLevel>-1) printInitialParameters();
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

void CombinedProfileLikelihood::setWimpMass(double mass){
  if(printLevel>0) {
    cout<<"Setting mass for "<<Name<<": "<<mass<<endl;
  }
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    pl->setWimpMass(mass);
  }
}


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

bool CombinedProfileLikelihood::update(){
  sigToEvents=0;
  TRAVERSE_EXPERIMENTS(it) {
    ProfileLikelihood* pl=it->second;
    initializedOK= initializedOK && pl->update();
    sigToEvents += pl->nSignalPerCm2();
  }
  if(printLevel>0) {
    cout<<"Overall sigma to events for "<<Name<<": "<<sigToEvents<<endl;
  }
  return initializedOK;
}


//--------- Simple extension of ProfileLikelihood for signal + background

PLcountingSB::~PLcountingSB() {}
PLcountingSB::PLcountingSB(double v,double b, double db,double e, double de) : 
   SignalAndBackground(v,b,db,e,de) , ProfileLikelihood("PL Simple Count s+b") {
  addParameter(new SigmaParameter());
  //addParameter(new TSystBkgParameter());
  if(dBackground<=0.) ignoreParameter(PAR_SYST_BKG_TVALUE);
  addParameter(new TEfficiencyParameter());
  if(dEfficiency<=0.) ignoreParameter(PAR_EFFICIENCY_TVALUE);
}


void   PLcountingSB::printFlagsAndParameters()   {
  if(isExperimentSet()) cout<<"Experiment "<<experiment <<"; ";
  cout<<"Parameters of "<<Name<<":"<<endl;
  printSBParameters();
}

double PLcountingSB::computeTheLogLikelihood(){
  double tb=getParameter(PAR_SYST_BKG_TVALUE)->getCurrentValue();
  double bkg=background+tb*dBackground;
  double te=getParameter(PAR_EFFICIENCY_TVALUE)->getCurrentValue();
  double eff=efficiency+te*dEfficiency;
  double sig=getParameter(PAR_SIGMA)->getCurrentValue();
  double mu=exposure*(bkg+eff*sigToEvents*sig);
  double ll=PoissonDist::logPdf(nObserved,mu)-tb*tb/2-te*te/2;
  if(printLevel>3) {
    cout<<"Computing log likelihood for "<<getName()<<endl;
    printFlagsAndParameters();
    cout<<"Sigma               : "<< sig 
        <<" (expected events   : "<<nSignalEvents(sig)<<") "<<endl
        <<"Efficiency t-value  : "<< te  <<endl
        <<"Real efficiency     : "<< eff <<endl
        <<"Background t-value  : "<< tb  <<endl
        <<"Real background     : "<< bkg <<endl
        <<"Expected total evts : "<< mu  <<endl
        <<"Log Likelihood is   : "<< ll  <<endl;
  }
  return ll;
}

void   PLcountingSB::setWimpMass(double ) {}
double PLcountingSB::eventsUpperLimit()   {
  return Exclusion::typicalUpperLimit(nObserved);
}
double PLcountingSB::nSignalPerCm2()     {return nSignalEvents(1.);}
bool   PLcountingSB::update(){
  if(printLevel>0){
    cout<<"Sigma to events for "<<Name<<" :"<<sigToEvents*exposure<<endl;
  }
  return true;
}


//------------- Exclusion a la Yellin ------------------------------

YellinPValue::~YellinPValue() {}
YellinPValue::YellinPValue(DataSet* data, XeDist* sig, double sToE, int c)
  : PValue() {
  setName("Yellin limit");
  initializeIt(NO_ANALYSIS);
  dataSet=data;
  column=c;
  signal=sig;
  lowerLimit=2.;
  sigToEvents=sToE;
}

YellinPValue::YellinPValue(OneDimSimulatedDataSet* data, double sToE):PValue() {
  setName("Yellin limit for "+data->getName());
  initializeIt(NO_ANALYSIS);
  dataSet=data;
  column=0;
  lowerLimit=2.;
  signal=data->getSignal();
  sigToEvents=sToE;
}

bool YellinPValue::update(){
  maxGapForOne=0;
  XeValues values("Values of "+dataSet->getName(),dataSet,column);
  if(printLevel>2){
    cout<<"Yellin data set was reordered"<<endl;
    values.printIt(2);
  }
  map<double,int>::iterator it=values.begin();
  double t=it->first; 
  double c=signal->cdf(t);
  it++;
  for(; it!=values.end();it++){
    double u=it->first;
    double p=signal->cdf(u);
    maxGapForOne=max(maxGapForOne,p-c);
    c=p;
  } 
  if(printLevel>1){
    cout<<"Yellin data set was updated, maxGap="<<maxGapForOne<<endl;
  }
  return true; 
}

double   YellinPValue::eventsUpperLimit() {
  return Exclusion::typicalUpperLimit(dataSet->getNEvents());
}
double   YellinPValue::nSignalPerCm2()    {return sigToEvents;}
bool     YellinPValue::checkPValue()      {return true;}
int      YellinPValue::forceCLs()         {return FORCE_FALSE;}
DataSet* YellinPValue::getDataSet()       {return dataSet;}

double YellinPValue::pValueS(double sigma){
  double mu=sigToEvents*sigma;
  double x=maxGapForOne*mu;
  if(x<=0.) return 0.;
  double dm=mu/x;
  int m=(int)dm;
  if(m==dm) m--;
  double c0=0;
  double fac=1.;
  int k;
  if(printLevel>2){
      cout<<endl<<"Looping PValue for sigma="<<sigma<<", mu="<<mu<<", x="<<x
          <<", m="<<m<<endl;
  }
  for(k=1;k<=m;k++){
    fac*=k;
    double c=exp(-k*x)*pow(k*x-mu,k)/fac*(1+k/(mu-k*x));
    c0 -=c;
    if(printLevel>2) cout<<" ... k="<<k<<", c="<<c<<", c0="<<c0<<endl;
    if(abs(c)<1.e-9) break;
  }
  if(printLevel>1){
    cout<<"PValue for sigma="<<sigma<<", mu="<<mu<<", x="<<x<<", m="<<m
        <<", up to k="<<k<<", c0="<<c0<<endl;
  }
  if(c0<0.) c0=0.;
  else if(c0>1.) c0=1.;
  return c0;
}

void YellinPValue::printFlagsAndParameters(){
  cout<<"PValue a la Yellin from data set "<<dataSet->getName()<<endl
      <<"  Column          : "<<column<<endl
      <<"  nEvents         : "<<dataSet->getNEvents()<<endl
      <<"  Sigma to events : "<<sigToEvents<<endl
      <<"  Maximum gap     : "<<maxGapForOne<<endl;
}
