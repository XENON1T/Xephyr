#include "XeMath.h"
#include "XeCore.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

XeMath::~XeMath() {}
XeMath::XeMath() : XeCore() {}

string XeMath::tValueName(double t){return "t=" +formatF(t,4,2);}

bool XeMath::printVector(double *v, int n, string header){
  STATIC_CONST  int N_COLS=5
                 , WIDTH=72/N_COLS-1
                 ;
  if(header.size()>0) cout<<header<<endl;
  for(int i=0;i<n;i++) {
    if(i%N_COLS==0) cout<<setw(3)<<i<<": ";
    cout<<" "<<setw(WIDTH)<<v[i];
    if((i+1)%N_COLS==0) cout<<endl;
  }
  if(n%N_COLS!=0) cout<<endl;
  return n>0;
}

bool XeMath::printLogVector(double *v, int n, string header){
  STATIC_CONST  int N_COLS=5
                 , WIDTH=72/N_COLS-1
                 ;
  if(header.size()>0) cout<<header<<endl;
  for(int i=0;i<n;i++) {
    if(i%N_COLS==0) cout<<setw(3)<<i<<": ";
    cout<<" "<<setw(WIDTH)<<log10(v[i]);
    if((i+1)%N_COLS==0) cout<<endl;
  }
  if(n%N_COLS!=0) cout<<endl;
  return n>0;
}

bool XeMath::printVector(vector<double>&v,string header){
  return printVector(&(v[0]),v.size(),header);
}

bool XeMath::printVector(int *v, int n, string header){
  STATIC_CONST  int N_COLS=10
                 , WIDTH=72/N_COLS-1
                 ;
  if(header.size()>0) cout<<header<<endl;
  for(int i=0;i<n;i++) {
    if(i%N_COLS==0) cout<<setw(3)<<i<<": ";
    cout<<" "<<setw(WIDTH)<<v[i];
    if((i+1)%N_COLS==0) cout<<endl;
  }
  if(n%N_COLS!=0) cout<<endl;
  return n>0;
}

bool XeMath::printVector(vector<int>&v,string header){
  return printVector(&(v[0]),v.size(),header);
}



double  XeMath::vectorSum(double *v, int n){
  double total=0.;
  for(int i=0;i<n;i++) total += v[i];
  return total;
}

bool XeMath::equalizeVector(double *v, int n){
  if(n<2) return false;
  double total=vectorSum(v,n);
  for(int i=0;i<n;i++) v[i]=total/n;
  return true;
}

bool XeMath::normalizeVector(double *v, int n,double *w){
  if(n<1) return false;
  double total=vectorSum(v,n);
  if(total==0.) return false;
  if(w==NULL) w=v;
  for(int i=0;i<n;i++) w[i]=v[i]/total;
  return true;
}

bool XeMath::equalizeVector(vector<double>&v) {
  return equalizeVector(&(v[0]),v.size());
}

bool XeMath::normalizeVector(vector<double>&v){
  return normalizeVector(&(v[0]),v.size());
}

bool XeMath::normalizeVector(vector<double>&v, vector<double>&w) {
  return normalizeVector(&(v[0]),v.size(),&(w[0]));
}

double XeMath::vectorSum(vector<double>&v)    {
  return vectorSum(&(v[0]),v.size());
}


int XeMath::indexInterval(double v,double vMi,double vMa,int nIntervals,bool l){
  if(l){
    if(vMi<=0.) vMi=1.;
    return indexInterval(log10(v),log10(vMi),log10(vMa),nIntervals);
  }
  if(v<vMi || v>vMa ||nIntervals<1) return -1;
  int i=(int)((v-vMi)/(vMa-vMi)*nIntervals);
  if(i==nIntervals) i--;
  return i;
}

double XeMath::computePolynomial(double x, int n, double *coef){
  if(x==0. || n==0) return coef[0];
  double s=0;
  for(int i=n;i>=0;i--) s=s*x+coef[i];
  return s;
}

//------------------- Tolerances --------------------------------

XeTolerance::~XeTolerance(){}
XeTolerance::XeTolerance(double eps)  {epsilon=eps;}

bool XeTolerance::check(double d1,double d2,bool print){
  double d[2];
  d[0]=d1;
  d[1]=d2;
  return check(d,2,print);
}

bool XeTolerance::check(double d1,double d2,double d3,bool print){
  double d[3];
  d[0]=d1;
  d[1]=d2;
  d[2]=d3;
  return check(d,3,print);
}

bool XeTolerance::check(double d1,double d2,double d3,double d4,bool print){
  double d[4];
  d[0]=d1;
  d[1]=d2;
  d[2]=d3;
  d[3]=d4;
  return check(d,4,print);
}

bool XeTolerance::check(vector<double> &d, bool print){
  return check(&(d[0]),d.size(),print);
}

bool XeTolerance::check(double* d, int n,bool print){
  bool ok=true;
  for(int i=0;i<n;i++) {
    for(int j=0;j<i;j++) if(abs(d[i]-d[j])>epsilon) ok=false;
  }
  if(print && !ok) {
    cout<<"Mismatch ";
    for(int i=0;i<n;i++) {
      cout<<d[i];
      if(i<n-1) cout<<"/";
    }
    cout<<endl;
  }
  return ok;
}

//-------------- Interpolate ----------------------------------

XeInterpolation::XeInterpolation() {}
XeInterpolation::~XeInterpolation(){}


double  XeInterpolation::interpolate(int mode, double *x, double *y
                                    , int n, double z){
  switch(mode) {
    case LINEAR    : return linear(x,y,n,z);
    case PARABOLIC : return parabolic(x,y,n,z);
  }
  cout<<"Undefined mode for interpolation :"<<mode<<endl;
  return UNDEFINED;
}

double XeInterpolation::linear(double *x, double *y, int n, double z){
  if(z<x[0])   return y[0];
  if(z>x[n-1]) return y[n-1];
  int index=0;
  for(;index<n-2;index++) if(z<x[index+1]) break;
  return linear(&x[index],&y[index],z);
}

double XeInterpolation::linear(double *x, double *y, double z){
  return y[0]+ (y[1]-y[0])*(z-x[0])/(x[1]-x[0]);
}

double XeInterpolation::parabolic(double *x, double *y, int n, double z){
  if(z<x[0])   return y[0];
  if(z>x[n-1]) return y[n-1];
  int index=0;
  for(;index<n-3;index++) if(z<x[index+1]) break;
  return parabolic(&x[index],&y[index],z);
}

double XeInterpolation::parabolic(double *xorg, double *y, double z){

  double x[3];
  double xa=(xorg[0]+xorg[1]+xorg[2])/3.;
  for(int i=0;i<3;i++) x[i]=xorg[i]-xa;

  double a =
    -(x[0]*y[1] - x[1]*y[0] - x[0]*y[2] + x[2]*y[0] + x[1]*y[2] - x[2]*y[1])
     /((x[0] - x[1])*(x[0]*x[1] - x[0]*x[2] - x[1]*x[2] + x[2]*x[2]));

  double b =
      (x[0]*x[0]*y[1] - x[1]*x[1]*y[0] - x[0]*x[0]*y[2]
     + x[2]*x[2]*y[0] + x[1]*x[1]*y[2] - x[2]*x[2]*y[1])
     /((x[0] - x[1])*(x[0]*x[1] - x[0]*x[2] - x[1]*x[2] + x[2]*x[2]));

  double c =
     -(- y[2]*x[0]*x[0]*x[1] + y[1]*x[0]*x[0]*x[2] + y[2]*x[0]*x[1]*x[1]
     - y[1]*x[0]*x[2]*x[2] - y[0]*x[1]*x[1]*x[2] + y[0]*x[1]*x[2]*x[2])/
     ((x[0] - x[1])*(x[0]*x[1] - x[0]*x[2] - x[1]*x[2] + x[2]*x[2]));

   return a*(z-xa)*(z-xa)+b*(z-xa)+c;
}

double XeInterpolation::exponential(double fraction, double y0, double y1){
  if(y0<=0.|| y1<=0.) {
    cout<<"Can't perform a log interpolation beween "<<y0<<" and "<<y1<<endl;
    return 0.;
  }
  double l0=log(y0);
  double l1=log(y1);
  return exp(fraction*l1+(1.-fraction)*l0);
}

double XeInterpolation::linear(double fraction, double y0, double y1){
  return fraction*y1 + (1.-fraction)*y0;
}

double XeInterpolation::interpolate(int mode, double fraction
                                   , double y0, double y1){
  switch(mode) {
    case LINEAR      : return linear(fraction,y0,y1);
    case EXPONENTIAL : return exponential(fraction,y0,y1);
  }
  cout<<"Invalid interpolation mode:"<<mode<<endl;
  return UNDEFINED;
}

//---------------------------  General Interval class ------------------- 

XeBins::~XeBins(){}
XeBins::XeBins()            : XeObject(){}
XeBins::XeBins(string name) : XeObject(name){}

XeBins::XeBins(string name, vector<double>& e) : XeObject(name){
  edges=e;
  establishTheMap();
}

XeBins* XeBins::newLog10(){
  vector<double> l(nBins+1);
  for(int b=0;b<=nBins;b++) l[b]=log10(edges[b]);
  return new XeBins("Log "+Name,l);
}


bool XeBins::printIt(int level){
  cout<<getName()<<" has "<<nBins<<" bins"<<endl;
  if(level>0) {
    for(int b=0;b<nBins;b++){
      printf("%2d : from %.2f to %.2f\n",b,getLowerEdge(b),getUpperEdge(b));
    }
  }
 return nBins>0;
}

void XeBins::establishTheMap(){
  nBins=edges.size()-1;
  vMin=edges[0];
  vMax=edges[nBins];
  double lastV=vMin-1;
  bins.clear();
  for(int i=0;i<nBins;i++){
    double v=edges[i];
    if(v<=lastV) {
      cout<<"Error in "<<getName()
          <<": edges of the interval aren't in strictly ascending order"
          <<endl; 
      return;
    }
    lastV=v;
    bins[v]=i;
  }
}

void XeBins::extend(double vmi, double vma){
  if(vmi<vMin) edges[0]=vmi;
  if(vma>vMax) edges[nBins]=vma;
  establishTheMap();
}

bool XeBins::contains(double v) {return ok && v>=vMin && v<vMax;}
bool XeBins::isOk()             {return ok;}                      
int  XeBins::getNBins()         {return ok? nBins:0;}
int  XeBins::getBin(double v)   {
  if(!contains(v)) return -1;
  map<double,int>::iterator it=bins.lower_bound(v);
  if(it!=bins.begin()) it--;
  int b=it->second;
  if(b>=nBins || b<0) cout<<"Houston, we have a serious problem"<<endl;
  return b;
}

pair<int,double> XeBins::getBinAndFraction(double x){
  if(x<=vMin)     return pair<int,double>(0,0.);
  else if(x>=vMax) return pair<int,double>(edges[nBins-1],1.);
  int b=getBin(x);
  double f= (x-edges[b])/(edges[b+1]-edges[b]);
  double epsilon=1.e-6;
  if(f<epsilon) f=0.;
  else if(f>1.-epsilon) {f=0.; b++;}
  return pair<int,double>(b,f);
}

map<int,double> XeBins::getBinsAndWeights(double y1, double y2){
  double x1=max(y1,vMin);
  double x2=min(y2,vMax);
  double overall=(x2-x1)/(y2-y1);
  pair<int,double> p1=getBinAndFraction(x1);
  pair<int,double> p2=getBinAndFraction(x2);
  int b1=p1.first;
  int b2=p2.first;

  map<int,double> result;
  if(b1==b2) result[b1]= overall;
  else{
    double f1=p1.second;
    double f2=p2.second;
    overall=overall/(b1+f1-b2-f2);
    for(int b=b1;b<=b2;b++){
      double f=1.;
      if(b==b1) f=1.-f1;
      else if(b==b2) f=f2;
      result[b]=f*overall;
    }
  }

  return result;
}

double  XeBins::getCenter(int i)    {return .5*(edges[i]+edges[i+1]); }
double  XeBins::getLowerEdge(int i) {return edges[i];}
double  XeBins::getUpperEdge(int i) {return edges[i+1];}
double  XeBins::getLowerEdge()      {return edges[0];}
double  XeBins::getUpperEdge()      {return edges[nBins];}
double* XeBins::getLowerEdges()     {return &edges[0];}
double* XeBins::getUpperEdges()     {return &edges[1];}

EquidistantBins::EquidistantBins() : XeBins() {} 
EquidistantBins::EquidistantBins(string what, int n,double a, double b) 
               : XeBins(getEquidistantName(what,n,a,b)) {
  ok=a<b && n>0;
  if(!ok) return;
  vMin=a;
  vMax=b;
  step=(b-a)/n;
  edges.push_back(a);
  for(int i=1;i<=n;i++) edges.push_back(a+step*(i));
  establishTheMap();
}
EquidistantBins::~EquidistantBins(){}
int EquidistantBins::getBin(double value){
  if(contains(value)) return (int)((value-vMin)/step);
  return -1;
}
string EquidistantBins::getEquidistantName(string w,int n, double a, double b){
  char title[200];
  sprintf(title,"%d equidistant bins of %s from %f to %f",n,w.c_str(),a,b);
  return string(title);
}

EquiContentBins::~EquiContentBins(){}
EquiContentBins::EquiContentBins() : XeBins() {}

EquiContentBins::EquiContentBins(string what, int nB, vector<double>* data
           , double vMi, double vMa) : XeBins(getEquiContentName(what,nB)) {
  map<double,int> datamap;
  int n=data->size();
  for(int i=0;i<n;i++) {
    double d=(*data)[i];
    if(d>=vMi && d<=vMa) datamap[d]=datamap[d]+1;
  }
  ok= fillFromMap(&datamap,nB,n);
}

bool EquiContentBins::fillFromMap(map<double,int>* data,int nB,int m){
  int n=m;
  if(n==0){
    for(map<double,int>::iterator it=data->begin();it!=data->end();it++){
      n+=it->second;
    }
  }
  nBins=nB;
  if(n<2*nB){
    cout<<"Only "<<n<<" events for "<<nBins<<" bins ... giving up"<<endl;
    return false;
  }
  int current=0;
  int bin=0;
  int bound=0;
  int previousW=0;
  double v=0.;
  double previousV=0.;
  for(map<double,int>::iterator it=data->begin();it!=data->end();it++){
    v=it->first;
    int w=it->second;
    if(current>=bound){
      double e=bound==0? v : (v*w+previousV*previousW)/(w+previousW);
      edges.push_back(e);
      bin++;
      bound=(int)(.5+(double)(bin)*n/nBins);
    }
    previousW=w; 
    previousV=v; 
    current+=w;
  } 
  edges.push_back(v);
  int s=edges.size();
  if(s!=nBins+1){
    cout<<"Internal mismatch: data.size="<<n<<" returns edges.size()="<<s
        <<" for "<<nBins<<" requested bins"<<endl;
    return false;
  }
  establishTheMap(); 
  return true;
}

string EquiContentBins::getEquiContentName(string w,int n){
  char title[200];
  sprintf(title,"%d equal contents bins of %s",n,w.c_str());
  return string(title);
}

//------------------- General range classes--------------------------------

S1Range*     XeRange::defaultS1Range=NULL;
ErRange*     XeRange::defaultErRange=NULL;
YRange*      XeRange::defaultYRange=NULL;
MassRange*   XeRange::defaultMassRange=NULL;
VEscRange*   XeRange::defaultVEscRange=NULL;
SigmaRange*  XeRange::defaultSigmaRange=NULL;
EventRange*  XeRange::defaultEventRange=NULL;
tValueRange* XeRange::defaultTValueRange=NULL;

bool XeRange::operator==(const XeRange& other){
  if(mode==GENERAL) return xs==other.xs;

  return mode==other.mode
      && nBins==other.nBins
      && lowest==other.lowest
      && highest==other.highest
       ;
}


XeRange::XeRange(): XeObject() {}
XeRange::~XeRange() {}

XeRange::XeRange(string what,int n,double v0,double v1,int m): 
  XeObject(getTheName(what,n,v0,v1)){
  lowest=v0;
  highest=v1;
  nPoints=n;
  nBins=nPoints-1;
  mode=m;
  referenceIndex=UNDEFINED_INT;
}

bool XeRange::printIt(int level){
  if(level<1) return true;
  cout<<getName()<<" has "<<nPoints<<" points"<<endl;
  if(level<2) return true;
  printVector(xs);
  return nBins>0;
}

string XeRange::getTheName(string what,int n,double v0, double v1){
  char title[150];
  sprintf(title,"%s range %d points from %.1f to %.1f",what.c_str(),n,v0,v1);
  return string(title);
}

void XeRange::setReferenceIndex(double def) {
  referenceIndex=getClosestIndex(def);
}
int  XeRange::getClosestIndex(double val){
  int index=0;
  double m=VERY_LARGE;
  for(int i=1;i<nPoints;i++){
    double d=abs(val-getValue(i));
    if(d<m) {m=d; index=i;}
  }
  return index;
}

double          XeRange::getMin()              {return lowest;}
double          XeRange::getMax()              {return highest;}
double          XeRange::getStep ()            {return step;}
int             XeRange::getReferenceIndex()   {return referenceIndex;}
int             XeRange::getNPoints()          {return nPoints;}
int             XeRange::getNBins()            {return nBins;}
int             XeRange::getMode()             {return mode;}
bool            XeRange::isLinear()            {return mode==LINEAR;}
bool            XeRange::isGeneral()           {return mode==GENERAL;}
bool            XeRange::isLogarithmic()       {return mode==LOG;}
double*         XeRange::getValues()           {return &(xs[0]);}
double          XeRange::getValue(int i)       {return xs[i];}
double          XeRange::getNextValue(int i)   {return xs[i+1];}
vector<double>* XeRange::getTheValues()        {return &xs;}

double          XeRange::centerOfBin(int b)    {
  return .5*(getValue(b)+getNextValue(b));
}

pair<double,double> XeRange::getInterval(int i){
  if(i<0 || i>nPoints-2) {
    cout<<"Invalid interval :"<<i<<endl;
    return  pair<double,double>(UNDEFINED,UNDEFINED);
  }
  return pair<double,double>(getValue(i),getNextValue(i));
}

void XeRange::completeRange(){
  xs.resize(nPoints);
  for(int i=0;i<nPoints;i++) xs[i]=computeValue(i);
}

S1Range* XeRange::getDefaultS1Range() {
  if(defaultS1Range==NULL) defaultS1Range=new S1Range();
  return defaultS1Range;
}

ErRange* XeRange::getDefaultErRange() {
  if(defaultErRange==NULL) defaultErRange=new ErRange();
  return defaultErRange;
}

YRange* XeRange::getDefaultYRange() {
  if(defaultYRange==NULL) defaultYRange=new YRange();
  return defaultYRange;
}

MassRange* XeRange::getDefaultMassRange() {
  if(defaultMassRange==NULL) defaultMassRange=new MassRange();
  return defaultMassRange;
}

VEscRange* XeRange::getDefaultVEscRange() {
  if(defaultVEscRange==NULL) defaultVEscRange=new VEscRange();
  return defaultVEscRange;
}

SigmaRange* XeRange::getDefaultSigmaRange() {
  if(defaultSigmaRange==NULL) defaultSigmaRange=new SigmaRange();
  return defaultSigmaRange;
}

EventRange* XeRange::getDefaultEventRange() {
  if(defaultEventRange==NULL) defaultEventRange=new EventRange();
  return defaultEventRange;
}

tValueRange* XeRange::getDefaultTValueRange() {
  if(defaultTValueRange==NULL) defaultTValueRange=new tValueRange();
  return defaultTValueRange;
}

XeRange* XeRange::findDefaultRange(){
  if(*this==*getDefaultS1Range()) return getDefaultS1Range();
  if(*this==*getDefaultErRange()) return getDefaultErRange();
  if(*this==*getDefaultTValueRange()) return getDefaultTValueRange();
  if(*this==*getDefaultVEscRange()) return getDefaultVEscRange();
  if(*this==*getDefaultMassRange())  return getDefaultMassRange();
  if(*this==*getDefaultSigmaRange()) return getDefaultSigmaRange();
  if(*this==*getDefaultEventRange()) return getDefaultEventRange();
  return NULL;
}

GeneralRange::GeneralRange(string what,int n,double* v)
  : XeRange(what,n,v[0],v[n-1],GENERAL) {
  step=0.;
  for(int i=0;i<n;i++) {xs.push_back(v[i]);bins[i]=v[i];}
}
GeneralRange::GeneralRange() : XeRange(){}
GeneralRange::~GeneralRange(){}
double GeneralRange::computeValue(double i) {
  if(i<0) return xs[0];
  else if (i>=nPoints-1) return xs[nPoints-1];
  int index=(int)i;
  double f=i-index;
  return f*xs[index+1]+(1.-f)*xs[index];
}

int  GeneralRange::getIndex(double v)   {
  map<double,int>::iterator it=bins.lower_bound(v);
  if(it==bins.end()) return -1;
  if(it!=bins.begin()) it--;
  int b=it->second;
  if(b>=nPoints || b<0) cout<<"Houston, we have a serious problem"<<endl;
  return b;
}
double  GeneralRange::computeValue(int i)        {return xs[i];}
void    GeneralRange::setValue(int i, double v)  {xs[i]=v;}

LinearRange::LinearRange(string what,int n,double v0,double v1)
  : XeRange(what,n,v0,v1,LINEAR) {
  if(nPoints==1) step=0.;
  else           step=(highest-lowest)/(nPoints-1); 
  completeRange();
}
LinearRange::LinearRange(): XeRange(){}
LinearRange::~LinearRange(){}
double LinearRange::computeValue(int i)    {return lowest+i*step;}
double LinearRange::computeValue(double i) {return lowest+i*step;}
int    LinearRange::getIndex(double v)   {
  int b=(v-lowest)/step;
  if(b<0 || b>=nPoints) return -1;
  return b;
 }
pair<int,double> LinearRange::getIndexAndFraction(double x){
  double bin=(x-lowest)/step;
  if(bin<0 || bin>=nBins) return pair<int,double>(-1,0.);
  int b=(int)(bin+.001);
  return pair<int,double>(b,x-lowest-b*step);
}

pair<int,double> LinearRange::getTolerantIndexAndFraction(double x){
  double bin=(x-lowest)/step;
  if(bin<0) return pair<int,double>(0,0.);
  else if(bin>=nBins) return pair<int,double>(nBins-1,1.);
  int b=bin;
  return pair<int,double>(b,x-lowest-b*step);
}

LogRange::LogRange(string what,int n,double v0,double v1)
  : XeRange(what,n,v0,v1,LOG) {
  if(nPoints==1) step=0.;
  else           step=log(highest/lowest)/(nPoints-1); 
  completeRange();
}
LogRange::LogRange():XeRange(){}
LogRange::~LogRange(){}
double LogRange::computeValue(double i)    {return lowest*exp(i*step);}
double LogRange::computeValue(int i)       {return lowest*exp(i*step);}
int    LogRange::getIndex(double v)   {
  int b=log(v/lowest)/step;
  if(b<0 || b>=nPoints) return -1;
  return b;
}


double PRLMassRange::prlMasses[PRLMassRange::N_PRL_MASSES]= {
      6. , 7.   , 8.   , 9.   , 10.  , 12.  , 13.  , 14.   , 16.  , 18.  , 
     20. , 22.  , 24.  , 26.  , 28.  , 30.  , 35.  , 40.   , 45.  , 50.  ,
     53. , 55.  , 57.  , 60.  , 65.  , 70.  , 80.  , 90.   , 100. , 200. ,
    300. , 400. , 500. , 600. , 700. , 800. , 900. , 1000. };
PRLMassRange::~PRLMassRange(){}
PRLMassRange::PRLMassRange() : 
  GeneralRange("PRL Mass Range",N_PRL_MASSES,prlMasses) {}


MassRange::MassRange(int n, double m0, double m1) 
   : LogRange("Mass",n,m0,m1) {
  setReferenceIndex(DEFAULT_MASS);
}
MassRange::~MassRange(){}

SigmaRange::SigmaRange(int n, double s0, double s1):LogRange("Sigma",n,s0,s1){
  setReferenceIndex(DEFAULT_SIGMA);
}
SigmaRange::~SigmaRange(){}

S1Range::S1Range(int n, double p0, double p1):LinearRange("S1",n,p0,p1){}
S1Range::~S1Range(){}

VEscRange::VEscRange(int n,double v0,double v1):LinearRange("V_{esc}",n,v0,v1){
  setReferenceIndex(DEFAULT_VESC);
}
VEscRange::~VEscRange(){}

tValueRange::tValueRange(int n,double v0,double v1)
   : LinearRange("t-value",n,v0,v1){
  setReferenceIndex(DEFAULT_TVAL);
}
tValueRange::~tValueRange(){}

ErRange::ErRange(int n, double e0, double e1) : LinearRange("ERecoil",n,e0,e1){}
ErRange::~ErRange(){}

EventRange::EventRange(int n, double s0, double s1)
		  : LinearRange("Number of events",n,s0,s1){}
EventRange::~EventRange(){}

YRange::YRange(int n, double y0, double y1) : LinearRange("Y",n,y0,y1){}
YRange::~YRange(){}

//----------------------------- Integration interface --------------------

static       integrable*  fi;
static int   intMode;

double integrating(double x){return fi->getValue(intMode,x);}

integrable::integrable()    {}
integrable::~integrable()   {}

double integrable::integrate(int mode ,double x0, double x1) {
  if(x0>=x1) return 0.;
  fi=this;
  intMode=mode;
  ROOT::Math::Functor1D wf(&integrating);
  ROOT::Math::GaussLegendreIntegrator ig;
  ig.SetFunction(wf);
  ig.SetRelTolerance(0.0001);
  ig.SetNumberPoints(200);
  return ig.Integral(x0,x1);
}


//------------ Wrapper to solvable interface -------------------------

solvable::solvable()  {setLimits(0.,1.);}

solvable::solvable(double xmi,double xma)  {setLimits(xmi,xma);}

void solvable::setLimits(double xmi,double xma)  {
  setSolverLimits(xmi,xma);
  theFunction=new TF1("solvable",this,&solvable::Evaluate,0.,1.,0
                    ,"solvable","Evaluate");
}

solvable::~solvable()              {delete theFunction;}

double solvable::solve(double w)   {
  if(getDebugLevel()>1) {
    cout<<"===> Solving f(x)="<<w<<", xMin="<<xmin<<", xMax="<<xmax<<endl;
  }
  double raw=theFunction->GetX(w,0.001,1.001,1E-3,10,false);
  double x=raw*range+xmin;
  if(getDebugLevel()>1) {
    cout<<"===> Solution: raw="<<raw<<" x="<<x<<endl;
  }
  return x; 
}

double solvable::Evaluate(double *raw, double*) {
  double x=xmin+range*raw[0];
  if(getDebugLevel()>2) {
    cout<<endl<<endl<<"===> Calling solver for r="<<raw[0]<<", x="<<x<<endl;
  }
  double r=getValue(x);
  if(getDebugLevel()>2) {
    cout<<"===> Solver called for r="<<raw[0]<<", x="<<x
        <<" and returned "<<r<<endl;
  }
  return r;
}

void solvable::setSolverLimits(double xmi,double xma) {
  xmin=xmi;
  xmax=xma;
  range=xmax-xmin;
}
