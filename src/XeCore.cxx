#include "XeCore.h"
#include "TGaxis.h"
// -- A static class to count method invocation ----

XeTool::XeTool()  {}
XeTool::~XeTool() {}

map<string,int>  MethodCounter::counts;

MethodCounter::MethodCounter() : XeTool() {}
MethodCounter::~MethodCounter()         {}
void MethodCounter::count(string method){counts[method]++;}
void MethodCounter::reset()             {counts.clear();}
void MethodCounter::list()              {
 cout<<endl<<"Breakdown of method invocation"<<endl;
 map<string,int>::iterator it;
 for(it=counts.begin();it!=counts.end();it++){
//   cout<<setw(6)<<it->second<<"  "<<it->first<<endl;
 }
}

//-- a Simple wrapper to stopper

stopper::stopper(string n):XeTool(), TStopwatch() {name=n;}
stopper::~stopper(){}

void stopper::restart() {Start(true);}

void stopper::print() {
  double cpu=CpuTime();
  cout<<name<<" performed in "<<formatF(cpu,6,3)<<" CPU seconds"<<endl;
}

void stopper::lapse() {
  double cpu=CpuTime();
  cout<<name<<" @ "<<formatF(cpu,6,3)<<"s"<<endl;
}

// -------------- A purely static class for general handling ------------

int      XeCore::debugLevel      =     0 ;
int      XeCore::traceLevel      =     0 ;
bool     XeCore::withWarnings    = false ;
bool     XeCore::savedWarnings   = false ;
bool     XeCore::coreInitialized = false ;
stopper* XeCore::mainStopper     =  NULL ;

     XeCore::XeCore()                      {}
     XeCore::~XeCore()                     {}
void XeCore::setDebugLevel(int l)          {debugLevel=l;}
void XeCore::setTraceLevel(int l)          {traceLevel=l;}
void XeCore::showTheReferences(bool s)     {XeObject::showTheReferences(s);}
void XeCore::suppressWarnings(bool sup)    {withWarnings=!sup;}
void XeCore::showWarnings(bool show)       {withWarnings=show;}
void XeCore::resetPrecision()              {cout<<RESET_PRECISION;}
void XeCore::saveWarnings()                {savedWarnings=withWarnings;}
void XeCore::restoreWarnings()             {withWarnings=savedWarnings;}
int  XeCore::getDebugLevel()               {return debugLevel;}
int  XeCore::getTraceLevel()               {return traceLevel;}
bool XeCore::doWarn()                      {return withWarnings;}

bool XeCore::initializeCore(){
  if(coreInitialized) return true;
  coreInitialized=true;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,NULL);
  XeTopObject::newInstance(); 
  mainStopper=new stopper("");
  return true;
}

bool XeCore::fileExists(string path){
  FileStat_t stat;
  return gSystem->GetPathInfo(path.c_str(),stat)==0;
}

string XeCore::addASpace(string t) {return t==""? t: t+" ";}

string XeCore::format0I(int v,int w){
  char f[10];
  char t[10];
  sprintf(f," 0%dd",w);
  f[0]='%';
  sprintf(t,f,v);
  return string(t);
}

string XeCore::formatI(int v,int w, int trailer){
  char f[10];
  char t[50];
  sprintf(f," %dd",w);
  f[0]='%';
  sprintf(t,f,v);
  return string(t)+string(trailer,' ');
}

string XeCore::formatLI(int v,int w, int p){
  return leftJustify(formatI(v,w,p),w,true);
}

string XeCore::formatRI(int v,int w, int p){
  return rightJustify(formatI(v,w,p),w,true);
}

string XeCore::formatF(double v,int w, int p){
  char f[10];
  char t[50];
  sprintf(f," %d.%df",w,p);
  f[0]='%';
  sprintf(t,f,v);
  return string(t);
}

TString XeCore::formatR(double v,int w, int p){
// this doesn't work! WTF
  char t[10];
  cout << v << endl;
  int laMadonna = int(fabs(v) * pow(10,w));
  cout << "WARNING!!!:::: XeCore::formatR is not working properly. don't use it!" << endl; 
  cout << ((int)fabs(v) ) * ((int)pow(10,w)) << endl;
  if(v>=0) sprintf(t,"%d", laMadonna );
  else    sprintf(t,"m%d", abs(laMadonna) );

  cout << "lll " << t << endl;

  return TString(t);
}

string XeCore::formatLG(double v,int w, int p){
  return leftJustify(formatG(v,w,p),w,true);
}

string XeCore::formatRG(double v,int w, int p){
  return rightJustify(formatG(v,w,p),w,true);
}

string XeCore::formatLF(double v,int w, int p){
  return leftJustify(formatF(v,w,p),w,true);
}

string XeCore::formatRF(double v,int w, int p){
  return rightJustify(formatF(v,w,p),w,true);
}


string XeCore::formatG(double v,int w, int p){
  char f[10];
  char t[50];
  sprintf(f," %d.%dg",w,p);
  f[0]='%';
  sprintf(t,f,v);
  return string(t);
}

string XeCore::formatLatex(double v, int p){
  string g=formatG(v,0,p);
  if(g[0]=='1' && g[1]=='e') g=g.substr(1);
  size_t i=g.find('e');
  if(i==string::npos) return g;
  g=g.substr(0,i)+"\\,10^{"+g.substr(i+1)+"}";
  if(g[0]=='\\') g=g.substr(2);
  return g;
}

vector<double> XeCore::decodeVector(string text){
  vector<double> v;
  istringstream is(text);
  while (1) {
    double t=CRAZY;
    is>>t;
    if(t==CRAZY) break;
    v.push_back(t);
  }
  return v;
}

string XeCore::yesOrNo(bool b)  {return b? "yes":" no";}
string XeCore::doOrDont(bool b) {return b? "":"don't";}

string XeCore::justify(string s, int w, bool left, bool trm){
  if(trm) s=trim(s);
  int l=s.length();
  if(l>=w) return s;
  if(left) return s+string(w-l,' ');
  else return string(w-l,' ')+s;
}

string XeCore::leftJustify(string s,int w,bool trm) {
  return justify(s,w,true,trm);
}

string XeCore::rightJustify(string s, int w, bool trm){
  return justify(s,w,false,trm);
}

string XeCore::removeHeader(string original, string header){
  int lh=header.size();
  if(original.substr(0,lh)!=header) return original;
  return original.substr(lh,string::npos);
}

string XeCore::trim(string s){
  unsigned long p=s.find_first_not_of(' ');
  if(p==string::npos) return "";
  int q=s.find_last_not_of(' ');
  return s.substr(p,q-p+1);
}

string XeCore::replaceSlashes(string s){
  string o=s;
  int l=o.size();
  for(int i=0;i<l;i++) if(o[i]=='/') o[i]='_';
  return o;
}

string XeCore::replaceBlanks(string s){
  string o=s;
  int l=o.size();
  for(int i=0;i<l;i++) if(o[i]==' ') o[i]='_';
  return o;
}

string XeCore::replaceParentheses(string s){
  string o=s;
  int l=o.size();
  for(int i=0;i<l;i++) {
    if(o[i]=='(' || o[i]==')' || o[i]=='{' || o[i]=='}' ) o[i]='_' ;
  }
  return o;
}

string XeCore::makeItALaTeXName(string s){
  if(!isALaTeXName(s)) return s;
  int l=s.size();
  string r="\\mathrm{";
  for(int i=0;i<l;i++){
    bool space=false;
    if(s[i]==' '){
       space=true;
       if(i<l-1){
         if(isALaTeXChar(s[i+1])) space=false;
       }
       if(i>2){
         if(s[i-2]=='\\' && s[i-1]==',') space=false;
       }
    }
    if(space) r+="\\, ";
    else       r+=s[i];
  }
  r+="}";
  return r;
} 

string XeCore::makeItAFileName(string s){
  return replaceParentheses(replaceBlanks(replaceSlashes(s)));
}

string XeCore::positiveFlag(int flag){
  if(flag>=0) return formatI(flag,4);
  switch(flag) {
    case UNDEFINED_INT : return "Und.";
    case SAME          : return "Same";
    case NEXT          : return "Next";
    case AUTO          : return "Auto";
    case ALL           : return " All";
    case NONE          : return "None";
    case GENERAL       : return "Gen.";
    case LINEAR        : return "Lin.";
    case LOG           : return "Log.";
    case EXPONENTIAL   : return "Exp.";
    case LOGX_LINY     : return "LogL";
    case PARABOLIC     : return "Par.";
  }
  return "????";
}

bool XeCore::isALaTeXChar(char s){
  STATIC_CONST  int N_LATEX_CHAR=4;
  STATIC_CONST  char lc[N_LATEX_CHAR]={'\\','$','_','^'}; 
  for(int c=0;c<N_LATEX_CHAR;c++) if(s==lc[c]) return true;
  return false;
}

bool XeCore::isALaTeXName(string s){
  int l=s.size();
  for(int c=0;c<l;c++) if(isALaTeXChar(s[c])) return true;
  return false;
}

string XeCore::upperCase(string s){
 string str=s;
 for(unsigned i=0;i<str.length();i++) str[i] = toupper(str[i]);
 return str;
}

string XeCore::lowerCase(string s){
 string str=s;
 for(unsigned i=0;i<str.length();i++) str[i] = tolower(str[i]);
 return str;
}
void XeCore::displayFlag(string title, int tab, string what, bool def){
  string d= def ? "   (default)"
                : "   (not the default)";
  cout<<string(tab,' ')<<" - " <<leftJustify(title,25)<<":"
       <<rightJustify(what,6)<<d<<endl;
}

void XeCore::displayFlag(string title, int tab, bool value,bool def){
  displayFlag(title,tab,yesOrNo(value),value==def);
}

void XeCore::displayFlag(string title, int tab, int value, int def){
  char w[12];
  sprintf(w,"%d",value);
  displayFlag(title,tab,string(w),value==def);
}

void XeCore::displayFlag(string title, int tab, double value, double def){
  char w[12];
  sprintf(w,"%.2f",value);
  displayFlag(title,tab,string(w),value==def);
}

string XeCore::getEdgeName(int edge) {
  switch(edge) {
    case NO_BAND     : return "";
    case LOWER_EDGE  : return "lower";
    case BAND_CENTER : return "center";
    case UPPER_EDGE  : return "upper";
  }
  return UNDEFINED_STRING;
}

void XeCore::fillHist(TH1* h, vector<double> &v){
  int n=v.size();
  for(int i=0;i<n;i++) h->Fill(v[i]);
}

TH1F* XeCore::getCumulated(TH1F* h){
  TH1F* hist = (TH1F*) h->Clone();
  int nBins=h->GetNbinsX();
  double sum=0;
  double total=h->Integral()+h->GetBinContent(0)+h->GetBinContent(nBins+1);
  for(int b=0;b<=nBins+1;b++){
    sum += h->GetBinContent(b);
    hist->SetBinContent(b,sum/total);
  }
  return hist;
}

void XeCore::saveObjects(string fileName){
  top->writeToSimpleFile(fileName,"Xephyr objects");
}

void XeCore::restoreObjects(string fileName){
  top->readFromSimpleFile(fileName,"Xephyr objects");
}

//--------------------   XeGraphics ----------------------------
TCanvas*   XeGraphics::canvas               = NULL;
TPad*      XeGraphics::pad                  = NULL;
int        XeGraphics::font                 = DEFAULT_FONT;
int        XeGraphics::xCanvas              = 1;
int        XeGraphics::yCanvas              = 1;
int        XeGraphics::wCanvas              = 0;

XeGraphics::XeGraphics()   {}
XeGraphics::~XeGraphics()  {}

void XeGraphics::setFont(int f)     {font=f;}
int  XeGraphics::getFont()          {return font;}

void XeGraphics::updateGraphics() {
  if(gPad!=NULL) gPad->Update(); 
  if(canvas!=NULL) canvas->Update();
}

int  XeGraphics::isLog(int scale)                {return scale==LOG? 1:0;}
bool XeGraphics::isScaleError(int flag,bool war) {return !isScaleOK(flag,war);}
bool XeGraphics::isScaleOK(int flag, bool warn)  {
  if(flag==LINEAR || flag==LOG || flag==AUTO) return true;
  if(warn) cout<<"Invalid LINEAR/LOG flag: "<<flag<<endl;
  return false;
}

void XeGraphics::newFrame(string title, double xmi,double xma
                       ,double ymi,double yma,string tx,string ty
                       ,int xMode, int yMode){
  if(isScaleError(xMode)) return;
  if(isScaleError(yMode)) return;
  string lt=makeItALaTeXName(title);
  TH1I *hFrame = new TH1I(lt.c_str(),lt.c_str(),10,xmi,xma);
  hFrame->SetMinimum(ymi);
  hFrame->SetMaximum(yma);
  hFrame->SetDirectory(0);
  hFrame->SetLineColor(0);
  gPad->SetLogx(isLog(xMode));
  gPad->SetLogy(isLog(yMode));
  drawHist(hFrame,tx,ty);
}

void XeGraphics::newFrameZ(string title, double xmi,double xma
                        ,double ymi,double yma,double zmi, double zma 
                        ,int zMode, string tx,string ty){
  if(isScaleError(zMode)) return;
  string lt=makeItALaTeXName(title);
  TH2F *hFrame = new TH2F(lt.c_str(),lt.c_str(),10,xmi,xma,10,ymi,yma);
  hFrame->SetMinimum(zmi);
  hFrame->SetMaximum(zma);
  hFrame->Fill(xmi,ymi,zmi/2);
  hFrame->SetDirectory(0);
  hFrame->SetLineColor(0);
  drawHist(hFrame,tx,ty,"COLZ");
  gPad->SetLogz(isLog(zMode));
}

void XeGraphics::drawHist(TH2* h, string tx,string ty,int zMode){
  if(isScaleError(zMode,false)) return;
  drawHist(h,tx,ty,"COLZ");
  gPad->SetLogz(isLog(zMode));
}

void XeGraphics::drawHist(TH1* h,int plot, string tx,string ty,string opt){
  if(isScaleError(plot,false)) return;
  drawHist(h,tx,ty,opt);
  if(plot!=AUTO) gPad->SetLogy(isLog(plot));
}

void XeGraphics::drawHist(TH1* h, string tx,string ty,string opt){

  h->SetStats(0);

  if(tx.size()>0) h->GetXaxis()->SetTitle(tx.c_str());
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetTitleOffset(1.10);

  if(ty.size()>0) h->GetYaxis()->SetTitle(ty.c_str());
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetYaxis()->SetTickLength(0.02);

  if(font>0){
    h->GetXaxis()->SetLabelFont(font);
    h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetLabelFont(font);
    h->GetYaxis()->SetTitleFont(font);
  }

  gStyle->SetTitleFont(font,"");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  h->Draw(opt.c_str());

}

void XeGraphics::saveCanvas(string header){
   string fName=canvas->GetName();
   if(header!="") fName=header+"_";
   fName=makeItAFileName(fName);
   string dir=".";
   STATIC_CONST  int N_DIRS=2;
   const string dirs[N_DIRS]={"../pictures","pictures"};
   for(int d=0;d<N_DIRS;d++) if(fileExists(dirs[d])) dir=dirs[d];
   if(dir==".") {
     cout<<"No picture directory, save picture in current directory"<<endl;
   }
   fName=dir+"/"+fName+".png";
   canvas->SaveAs(fName.c_str(),"replace");
}

void XeGraphics::newCanvas(string name,int nx, int ny, bool square) {
  if(name.size()==0) name="Xenon100";
  canvas= new TCanvas(name.c_str(),name.c_str(),800,square?800:600);
  if(square) canvas->SetFixedAspectRatio();
  canvas->cd();
  canvas->SetFillStyle(4000);
  if(nx!=1 || ny!=1) subCanvas(nx,ny);
}

void XeGraphics::subCanvas(double xmi,double xma,double ymi,double yma){
  canvas->cd();
  pad=new TPad("","",xmi,ymi,xma,yma);
  pad->SetFillStyle(4000);
  pad->Draw();
  pad->cd();
}

void XeGraphics::nextSubCanvas()       {subCanvas(NEXT);}
void XeGraphics::subCanvas(int window) {subCanvas(SAME,SAME,window);}

void XeGraphics::subCanvas(int nx,int ny, int w){

  if(nx!=SAME) xCanvas=nx;
  if(ny!=SAME) yCanvas=ny;
  if(w==NEXT) wCanvas++;
  else if(w!=SAME) wCanvas=w;

  wCanvas=wCanvas%(xCanvas*yCanvas);
  int ix=wCanvas%xCanvas;
  int iy=yCanvas-1-(wCanvas-ix)/xCanvas;
  double wx=.9/xCanvas;
  double wy=.9/yCanvas;
  subCanvas(.05+ix*wx,.05+(ix+1)*wx,.05+iy*wy,.05+(iy+1)*wy);
}

TH1F*  XeGraphics::drawCumulated(TH1F* h, int color){

  TH1F* hist=getCumulated(h);

// from ROOTSYS/tutorials/hist/twoscales.C
  gPad->Update();
  float rightmax = 1.1*hist->GetMaximum();
  float scale = gPad->GetUymax()/rightmax;
  hist->SetLineColor(color);
  hist->Scale(scale);
  hist->Draw("same");
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
          gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(color);
  axis->SetLabelColor(color);
  axis->Draw();

  return hist;
}

//------------------------- Basic Object in XeCore -----------------------

bool          XeObject::showReference = false ;
  
XeObject::XeObject() : TObject()            {initializeObject();}
XeObject::XeObject(string n) : TObject()    {initializeObject(n);}
XeObject::~XeObject()                       {if(up!=NULL) up->detach(this);}

void   XeObject::setName(const char* n)     {setName(string(n));}
void   XeObject::setReference(string ref)   {Reference=ref;}
void   XeObject::setLegend(string leg)      {legend= leg=="" ? getName():leg;}
void   XeObject::showTheReferences(bool s)  {showReference=s;}
void   XeObject::linkUp(XeObject* u)        {up=u;}
void   XeObject::setFrameName(string n)     {FrameName=n;}
void   XeObject::setLaTeXName(string n)     {LaTeXName=n;}
void   XeObject::addName(string n)          {setName(Name+" "+n);}
void   XeObject::doNotMarkAsChanged()       {}
void   XeObject::markFromCheck(bool ok)     {if(ok) markOK(); else markError();}
void   XeObject::markOK()                   {status=CHECKED_OK;}
void   XeObject::markUnchecked()            {status=UNCHECKED;}
bool   XeObject::printIt(int)               {return true;}
bool   XeObject::print(int level)           {return printIt(level);}
bool   XeObject::isChanged()                {return changed;}
bool   XeObject::isOK()                     {return status!=IN_ERROR;}
bool   XeObject::isCheckedOK()              {return status==CHECKED_OK;}
bool   XeObject::isError()                  {return status==IN_ERROR;}
bool   XeObject::check(bool recompute)      {return checkIt(recompute);}
bool   XeObject::isReferenced()             {return Reference.size()>0;}
bool   XeObject::withReference()            {return showReference;}
bool   XeObject::hasName()                  {return Name!=UNDEFINED_STRING;}
bool   XeObject::hasLegend()                {return getLegend()!="";}
int    XeObject::getStatus()                {return status;}
string XeObject::getStatusName()            {return getStatusName(status);}
string XeObject::getMiniStatusName()        {return getMiniStatusName(status);}
string XeObject::getLegend()                {return legend;}
string XeObject::getName()                  {return Name;}
string XeObject::getNameSpace()             {return addASpace(Name);}
string XeObject::getLaTeXName()             {return LaTeXName;}
string XeObject::getFrameName()             {return FrameName;}
string XeObject::getReference()             {return "Ref: "+Reference;}

bool   XeObject::checkIt(bool) {
  if(status==IN_ERROR) return false;
  status=CHECKED_OK;
  return true;
}
 
void   XeObject::markError() {
  status=IN_ERROR;
  cout<<"--> "<<getName()<<" marked in error"<<endl;
}

void   XeObject::openCanvas(int nx,int ny, bool s) { 
  XeGraphics::newCanvas(getName(),nx,ny,s);
}

const char*        XeObject::getNameChar()  {return Name.c_str();}
XeObject*          XeObject::getUp()        {return up;}
XeObject*          XeObject::getDown(int g) {return downs[g];}
int                XeObject::getNDowns()    {return downs.size();}
vector<XeObject*>* XeObject::getDowns()     {return &downs;}

string XeObject::getStatusName(int s){
  switch(s) {
    case  UNCHECKED  : return "Unchk";
    case  CHECKED_OK : return "   OK";
    case  IN_ERROR   : return "Error";
  }
  return UNDEFINED_STRING;
}

string XeObject::getMiniStatusName(int s){
  switch(s) {
    case  UNCHECKED  : return " ";
    case  CHECKED_OK : return "+";
    case  IN_ERROR   : return "-";
  }
  return UNDEFINED_STRING;
}

void XeObject::traceFlags() {traceTheFlags();}
void XeObject::traceTheFlags() {}

void XeObject::trace(trackingPoint tp, string method){
  if(traceLevel<1) return;
  if(tp==ENTRY) {
    cout<<endl<<"Status at entry of "<<method<<" ";
    mainStopper->lapse();
  } 
  if(traceLevel>1) traceTheFlags();
  if(tp==EXIT) {
    cout<<"was the Status at exit of "<<method<<" ";
    mainStopper->lapse();
    cout<<endl;
  } 
}

void XeObject::setName(string n)  {
  Name=n;
  setFrameName(makeItAFileName(Name));
  setLaTeXName(makeItALaTeXName(Name));
}

void  XeObject::setXenon100Reference(string ref){
  setReference(XENON_100_REFERENCE+":"+ref);
}

void XeObject::printReference(int ntab){
  if(!(isReferenced() && withReference())) return;
  string header=ntab==0?"": string(ntab,' ');
  cout<<header<<"["<<getReference()<<"]"<<endl;
}

void XeObject::initializeObject(string n)  {
  initializeCore();
  up=NULL;
  setName(n);
  setReference("");
  experiment=NONE;
  changed=true; 
  markUnchecked();
  if(top!=NULL) top->attach(this);
}
 
void XeObject::markAsChanged(bool propagateUp) {
  changed=true;
  if(propagateUp && up!=NULL && up!=this) up->markAsChanged();
}

void XeObject::markAsRecomputed(){
  changed=false;
  int n=downs.size(); 
  for(int i=0;i<n;i++) {
    XeObject* d=downs[i];
    if(d!=NULL) d->markAsRecomputed();
  }
}
 
bool XeObject::checkAll() {
 bool ok=true;
 int n=downs.size(); 
  for(int i=0;i<n;i++) {
    XeObject* d=downs[i];
    if(d!=NULL)  ok = ok && d->checkAll();
  }
  if(ok) ok=checkIt();
  if(!ok) {
    cout<<" ===> Problem with "<<Name<<endl;
    return false;
  }
  if(traceLevel>0) cout<<getName()<<"::checkAll returned "<<ok<<endl;
  return ok;
}

void XeObject::attach(XeObject* d){
  if(d==NULL) {
    cout<<getName()<<" can not attach a null object"<<endl;
    return;
  }
  if(traceLevel>1) cout<<getName()<<" attaches "<<d->getName()<<endl;
  markAsChanged();
  XeObject* dup=d->getUp();
  if(dup==this) return;  
  if(this!=top){
    if(dup!=NULL && dup!=top) {
      cout<<"XeObject '"<<getName()<<"' can not attach '"<<d->getName()
          <<"' which is already attached to '"<<dup->getName()<<"'"<<endl;
    }
    top->detach(d); 
  }
  downs.push_back(d);
  d->linkUp(this);
}

void XeObject::detachItself(){
  if(up==NULL) {
    cout<<"XeObject '"<<getName()<<"' has no father to detach from!"<<endl;
  }
  else if(up!=top) up->detach(this);
}

void XeObject::detach(XeObject* d){
  if(d==NULL) return;
  markAsChanged();
  for(vector<XeObject*>::iterator it=downs.begin();it!=downs.end();it++){
    XeObject* obj=*it;
    if(obj==d) {
      downs.erase(it);
      if(this!=top) top->attach(d);
      break;
    }
  }
}

void XeObject::printTree(int depth, int offset){
  string header(5*offset,' ');
  cout<<header<<"-> "<<getMiniStatusName()<<" "<<getName()<<endl;
  if(depth==0) return;
  int n=downs.size();
  for(int i=0;i<n;i++) downs[i]->printTree(depth-1, offset+1);
}


string XeObject::resultsDirectory      = "../Results" ;
string XeObject::getResultsDirectory() {return resultsDirectory+"/";}
void   XeObject::setResultsDirectory(string dir){
  if(dir!="") resultsDirectory=dir;
  cout<<"Results directory set to "<<resultsDirectory<<endl;
}

void XeObject::writeToSimpleFile(string fName,string oName){
  if(fName=="") fName=getFrameName();
  if(oName=="") oName=fName;
  string full=XeObject::getResultsDirectory()+fName+".root";
  TFile tf(full.c_str(),"recreate");
  tf.cd();
  write(oName);
  tf.Close();
}
void XeObject::readFromSimpleFile(string fName,string oName){
  if(fName=="") {
    fName=Name;
    if(fName=="") {
      cout<<"Can't read an object with no name in automatic file name mode"
          <<endl;
      return;  
    }
  }
  if(oName=="") oName=fName;
  string full=XeObject::getResultsDirectory()+fName+".root";
  TFile tf(full.c_str(),"old");
  tf.cd();
  read(oName);
  tf.Close();
}
void XeObject::write(string nam) {
  if(nam=="") nam=getName(); 
  Write(nam.c_str());
}
void XeObject::read(string nam)  {
  if(nam=="") nam=getName(); 
  Read(nam.c_str());
}

void XeObject::setTValueLegend(double t) {
  setLegend("t-value = "+formatF(t,5,1));
}

void XeObject::setMassLegend(double mass) {
  setLegend("M ="+formatF(mass,6,1)+" GeV/c^{2}");
}

void XeObject::setBandLegend(int band) { setLegend("Band "+formatI(band,3)); }

int  XeObject::getExperiment()      {return experiment;}
bool XeObject::isExperimentSet()    {return experiment!=NONE;}

void XeObject::setExperiment(int e) {
  experiment=e;
  int n=downs.size(); 
  for(int d=0;d<n;d++) downs[d]->setExperiment(experiment);
}


XeTopObject* XeObject::getTop()       {return top;}

// the top object to attach evrything

void XeTopObject::newInstance()              {top=new XeTopObject();}
XeTopObject::~XeTopObject()                  {}
XeTopObject::XeTopObject() : XeObject("top") {};


// --------------- An extension to trivial map ---------------------

XeTable::XeTable()         : map<int,double>(), XeObject()  {}
XeTable::XeTable(string n) : map<int,double>(), XeObject(n) {
  first=VERY_LARGE_INT;
  last=VERY_SMALL_INT;
} 

XeTable::~XeTable() {}
void XeTable::add(int i, double value){
  if(find(i)==end()) {
    (*this)[i]=value;
    first=min(first,i);
    last=max(last,i);
  }
  else cout<<"In table "<<Name<<", entry "<<i<<" already exists"<<endl;
}
bool XeTable::printIt(int ){
  cout<<"Table "<<Name<<endl;
  for(XeTableIterator it=begin();it!=end();it++) {
    cout<<setw(4)<<it->first<<" : "<<formatG(it->second,10,5)<<endl;
  }
  return true;
}

int XeTable::getFirst() {return first;}
int XeTable::getLast()  {return last;}

// --------------- S1/S2 flattener ---------------------------------------


Flattener::~Flattener()                     {}
Flattener::Flattener()         : XeObject() {}
Flattener::Flattener(string n) : XeObject(n){}
void Flattener::unflatten(TGraph* orig, TGraph* added){
  int nOrig=orig->GetN();  
  int nAdded=added->GetN();
  added->Set(nOrig+nAdded);
  double* S1=orig->GetX();
  double* flat=orig->GetY();
  for(int i=0;i<nOrig;i++) {
    double S2=unflatten(S1[i],flat[i]);
    added->SetPoint(i+nOrig,S1[i],S2);
  }
}

//------------------- Fully static class to handle S1S2 display ------------


int        S1S2Display::displayMode          = UNKNOWN_DISPLAY;
int        S1S2Display::nBands               = DEFAULT_N_BANDS;
int        S1S2Display::nSlices              = DEFAULT_N_SLICES;
double     S1S2Display::xmin                 = UNDEFINED;
double     S1S2Display::xmax                 = UNDEFINED;
double     S1S2Display::ymin                 = UNDEFINED;
double     S1S2Display::ymax                 = UNDEFINED;
bool       S1S2Display::cutsAreHatched       = true;    
Flattener* S1S2Display::currentFlattener     = NULL;

           S1S2Display::S1S2Display()       {}
           S1S2Display::~S1S2Display()      {}

void       S1S2Display::setCurrentFlattener(Flattener* f) {currentFlattener=f;}

void       S1S2Display::setNBands(int nb)   {nBands=nb;}
void       S1S2Display::setNSlices(int ns)  {nSlices=ns;}
void       S1S2Display::setXmin(double xmi) {xmin=xmi;}
void       S1S2Display::setXmax(double xma) {xmax=xma;}
void       S1S2Display::setYmin(double ymi) {ymin=ymi;}
void       S1S2Display::setYmax(double yma) {ymax=yma;}
void       S1S2Display::hatchTheCuts(bool h){cutsAreHatched=h;}
double     S1S2Display::getXmin()           {return xmin;}
double     S1S2Display::getXmax()           {return xmax;}
double     S1S2Display::getYmin()           {return ymin;}
double     S1S2Display::getYmax()           {return ymax;}
bool       S1S2Display::areCutsHatched()    {return cutsAreHatched;}
bool       S1S2Display::isModeDefined()     {
  return displayMode>=0 && displayMode<N_DISPLAYS;
}
Flattener* S1S2Display::getCurrentFlattener(){return currentFlattener;}

double S1S2Display::getY(double x, double y)   {
  switch(getMode()){
    case BAND_VS_SLICE       :
    case S2_VS_S1            : return y;
    case S2_OVER_S1_VS_S1    : return y/x;
    case FLATTENED_S2_VS_S1  : return currentFlattener->flatten(x,y);
  }
  cout<<"S1S2Display::getY called in invalid mode: "<<getName()<<endl;
  return UNDEFINED;
}

double S1S2Display::undoY(double x, double y)   {
  switch(getMode()){
    case BAND_VS_SLICE       :
    case S2_VS_S1            : return y;
    case S2_OVER_S1_VS_S1    : return y*x;
    case FLATTENED_S2_VS_S1  : return currentFlattener->unflatten(x,y);
  }
  cout<<"S1S2Display::undoY called in invalid mode: "<<getName()<<endl;
  return UNDEFINED;
}

bool   S1S2Display::isFullyDefined() {
  return isModeDefined() && xmax!=UNDEFINED && xmin!=UNDEFINED
                         && ymax!=UNDEFINED && ymin!=UNDEFINED ;
}

void S1S2Display::setMode(int m){
  if(m==CURRENT_DISPLAY) {}
  else if(m<0 || m>=N_DISPLAYS){
    cout<<"Undefined mode "<<m<<endl;
    displayMode=UNKNOWN_DISPLAY;
    return;
  }
  else displayMode=m;
  if(getDebugLevel()>1) {
    cout<<"Display mode set to "<<getName()<<endl;
  }
}

string S1S2Display::getXLabel() {
  switch(getMode()) {
    case S2_VS_S1            :
    case FLATTENED_S2_VS_S1  :
    case S2_OVER_S1_VS_S1    : return LABEL_S1;
    case BAND_VS_SLICE       : return LABEL_SLICE;
  }
  return UNDEFINED_STRING;
}

string S1S2Display::getYLabel() {
  switch(getMode()) {
    case S2_VS_S1            : return LABEL_S2;
    case S2_OVER_S1_VS_S1    : return LABEL_S2_OVER_S1;
    case FLATTENED_S2_VS_S1  : return LABEL_FLATTENED;
    case BAND_VS_SLICE       : return LABEL_BAND;
  }
  return UNDEFINED_STRING;
}

string S1S2Display::getName() {
  if(getMode()>=0) return getYLabel()+" vs "+getXLabel();
  return "";
}

string S1S2Display::getCanvasName() {
  if(getMode()>=0)return makeItAFileName(getName());
  return "";
}

int S1S2Display::getMode(int m){
  if(m!=CURRENT_DISPLAY) setMode(m);
  if(isModeDefined()) return displayMode;
  cout<<"Invalid S1S2Display mode :"<<displayMode<<endl;
  return UNKNOWN_DISPLAY;
}

int S1S2Display::getDefaultLogMode(int m){
  int mod=getMode(m);
  switch(mod){
    case S2_VS_S1            : 
    case S2_OVER_S1_VS_S1    : return 1;
    case BAND_VS_SLICE    :
    case FLATTENED_S2_VS_S1  : return 0;
  }
  return UNDEFINED_INT;
}

bool S1S2Display::requiresComputation(int m){
  int mod=getMode(m);
  switch(mod){
    case S2_OVER_S1_VS_S1    : 
    case FLATTENED_S2_VS_S1  : return true;
  }
  return false;
}

bool S1S2Display::checkFull(){
  if(getMode()>=0) {
    if(isFullyDefined()) return true;
    cout<<"Missing X or Y limits for S1S2Display"<<endl;
  }
  return false;
}

void S1S2Display::setLimits(double xmi,double xma,double ymi,double yma){
  S1S2Display::setXmin(xmi);
  S1S2Display::setXmax(xma);
  S1S2Display::setYmin(ymi);
  S1S2Display::setYmax(yma);
  if(getDebugLevel()>1) {
    printf("Display limits set to %f<X<%f, %f<Y<%f \n",xmin,xmax,ymin,ymax);
  }
}

void S1S2Display::setTheDefaultLimits(){
  switch(getMode()){
    case S2_VS_S1            : setLimits(1.,31.,100.,20000.)  ; return;
    case BAND_VS_SLICE    : setLimits(0.,nSlices,0.,nBands); return;
    case S2_OVER_S1_VS_S1    : setLimits(1.,31.,10.,1000.)    ; return;
    case FLATTENED_S2_VS_S1  : setLimits(1.,31.,-1.5,1.5)     ; return;
  }
  cout<<"Don't know how to set the default limits"<<endl;
}

void S1S2Display::openCanvas(int nx, int ny,bool s) {
   newCanvas(getCanvasName(),nx,ny,s);
}

void S1S2Display::drawFrame(string title){
  if(!checkFull()) return;
  newFrame(title,xmin,xmax,ymin,ymax,getXLabel(),getYLabel());
}

void S1S2Display::drawFrameZ(string title,double zmin,double zmax,int zMode){
  if(!checkFull()) return;
  newFrameZ(title,xmin,xmax,ymin,ymax,zmin,zmax,zMode,getXLabel(),getYLabel());
}

void S1S2Display::drawBox(double x0,double x1,double y0,double y1
                        , int color, bool ratio){
  STATIC_CONST  int N_POINTS=13;
  double x[N_POINTS];
  double y[N_POINTS];
  int    nPoints=N_POINTS;

  if(getMode()==BAND_VS_SLICE){
    nPoints=5;
    x[0]=x0; x[1]=x0; x[2]=x1; x[3]=x1; x[4]=x0;
    y[0]=y0; y[1]=y1; y[2]=y1; y[3]=y0; y[4]=y0;
  }
  else {
    for(int i=1;i<7;i++){
      x[i]=x0+.2*(i-1)*(x1-x0);
      y[i]=getY(x[i],y1* (ratio? x[i]:1.));
      x[13-i]=x[i];
      y[13-i]=getY(x[13-i],y0* (ratio? x[13-i]:1.));
    }
    x[0]=x[12];
    y[0]=y[12];
  }  

  TPolyLine *pl=new TPolyLine(nPoints,x,y);
  if(color>=0) pl->SetFillColor(gStyle->GetColorPalette(color));
  pl->Draw("F");
}



//------- Stuff which can be drawn with frames and default limits in S1S2

S1S2Object::S1S2Object() : XeObject() {
  automaticLimits=false;
}

S1S2Object::S1S2Object(string name) :  XeObject(name) {
  setName(name);
  automaticLimits=false;
}

S1S2Object::~S1S2Object(){}

void S1S2Object::setAutomaticLimits(bool var) {automaticLimits=var;}
void S1S2Object::setS1S2Mode(int displayMode) {
  S1S2Display::setMode(displayMode);
}
int  S1S2Object::getS1S2Mode()              {return S1S2Display::getMode();}
bool S1S2Object::isDrawable(int displayMode){return displayMode!=BAND_VS_SLICE;}

void S1S2Object::drawS1S2(int displayMode,bool stayOnPrevious) {
  S1S2Display::setMode(displayMode);
  displayMode=S1S2Display::getMode();
  if(!isDrawable(displayMode)){
    cout<<"Can't draw "<<getName()<<" in displayMode "<<S1S2Display::getName()
        <<endl;
    return;
  }
  if(getDebugLevel()>1) {
    cout<<endl<<"Master draw of "<<Name<<endl;
  }
  if(gPad==NULL) S1S2Display::newCanvas();
  if(!stayOnPrevious) {
    if(displayMode!=BAND_VS_SLICE){
      drawTheFrame();
      gPad->SetLogy(S1S2Display::getDefaultLogMode());
    }
  } 
  draw();
}

void S1S2Object::drawTheFrame(){drawFrame();}

void S1S2Object::drawFrame(){
    if(automaticLimits) S1S2Display::setLimits(minX(),maxX(),minY(),maxY());
    else                S1S2Display::setTheDefaultLimits();
    gPad->SetLogy(S1S2Display::getDefaultLogMode());
    S1S2Display::drawFrame(Name);
}

void S1S2Object::drawFrameWithZ(double zmin, double zmax, int scaleZ){
    if(automaticLimits) S1S2Display::setLimits(minX(),maxX(),minY(),maxY());
    else                S1S2Display::setTheDefaultLimits();
    S1S2Display::drawFrameZ(Name,zmin,zmax,scaleZ);
}

//------------- Drawing style -------------------------------------

XeStyle::~XeStyle(){}
XeStyle::XeStyle(int lCol,int lWidth,int lStyle,int mCol,double mSize,int mStyle
                ,int fCol,int fStyle) : XeObject()  {
  setLineProperties(lCol,lWidth,lStyle);
  setMarkerProperties(mCol,mSize,mStyle);
  setFillProperties(fCol,fStyle);
}

void XeStyle::setLineProperties(int col, int width, int style) {
  setLineColor(col);
  setLineWidth(width);
  setLineStyle(style);
}

void XeStyle::setMarkerProperties(int col, double size, int style) {
  setMarkerColor(col);
  setMarkerSize(size);
  setMarkerStyle(style);
}

void XeStyle::setFillProperties(int col, int style) {
  setFillColor(col);
  setFillStyle(style);
}

int XeStyle::getStyleColor(int index){
  STATIC_CONST  int N_COLORS=10;
  STATIC_CONST  int COLORS[N_COLORS]
                  ={kBlack,kRed,kBlue,darkGreen,kSpring,9,30,38,46,49};
  int i=abs(index);
  if(i>N_COLORS) return DEFAULT_LINE_COLOR;
  return COLORS[i];
}

int XeStyle::getRainbowColor(int index, int nC, int mode){
  double frac= 1.- index/(double)(nC-1);
  if(mode==SYMMETRIC) frac= 1.- 2.*abs(frac-.5);
  int color=frac*49;
  return gStyle->GetColorPalette(color);
}

int      XeStyle::getLineStyle()             {return lineStyle;}
int      XeStyle::getLineColor()             {return lineColor;}
int      XeStyle::getLineWidth()             {return lineWidth;}
void     XeStyle::setLineStyle(int style)    {lineStyle=style;}
void     XeStyle::setLineWidth(int w)        {lineWidth=w;}
void     XeStyle::setLineColor(int col)      {lineColor=col;}
void     XeStyle::setLineStyleColor(int i)   {lineColor=getStyleColor(i);}

int      XeStyle::getMarkerStyle()           {return markerStyle;}
int      XeStyle::getMarkerColor()           {return markerColor;}
double   XeStyle::getMarkerSize()            {return markerSize;}
void     XeStyle::setMarkerStyle(int style)  {markerStyle=style;}
void     XeStyle::setMarkerSize(double s)    {markerSize=s;}
void     XeStyle::setMarkerColor(int col)    {markerColor=col;}
void     XeStyle::setMarkerStyleColor(int i) {markerColor=getStyleColor(i);}

int      XeStyle::getFillStyle()             {return fillStyle;}
int      XeStyle::getFillColor()             {return fillColor;}
void     XeStyle::setFillStyle(int style)    {fillStyle=style;}
void     XeStyle::setFillColor(int col)      {fillColor=col;}

XeStyle* XeStyle::forEdge(int e) {
  return new XeStyle(kBlue,2,e==BAND_CENTER?1:2);
}
XeStyle* XeStyle::forTValue(double t) {
  return forEdge(t==0.?BAND_CENTER:LOWER_EDGE);
} 
bool     XeStyle::printIt(int)  {
  cout
   <<"Line color "<<lineColor<<" style "<<lineStyle<<" width "<<lineWidth
   <<endl
   <<"Marker color "<<markerColor<<" style "<<markerStyle<<" size "<<markerSize
   <<endl
   <<"Fill color "<<fillColor<<" style "<<fillStyle
   <<endl;
  return true;
}

//------- An object with style --------------------------

XeStylized::~XeStylized(){}
XeStylized::XeStylized()  : XeGraphics(), XeObject()         {}
XeStylized::XeStylized(string name, string x,string y, string z) 
   :XeGraphics(), XeObject(name) {
  setXYLabels(x,y,z);
  defaultXScale=LINEAR;
  defaultYScale=LINEAR;
}

void     XeStylized::setXYLabels(string x,string y,string z) {
  setXLabel(x);
  setYLabel(y);
  setZLabel(z);
}

void     XeStylized::applyStyle()            {}
void     XeStylized::setXLabel(string xL)    {xLabel=xL;}
void     XeStylized::setYLabel(string yL)    {yLabel=yL;}
void     XeStylized::setZLabel(string zL)    {zLabel=zL;}
void     XeStylized::setDefaultXScale(int s) {if(isScaleOK(s)) defaultXScale=s;}
void     XeStylized::setDefaultYScale(int s) {if(isScaleOK(s)) defaultYScale=s;}
int      XeStylized::getLineColor()          {return style.getLineColor();}
int      XeStylized::getLineStyle()          {return style.getLineStyle();}
int      XeStylized::getLineWidth()          {return style.getLineWidth();}
int      XeStylized::getMarkerColor()        {return style.getMarkerColor();}
int      XeStylized::getMarkerStyle()        {return style.getMarkerStyle();}
int      XeStylized::getFillStyle()          {return style.getFillStyle();}
int      XeStylized::getFillColor()          {return style.getFillColor();}
int      XeStylized::getDefaultXScale()      {return defaultXScale;}
int      XeStylized::getDefaultYScale()      {return defaultYScale;}
double   XeStylized::getMarkerSize()         {return style.getMarkerSize();}
string   XeStylized::getXLabel()             {return xLabel;}
string   XeStylized::getYLabel()             {return yLabel;}
string   XeStylized::getZLabel()             {return zLabel;}
string   XeStylized::generatedXLabel()       {return "";}
string   XeStylized::generatedYLabel()       {return "";}
XeStyle* XeStylized::getStyle()              {return &style;} 

void XeStylized::setRainbowColor(int index, int nC, int mode){
  int color=XeStyle::getRainbowColor(index,nC,mode);
  setLineColor(color);
  setMarkerColor(color);
  setFillColor(color);
}

void XeStylized::drawIt(int plot, string opt, int flag) {
  if(isScaleError(plot,false)) return;
  drawWithFrame(opt,flag);
  if(plot!=AUTO) gPad->SetLogy(isLog(plot));
}

void XeStylized::drawWithFrame(string opt, int flag) {
  drawFrame();
  draw(opt,flag);
}

void XeStylized::drawWithCanvasAndFrame(string opt, int flag) {
  openCanvas();
  drawWithFrame(opt,flag);
}

void XeStylized::setMarkerColor(int c)      {
  style.setMarkerColor(c);
  applyStyle();
}

void XeStylized::setMarkerStyleColor(int index){
  style.setMarkerStyleColor(index);
  applyStyle();
}

void XeStylized::setMarkerStyle(int s)      {
  style.setMarkerStyle(s);
  applyStyle();
}

void XeStylized::setMarkerSize(double s)    {
  style.setMarkerSize(s);
  applyStyle();
}

void XeStylized::setLineColor(int col)      {
  style.setLineColor(col);
  applyStyle();
}

void XeStylized::setLineStyleColor(int index){
  style.setLineStyleColor(index);
  applyStyle();
}

void XeStylized::setLineStyle(int sty)      {
  style.setLineStyle(sty);
  applyStyle();
}

void XeStylized::setLineWidth(int wid)      {
  style.setLineWidth(wid);
  applyStyle();
}

void XeStylized::setFillColor(int col)      {
  style.setFillColor(col);
  applyStyle();
}

void XeStylized::setFillStyle(int sty)      {
  style.setFillStyle(sty);
  applyStyle();
}

void     XeStylized::setStyle(XeStyle* sty)     {
  if(sty!=NULL) style=*sty; 
  applyStyle();
}

void XeStylized::setLineProperties(int lCol, int lWidth, int lStyle) {
  style.setLineProperties(lCol,lWidth,lStyle);
}

void XeStylized::setMarkerProperties(int mCol, double mSize, int mStyle) {
  style.setMarkerProperties(mCol, mSize, mStyle);
  setMarkerStyle(mStyle);
}

void XeStylized::setFillProperties(int fCol, int fStyle) {
  style.setFillProperties( fCol, fStyle);
}

void XeStylized::printStyle(){
  cout<<getName()<<" will be drawn with the labels: "<<endl
      <<" X  :"<<xLabel<<endl
      <<" Y  :"<<yLabel<<endl;
  style.printIt();
}

void XeStylized::drawFrame( double xmi, double xma, double ymi, double yma
                          , int xScale, int yScale ){
  if(xScale==AUTO) xScale=defaultXScale;
  if(yScale==AUTO) yScale=defaultYScale;
  if(isScaleError(xScale)) return;
  if(isScaleError(yScale)) return;
  if(xmi==AUTOMATIC)  xmi=getMinX();
  if(xma==AUTOMATIC)  xma=getMaxX();
  if(ymi==AUTOMATIC)  ymi=yScale==LOG? getMinYNotZero() : getMinY();
  if(yma==AUTOMATIC)  yma=getMaxY();
  string tx=getXLabel();
  string ty=getYLabel();
  if(tx=="") tx=generatedXLabel();
  if(ty=="") ty=generatedYLabel();
  double fa=yma>0? 1.1:.9;
  double fi=ymi>0? .9:1.1;
  XeGraphics::newFrame(getName(),xmi,xma,ymi*fi,yma*fa,tx,ty,xScale,yScale);
}


// ------------  Extension to XeGraph ------------------------

XeGraph::~XeGraph(){}
XeGraph::XeGraph() :  XeStylized() , graph() {}
XeGraph::XeGraph(string name, string Xlab, string Ylab) 
    :  XeStylized(name,Xlab,Ylab), graph(){
  setTheNames(name);
}

XeGraph::XeGraph(string fn,string name, int ,string Xlab,string Ylab,XeStyle* st
        ,string leg) : XeStylized(name,Xlab,Ylab)
        ,graph((getResultsDirectory()+fn).c_str()){
  if(name=="") name=fn;
  setStyle(st);
  setTheNames(name,leg);
} 

XeGraph::XeGraph(string name, int n,  string Xlab, string Ylab, XeStyle* st
 , string leg)  : XeStylized(name,Xlab,Ylab),graph(n) {
  setStyle(st);
  setTheNames(name,leg);
} 

XeGraph::XeGraph(string name, int n, double *x,double *y,  string Xlab
   , string Ylab , XeStyle* st,string leg ,double norm) 
  :  XeStylized(name,Xlab,Ylab), graph(n,x,y){
  setStyle(st);
  setTheNames(name,leg);
  if(norm!=1.) for(int i=0;i<n;i++) setPoint(i,x[i],y[i]*norm);
}

XeGraph::XeGraph(string name, TGraph& gr,  string Xlab, string Ylab,XeStyle* st
  ,string leg) : XeStylized(name,Xlab,Ylab), graph(gr) {
  setStyle(st);
  setTheNames(name,leg);
}

TGraph* XeGraph::getTGraph()                         {return &graph;}
double* XeGraph::getX()                              {return graph.GetX();}
double* XeGraph::getY()                              {return graph.GetY();}
double  XeGraph::getX(int i)                         {return getX()[i];}
double  XeGraph::getY(int i)                         {return getY()[i];}
int     XeGraph::getN()                              {return graph.GetN();}
void    XeGraph::setPoint(int i, double x, double y) {graph.SetPoint(i,x,y);}
void    XeGraph::set(int i)                          {graph.Set(i);}

void XeGraph::setTheNames(string name,string leg){
  setName(name);
  setLegend(leg);
  graph.SetNameTitle(getNameChar(),getNameChar());
}

void XeGraph::addGraph(TGraph& gr){
  int old=getN();
  int n=gr.GetN();
  double* x=gr.GetX();
  double* y=gr.GetY();
  set(old+n);
  for(int i=0;i<n;i++) setPoint(old+i,x[i],y[i]);
}

void XeGraph::draw(string options, int)  {graph.Draw(options.c_str());}
void XeGraph::drawS1S2(string opt){
  XeGraph* gr=this;
  if(S1S2Display::requiresComputation()) {
    double* x=gr->getX();
    double* y=gr->getY();
    int n=getN();
    gr=new XeGraph(getName());
    gr->setStyle(getStyle());
    for(int i=0;i<n;i++) {
      double z=S1S2Display::getY(x[i],y[i]);
      gr->setPoint(i,x[i],z);
    }
  }
  if(getDebugLevel()>1){
    cout<<"Modified XeGraph for S1S2:";
    gr->printIt(10);
  }
  gr->getTGraph()->Draw(opt.c_str());
}

void XeGraph::drawInsideS1S2(string opt){
  int n=getN();
  vector<double> vx;
  vector<double> vy;
  double* x=getX();
  double* y=getY();
  for(int i=0;i<n;i++) {
    if(x[i]>=S1S2Display::getXmin() && x[i]<=S1S2Display::getXmax()) {
      double z=max(
             min(S1S2Display::getY(x[i],y[i]),1.001*S1S2Display::getYmax())
           ,.999*S1S2Display::getYmin());
      vx.push_back(x[i]);
      vy.push_back(z);
    }   
  }
  int s=vx.size();
  XeGraph* gr=new XeGraph(getName()+" inside s1s2",s,&(vx[0]),&(vy[0]));
  gr->setStyle(getStyle());
  if(getDebugLevel()>1){
    cout<<"Modified XeGraph inside S1S2:";
    gr->printIt(10);
  }
  gr->getTGraph()->Draw(opt.c_str());
}

bool XeGraph::isCompatible(XeGraph* reference){
  int n=getN();
  if(reference->getN()!=n) return false;
  double* x=getX();
  double* xr=reference->getX();
  for(int i=0;i<n;i++) if(x[i]!=xr[i]) return false;
  return true;
}

void XeGraph::applyStyle(){
  graph.SetLineColor(getLineColor());
  graph.SetLineWidth(getLineWidth());
  graph.SetLineStyle(getLineStyle());
  graph.SetMarkerColor(getMarkerColor());
  graph.SetMarkerSize(getMarkerSize());
  graph.SetMarkerStyle(getMarkerStyle());
  graph.SetFillColor(getFillColor());
  graph.SetFillStyle(getFillStyle());
}

double XeGraph::getMinX(){
  double m=VERY_LARGE;
  int n=getN();
  double *x=getX();
  for(int i=0;i<n;i++) m=min(m,x[i]);
  return m;
}

double XeGraph::getMinY(){
  double m=VERY_LARGE;
  int n=getN();
  double *y=getY();
  for(int i=0;i<n;i++) m=min(m,y[i]);
  return m;
}

double XeGraph::getMinYNotZero(){
  double m=VERY_LARGE;
  int n=getN();
  double *y=getY();
  for(int i=0;i<n;i++) if(y[i]>0.) m=min(m,y[i]);
  return m;
}

double XeGraph::getMaxX(){
  double m=VERY_SMALL;
  int n=getN();
  double *x=getX();
  for(int i=0;i<n;i++) m=max(m,x[i]);
  return m;
}

double XeGraph::getMaxY(){
  double m=VERY_SMALL;
  int n=getN();
  double *y=getY();
  for(int i=0;i<n;i++) m=max(m,y[i]);
  return m;
}

double XeGraph::getSumX(){
  double s=0.;
  int n=getN();
  double *x=getX();
  for(int i=0;i<n;i++) s+=x[i];
  return s;
}

double XeGraph::getSumY(){
  double s=0.;
  int n=getN();
  double *y=getY();
  for(int i=0;i<n;i++) s+=y[i];
  return s;
}

bool XeGraph::printIt(int level) {
  if(level<1) return true;
  printHeader();
  if(level<2) return true;
  int n=getN();
  double *x=getX();
  double *y=getY();
  for(int i=0;i<n;i++) {
   cout<<formatI(i,4)<<formatG(x[i],12,5)<<formatG(y[i],12,5)<<endl;
  } 
  return true;
}

void XeGraph::printHeader(){
  int n=getN();
  if(hasName()) cout<<Name<<" has ";
  cout<<n<<" points "<<getLegend()<<" "<<getXLabel()<<" "<<getYLabel()<<endl;
}
  
bool XeGraph::printGraph(int m,bool withStat,double xNorm,double yNorm) {
  printHeader();
  int n=getN();
  if(withStat){
    double ix=getMinX()*xNorm;
    double ax=getMaxX()*xNorm;
    double sx=getSumX()*xNorm;
    double iy=getMinY()*yNorm;
    double ay=getMaxY()*yNorm;
    double sy=getSumY()*yNorm;
    cout<<ix<<" <= x <= "<<ax<<" (sum="<<sx<<") "
        <<iy<<" <= y <= "<<ay<<" (sum="<<sy<<") "
        <<endl;
  }
  double *x=getX();
  double *y=getY();
  if(m>1 && m<n) {
    n=m;
    cout<<"First "<<n<<" elements:"<<endl;
  }
  for(int i=0;i<n;i++) {
     cout<<setw(15)<<setprecision(6)<<x[i]*xNorm
         <<setw(15)<<setprecision(6)<<y[i]*yNorm<<endl;
  }
  resetPrecision();
  return n>0;
}

XeGraph* XeGraph::divide(XeGraph* reference){
 if(!isCompatible(reference)){
   cout<<"Can't divide incompatible graphs"<<endl;
   return NULL;
 }
  int n=getN();
  double *x=getX();
  double *y=getY();
  double* yr=reference->getY();
  XeGraph* res=new XeGraph(getName()+" / "+reference->getName(),n
                  ,getXLabel(),getYLabel(),getStyle());
  for(int i=0;i<n;i++){
    double ny=yr[i]==0.? 0.:y[i]/yr[i];
    res->setPoint(i,x[i],ny);
  }
  res->setLegend(getLegend());
  return res;
}

XeGraphPositive::~XeGraphPositive() {}
XeGraphPositive::XeGraphPositive(string name, int n, double *x, double *y
 , string aX , string aY, XeStyle* st ,string le) : XeGraph(name, aX,aY) { 
   build(n,x,y,st,le); 
}

void XeGraphPositive::build(int n, double *x, double *y, XeStyle* st,string l){
  int nz=0;
  for(int i=0;i<n;i++) if(y[i]>0.) nz++;
  set(nz);
  nz=0;
  for(int i=0;i<n;i++) if(y[i]>0.) setPoint(nz++,x[i],y[i]);
  setStyle(st);
  setLegend(l);
}

XeGraphPositive::XeGraphPositive(string name, TGraph* gr, string aX, string aY
  , XeStyle* st ,string l) : XeGraph(name, aX,aY){
  build(gr->GetN(),gr->GetX(),gr->GetY(),st,l);
}

XeGraphPositive::XeGraphPositive(string name, XeGraph* gr) : XeGraph(name) {
  build(gr->getN(),gr->getX(),gr->getY(),gr->getStyle()
       ,gr->getLegend());
  setName(gr->getName());
}

//----------- Group of XeGraphs -----------------------------------

XeMultiGraph::~XeMultiGraph(){}

XeMultiGraph::XeMultiGraph() : XeStylized(), vector<XeGraph*>() {}

XeMultiGraph::XeMultiGraph(string name, string xA, string yA, string zA) 
   : XeStylized(name,xA,yA,zA),  vector<XeGraph*>() {
  setName(name);
}

XeMultiGraph::XeMultiGraph(string name,int ng,XeGraph** g, double* zs, string xA
    ,string yA ,string zA) : XeStylized(name, xA,yA,zA),  vector<XeGraph*>() { 
  setName(name);
  for(int i=0;i<ng;i++) {
    double z=zs==NULL? i : zs[i];
    add(g[i],z);
  }
}

XeMultiGraph::XeMultiGraph(string name,vector<XeGraph*>& g,vector<double>* vz
 ,string xA,string yA,string zA) : XeStylized(name,xA,yA,zA),vector<XeGraph*>(){
  setName(name);
  int ng=g.size();
  for(int i=0;i<ng;i++) {
    double z=vz==NULL? i : (*vz)[i];
    add(g[i],z);
  } 
}

void     XeMultiGraph::setReferenceGraph(XeGraph* g){referenceGraph=g;}
void     XeMultiGraph::setReferenceGraph(int i)     {referenceGraph=(*this)[i];}
XeGraph* XeMultiGraph::getReferenceGraph()          {return referenceGraph;}

double          XeMultiGraph::getX(int g)           {return Vx[g];}
double*         XeMultiGraph::getX()                {return &Vx[0];}
vector<double>* XeMultiGraph::getVx()               {return &Vx;}
double          XeMultiGraph::getZ(int g)           {return Vz[g];}
double*         XeMultiGraph::getZ()                {return &Vz[0];}
vector<double>* XeMultiGraph::getVz()               {return &Vz;}

void XeMultiGraph::applyStyle() {
  int ng=size();
  for(int i=0;i<ng;i++) (*this)[i]->setStyle(&style);
}

XeMultiGraph* XeMultiGraph::normalize(XeGraph* refer){
  XeMultiGraph* mg=new XeMultiGraph(*this);
  int ng=size();
  for(int i=0;i<ng;i++) (*this)[i]= (*this)[i]->divide(refer);
  return mg;
}

XeMultiGraph* XeMultiGraph::normalize(int i) {
  if(i==AUTOMATIC) return normalize(referenceGraph);
  return normalize((*this)[i]);
}

void XeMultiGraph::modifyLineColor(int color) {
  int ng=size();
  for(int i=0;i<ng;i++) (*this)[i]->setLineColor(color);
}

void XeMultiGraph::modifyLineStyle(int s) {
  int ng=size();
  for(int i=0;i<ng;i++) (*this)[i]->setLineStyle(s);
}

void XeMultiGraph::modifyLineWidth(int width) {
  int ng=size();
  for(int i=0;i<ng;i++) (*this)[i]->setLineWidth(width);
}

void XeMultiGraph::setRainbowColors(int mode){
  int ng=size();
  for(int i=0;i<ng;i++) (*this)[i]->setRainbowColor(i,ng,mode);
}

void XeMultiGraph::setSymmetricRainbowColors()  {setRainbowColors(SYMMETRIC);}
void XeMultiGraph::setAsymmetricRainbowColors() {setRainbowColors(ASYMMETRIC);}

void XeMultiGraph::add(XeGraph* g,int color, string leg)  {
  g->setLineColor(color);
  if(leg!="") g->setLegend(leg);
  add(g);
}

void XeMultiGraph::add(XeGraph* g,double z)  {
  
  int n=g->getN();  
  double* x=g->getX();
  if(size()==0) {
    Vx.resize(n);
    for(int i=0;i<n;i++) Vx[i]=x[i];
  }
  else {
    int nX=Vx.size();
    if(nX==n){
      int i;
      for(i=0;i<n;i++) if(Vx[i]!=x[i]) break;
      if(i!=n){
        cout<<"Warning adding graph "<<g->getName()
            <<" with incompatible x-values"<<endl;
      }
    }
    else {
      cout<<"Warning adding graph "<<g->getName()<<" with incompatible length"
          <<endl;
    }
  }
  defaultXScale=g->getDefaultXScale();
  defaultYScale=g->getDefaultYScale();
  push_back(g);
  if(z==UNDEFINED) z=size()-1;
  Vz.push_back(z);
}

XeGraph* XeMultiGraph::getGraph(int i)  {
  if(i>=0 && i<(int)size()) {
    XeGraph *gr=(*this)[i];
    if(gr==NULL) {
      cout<<"XeGraph  #"<<i<<" in "<<getName()<<" isn't defined"<<endl;
    }
    return gr;
  }
  cout<<"Can't access XeGraph #"<<i<<" in "<<getName()<<endl;
  return NULL;
}

void XeMultiGraph::add(XeMultiGraph* mg) {
  int n=mg->size();
  for(int i=0;i<n;i++) add(mg->getGraph(i),mg->getZ(i));
}

XeMultiGraph* XeMultiGraph::transpose(){
  int nx=Vx.size();
  int nz=Vz.size();
  XeMultiGraph *mg=new XeMultiGraph(getName()+" transposed",getZLabel()
                                   ,getYLabel(),getXLabel());
  for(int g=0;g<nx;g++){
    XeGraph *xg=new XeGraph(getName(),nz);
    xg->setLegend(formatG(Vx[g],8,4));
    for(int p=0;p<nz;p++) {
      XeGraph* go=(*this)[p];
      xg->setPoint(p,Vz[p],go->getY(g));
    }
    mg->add(xg,Vx[g]);
  }
  return mg;
}


void XeMultiGraph::forceAxes(int which){
  XeGraph *gr=getGraph(which); 
  if(gr==NULL) return;
  setXYLabels(gr->getXLabel(),gr->getYLabel());
}

void XeMultiGraph::resetTheLegends() {theLegends.clear();}
bool XeMultiGraph::addLegend(string leg){
  if(theLegends.find(leg)==theLegends.end()){
    theLegends.insert(leg);
    return true;
  }
  return false;
}

string XeMultiGraph::generatedXLabel(){
  string res="";
  for(unsigned i=0;i<size();i++){
    XeGraph *gr=getGraph(i);
    if(gr==NULL) continue;
    string axis=gr->getXLabel();
    if(res=="") res=axis;
    else if(res!=axis){
      cout<<"Incompatibility in X axis names for XeMultigraph "<<getName()
          <<endl<<"consider using XeMultiGraph::forceAxes()"<<endl;
    }
  }
  return res;
}

string XeMultiGraph::generatedYLabel(){
  string res="";
  for(unsigned i=0;i<size();i++){
    XeGraph *gr=getGraph(i);
    if(gr==NULL) continue;
    string axis=gr->getYLabel();
    if(res=="") res=axis;
    else if(res!=axis){
      cout<<"Incompatibility in Y axis names for XeMultigraph "<<getName()
          <<endl<<"consider using XeMultiGraph::forceAxes()"<<endl;
    }
  }
  return res;
}

void XeMultiGraph::draw(string options,int flag){
  TLegend *leg=NULL;
  int nLegend=0;
  int ng=size();
  bool withLegend=flag==WITH_LEGEND;
  if(withLegend) {
    resetTheLegends();
    for(int i=0;i<ng;i++) {
      if((*this)[i]->hasLegend()){
        string l=(*this)[i]->getLegend();
        if(addLegend(l)) nLegend++;
      }
    }
    leg=new TLegend(.65,.88-0.04*nLegend,.88,.88);
    leg->SetFillStyle(0);
    leg->SetTextFont(getFont());
    leg->SetTextSize(0.025);
  } 
  resetTheLegends();
  for(int i=0;i<ng;i++) {
    (*this)[i]->draw(options);
    if(withLegend && (*this)[i]->hasLegend()){
      string l=(*this)[i]->getLegend();
      if(addLegend(l)) leg->AddEntry((*this)[i]->getTGraph(),l.c_str(),"lp");
    }
  }
  if(nLegend>0) leg->Draw();
}

double XeMultiGraph::getMaxX(){
  double m=VERY_SMALL;
  int ng=size();
  for(int i=0;i<ng;i++) m=max((*this)[i]->getMaxX(),m);
  return m;
}

double XeMultiGraph::getMaxY(){
  double m=VERY_SMALL;
  int ng=size();
  for(int i=0;i<ng;i++) m=max((*this)[i]->getMaxY(),m);
  return m;
}

double XeMultiGraph::getMinX(){
  double m=VERY_LARGE;
  int ng=size();
  for(int i=0;i<ng;i++) m=min((*this)[i]->getMinX(),m);
  return m;
}

double XeMultiGraph::getMinY(){
  double m=VERY_LARGE;
  int ng=size();
  for(int i=0;i<ng;i++) m=min((*this)[i]->getMinY(),m);
  return m;
}

double XeMultiGraph::getMinYNotZero(){
  double m=VERY_LARGE;
  int ng=size();
  for(int i=0;i<ng;i++) m=min((*this)[i]->getMinYNotZero(),m);
  return m;
}

bool XeMultiGraph::printIt(int level){
  int ng=size();
  bool ok=true;
  for(int i=0;i<ng;i++) ok=ok && (*this)[i]->printGraph(0,level>1);
  return ok;
}


errorHandler::errorHandler(TString name) : className(name) { }

void errorHandler::Error(TString functionName, TString message){

	cout << className + "::" + functionName + " - ERROR : " + message << endl;
	//this is a fatal error need to quit or will produce unwanted results
	exit(100); 
}
void errorHandler::Warning(TString functionName, TString message){
	cout << className + "::" + functionName + " - WARNING : " + message << endl;
}
void errorHandler::Info(TString functionName, TString message){
	cout << className + "::" + functionName + " - INFO : " + message << endl;
}

