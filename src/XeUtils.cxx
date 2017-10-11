#include "XeUtils.h"

//! global print level by default set to all.
//! 0 =  Debug-print all, 1 = Info + Warning + Error, 2 = Warning +  Error, 3 = Only error
int errorHandler::globalPrintLevel = 1;

errorHandler::errorHandler(TString name) : className(name) { 

    localPrintLevel = -1;
}

void errorHandler::Error(TString functionName, TString message){

	cout << className + "::" + functionName + " - ERROR : " + message << endl;
	//this is a fatal error need to quit or will produce unwanted results
	std::runtime_error(message.Data()); 
}
void errorHandler::Warning(TString functionName, TString message){
    if(getPrintLevel() < 3 ) 
        cout << className + "::" + functionName + " - WARNING : " + message << endl;
}
void errorHandler::Info(TString functionName, TString message){
    if(getPrintLevel() < 2 )
	    cout << className + "::" + functionName + " - INFO : " + message << endl;
}

void errorHandler::Debug(TString functionName, TString message){
    if(getPrintLevel() < 1 )
	    cout << className + "::" + functionName + " - DEBUG : " + message << endl;
}




/// ------------    Print Toools -------------- ///
printTools::printTools(){};

TString printTools::upperCase(TString s){
 TString str=s;
 str.ToUpper();
 return str;
}

TString printTools::lowerCase(TString s){
 TString str=s;
str.ToLower();
 return str;
}

TString printTools::format0I(int v,int w){
  char f[10];
  char t[10];
  sprintf(f," 0%dd",w);
  f[0]='%';
  sprintf(t,f,v);
  return TString(t);
}

TString printTools::formatI(int v,int w, int trailer){
  char f[10];
  char t[50];
  sprintf(f," %dd",w);
  f[0]='%';
  sprintf(t,f,v);
  return TString(t)+TString(trailer,' ');
}

TString printTools::formatLI(int v,int w, int p){
  return leftJustify(formatI(v,w,p),w,true);
}

TString printTools::formatRI(int v,int w, int p){
  return rightJustify(formatI(v,w,p),w,true);
}

TString printTools::formatF(double v,int w, int p){
  char f[10];
  char t[50];
  sprintf(f," %d.%df",w,p);
  f[0]='%';
  sprintf(t,f,v);
  return TString(t);
}

TString printTools::formatR(double v,int w, int p){
// this doesn't work! WTF
  char t[10];
  cout << v << endl;
  int laMadonna = int(fabs(v) * pow(10,w));
  cout << "WARNING!!!:::: printTools::formatR is not working properly. don't use it!" << endl; 
  cout << ((int)fabs(v) ) * ((int)pow(10,w)) << endl;
  if(v>=0) sprintf(t,"%d", laMadonna );
  else    sprintf(t,"m%d", abs(laMadonna) );

  cout << "lll " << t << endl;

  return TString(t);
}

TString printTools::formatLG(double v,int w, int p){
  return leftJustify(formatG(v,w,p),w,true);
}

TString printTools::formatRG(double v,int w, int p){
  return rightJustify(formatG(v,w,p),w,true);
}

TString printTools::formatLF(double v,int w, int p){
  return leftJustify(formatF(v,w,p),w,true);
}

TString printTools::formatRF(double v,int w, int p){
  return rightJustify(formatF(v,w,p),w,true);
}


TString printTools::formatG(double v,int w, int p){
  char f[10];
  char t[50];
  sprintf(f," %d.%dg",w,p);
  f[0]='%';
  sprintf(t,f,v);
  return TString(t);
}


TString printTools::trim(TString str){
  string s = str.Data();
  unsigned long p=s.find_first_not_of(' ');
  if(p==string::npos) return "";
  int q=s.find_last_not_of(' ');
  return TString(s.substr(p,q-p+1));
}


TString printTools::justify(TString s, int w, bool left, bool trm){
  if(trm) s=trim(s);
  int l=s.Length();
  if(l>=w) return s;
  if(left) return s+string(w-l,' ');
  else return string(w-l,' ')+s;
}



TString printTools::leftJustify(TString s,int w,bool trm) {
  return justify(s,w,true,trm);
}

TString printTools::rightJustify(TString s, int w, bool trm){
  return justify(s,w,false,trm);
}



TString printTools::doOrDont(bool b) {return b? "":"don't";}