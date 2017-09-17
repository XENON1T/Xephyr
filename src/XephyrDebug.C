int XephyrDebug(){
  gROOT->ProcessLine(".L XeVersion.cxx+g");
  gROOT->ProcessLine(".L XeCore.cxx+g");
  gROOT->ProcessLine(".L XeMath.cxx+g");
  gROOT->ProcessLine(".L XeStat.cxx+g");
  gROOT->ProcessLine(".L dataHandler.cxx+g");
  gROOT->ProcessLine(".L XePdfObjects.cxx+g");
  gROOT->ProcessLine(".L XeLikelihoods.cxx+g");
//  gROOT->ProcessLine(".L XeAnalysis.cxx+g");
//  gROOT->ProcessLine(".L HistLikelihood.cxx+g");
  gROOT->ProcessLine(".L AsymptoticExclusion.cxx+g");
  return 0;
}
