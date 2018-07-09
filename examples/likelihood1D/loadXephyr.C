int loadXephyr(){
 
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));

  gROOT->ProcessLine(".L " + xeDir +"/Xephyr/src/XeVersion.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/XeUtils.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/XeStat.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/dataHandler.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/XePdfObjects.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/XeLikelihoods.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/AsymptoticExclusion.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/ToyGenerator.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/plotHelpers.cxx+g");
  gROOT->ProcessLine(".L  " + xeDir +"/Xephyr/src/ToyFitterExclusion.cxx+g");
  gInterpreter->AddIncludePath( xeDir +"/Xephyr/src"); // in this case is just XEPHYR src from next dir.
  return 0;
}
