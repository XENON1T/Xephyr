int loadXephyr(){

  gROOT->ProcessLine(".L src/XeVersion.cxx+g");
  gROOT->ProcessLine(".L src/XeCore.cxx+g");
  gROOT->ProcessLine(".L src/XeMath.cxx+g");
  gROOT->ProcessLine(".L src/XeUtils.cxx+g");
  gROOT->ProcessLine(".L src/XeStat.cxx+g");
  gROOT->ProcessLine(".L src/dataHandler.cxx+g");
  gROOT->ProcessLine(".L src/XePdfObjects.cxx+g");
  gROOT->ProcessLine(".L src/XeLikelihoods.cxx+g");
  gROOT->ProcessLine(".L src/AsymptoticExclusion.cxx+g");
  return 0;
}
