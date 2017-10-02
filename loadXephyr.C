int loadXephyr(){
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/XeCore.cxx+g");
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/XeMath.cxx+g");
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/XeStat.cxx+g");
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/dataHandler.cxx+g");
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/XePdfObjects.cxx+g");
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/XeLikelihoods.cxx+g");
  gROOT->ProcessLine(".L $XEPHYR_DIR/src/AsymptoticExclusion.cxx+g");
  return 0;
}
