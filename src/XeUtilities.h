#ifndef XeUtilities_h
#define XeUtilities_h
#include "XeCore.h"

/**
   *  A set of utilities in random order and choice
*/

class XeUtilities : public XeCore {

  public :

 ~XeUtilities();
  XeUtilities();

  static void copyBands(string input, string output);

  void copyGraphs(string input, string output, vector<string> graphs);

};
#endif

