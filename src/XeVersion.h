#ifndef XeVersion_h
#define XeVersion_h

#include <iostream>
#include <string.h>
#include <RVersion.h>

#if ROOT_VERSION_CODE < ROOT_VERSION(5,34,20)
 #define STATIC_CONST  static const
#else
 #define STATIC_CONST static constexpr
#endif

using namespace std;


class XeVersion {
  public :
 ~XeVersion(){}
  XeVersion(){
    string tab="      ";
    string XephyrVersion = "6.0";
    string XephyrDate    = "10-Oct-2017";
    cout<<endl<<endl
        <<tab<<"=================== Welcome to Xephyr ====================="
        <<endl
        <<tab<<" Version "<<XephyrVersion<<" dated "<<XephyrDate
             <<", compiled with ROOT "<<ROOT_RELEASE<<endl
        <<tab<<"==========================================================="
        <<endl<<endl;
  }

};

#endif
