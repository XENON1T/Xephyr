#ifndef XE_UTILS
#define XE_UTILS

#include "TString.h"
#include <iostream>

using namespace std;

/**
  * \class errorHandler
  * \brief Helper Class to send consistently formatted messages. It allows different levels  of severity for the messages: Info, Warning and Error.

  * Error actually throws an  error with the std::runtime_error() which will 
  * abort (no so) gracefully and quit  Xephyr, whatever it is doing.
*/

class errorHandler {
   public:
    //! Constructor: @param name: Tag name, typically here you wanna put the
    //! name of the class for which you wanna handle errors.
	errorHandler(TString name);

    //! \brief Send an Error message to the user
    //
    //! Send a message with ERROR severity, Xephyr will quit.
    //! @param functionName: Tag identifier, function or in general reference place 
    //! in the code where this message has been produced. 
    //! @param message: message to the user.
    void Error(TString functionName, TString message);

    //! \brief Send a Warning message to the user
    //
    //! Send a message with WARNING severity, Xephyr will continue.
    //! @param functionName: Tag identifier, function or in general reference place 
    //! in the code where this message has been produced. 
    //! @param message: message to the user.
    void Warning(TString functionName, TString message);


    //! \brief Send an Info message to the user
    //
    //! Send a message with INFO severity, Xephyr will continue.
    //! @param functionName: Tag identifier, function or in general reference place 
    //! in the code where this message has been produced. 
    //! @param message: message to the user.
	void Info(TString functionName, TString message);
 
    //! \brief Send a Debug message to the user
    //
    //! Send a message with DEBUG severity, Xephyr will continue.
    //! @param functionName: Tag identifier, function or in general reference place 
    //! in the code where this message has been produced. 
    //! @param message: message to the user.
	void Debug(TString functionName, TString message);


    //! \brief Set the local printing level (for this instance only)
    //
    //! The print level is defined globally, if you want to modify, increase
    //! reduce the verbosity of this particular instance you can set its print level.
    //! @param level: the print level you wanna set. Default is < 0 which means all,
    //! 0 = Debug + Info + Warning + Error, 1 = Info + Warning + Error, 2 = Warning + Error, 3 = Error only.
    //! you can modify the global print level with errorHandler::globalPrintLevel = XX
    void setPrintLevel(int level) {localPrintLevel = level;};
    
    //! \brief Return the current print level for this instance. 
    //! remind that if you set with setPrintLevel this will override the global print level.
    int  getPrintLevel()  {return (localPrintLevel > 0 ? localPrintLevel : globalPrintLevel) ;};

   
    int localPrintLevel;              //! local print level for this instance modify it with setPrintLevel.
    
    static int globalPrintLevel;      //! global print level of all objects, modify this with errorHandler::globalPrintLevel.
                                      // and you'll modify all.
    
    TString className;	
};



class printTools {


public:

     
    printTools();

     static TString  upperCase(TString s);
     static TString  lowerCase(TString s);
     static TString  trim(TString s);
     static TString  justify(TString s, int w, bool left,bool trim);
     static TString  rightJustify(TString s, int w, bool trim=false);
     static TString  leftJustify(TString s,int w, bool trim=false);
     static TString  format0I(int v, int w=1);
     static TString  formatLI(int v, int w=1, int p=1);
     static TString  formatRI(int v, int w=1, int p=1);
     static TString  formatF(double v, int w=1, int p=1);
     static TString  formatR(double v, int w=1, int p=1);
     static TString  formatLF(double v, int w=1, int p=1);
     static TString  formatRF(double v, int w=1, int p=1);
     static TString  formatG(double v, int w=1, int p=1);
     static TString  formatLG(double v, int w=1, int p=1);
     static TString  formatRG(double v, int w=1, int p=1);
     static TString  formatI(int v, int w=1, int trailer=0);
     static TString  doOrDont(bool b) ;


};




#endif