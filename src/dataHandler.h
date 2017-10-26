#ifndef DATAHANDLER 
#define DATAHANDLER 
#include "XeStat.h"

#include "XeUtils.h"
#include <TString.h>
#include "TGraph2D.h"
#include "TNtuple.h"
#include <cmath>
#include <csignal>
#include <iostream>
#include <vector>
#include <utility>      // std::pair, std::make_pair

enum DATA_TYPE { 

  UNSPECIFIED_DATA = - 1
   , AM_BE_DATA               /*!< Neutron calibration data*/
   , E_GAMMA_DATA             /*!< Gamma calibration data*/
   , DM_DATA                  /*!< Dark matter data*/
   , ER_BACKGROUND            /*!< Gaussian ER background*/
   , ER_LEAKAGE               /*!< Anomalous leakage*/
   , NR_BACKGROUND            /*!< Neutron background*/
   , ALL_BACKGROUNDS          /*!< All backgrounds*/
   , DM_SIMULATED_DATA        /*!< Dark matter simulated data, generated from BKG only, 
				you have to fill this calling XeRun::generateData()*/
   , ASIMOV_DATA	     /*!, Asimov dataset generated with the function XeRun::generateAsimovData(mu)*/	
   } ;

using namespace std;




class dataHandler : public errorHandler{

    public:
	    dataHandler( TString name, TString fileName, TString dmTree);

	    // if you just want to run Sensitivity
	    dataHandler( TString name); 

	    //If you want to run with asimov as data
	    dataHandler( TString name, TH2F *h2pdf); 

	    // if you want ro run with fake data as data.
	    dataHandler( TString name, TH2F *h2pdf, int N); 

	    ~dataHandler();
	    
	    TTree *DMdata;

	    
	    //	    TH2F    *asimovData;

	    TFile   *file;

	    TString  Name;

	    TH2F    *grid;

	    TString FirstVarName;
	    TString SecondVarName;


	    float s1;
	    float s2;
	    float weight;

	    int dataType;

	    void  drawS1S2(TString opt="");

	    TGraph2D *gs1s2w;

	    double sumOfWeights;

	    vector<int> getSimulatedInfo(unsigned int size);
	    
	    double getS1(int N) { if (N>gs1s2w->GetN()) {printf ("ERROR %d larger than Entries (%d) \n",N,gs1s2w->GetN()); return 0;} else return (gs1s2w->GetX()[N]); }
	    double getS2(int N) { if (N>gs1s2w->GetN()) {printf ("ERROR %d larger than Entries (%d) \n",N,gs1s2w->GetN()); return 0;} else return (gs1s2w->GetY()[N]); }
	    double getW(int N) { if (N>gs1s2w->GetN()) {printf ("ERROR %d larger than Entries (%d) \n",N,gs1s2w->GetN()); return 0;} else return (gs1s2w->GetZ()[N]); }
	    
	    TGraph getS1S2();
		
	    Long64_t getEntries();

	    double getSumOfWeights();

	    void printSummary();

	    void     getEntry(Long64_t entry);

	    void     generateAsimov( TH2F *background);

	    void generateDataSet(TH2F *h2pdf, int N);
	
	    static double   integrate(TH2F *histo, double s1_min, double s1_max, double s2_min, double s2_max);
	    void setData(int Type) {dataType = Type; };



	    void fillDataHisto(TH2F *hist); 

	   double getValFromPdf( TH2F &histo );	   
	   
	   //change the tree pointer to data, assumed to be in same file.
	   void   setDataTree(TString nameTree); 

		//! \brief change the tree pointer to data.
		void   setDataTree(TTree *tree); 

	   //change the  file.
	   void   setFile(TString PathtoFile); 
	   //change the  file and tree in one go.
	   void   setFileAndTree(TString PathtoFile, TString nameTree); 

	   void addToDataSet(TH2F *h2pdf, int N);


};





#endif
