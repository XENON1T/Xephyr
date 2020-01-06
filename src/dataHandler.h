#ifndef DATAHANDLER 
#define DATAHANDLER 
#include "XeStat.h"

#include "XeUtils.h"
#include <TString.h>
#include "TGraph2D.h"
#include "TNtuple.h"
#include "TParameter.h"
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
	    dataHandler( TString name, TString fileName, TString dmTree,
					 TString x_name="x", TString y_name="y", TString z_name="z");

	    // if you just want to run Sensitivity
	    dataHandler( TString name); 

	    //If you want to run with asimov as data
	    dataHandler( TString name, TH3F *h3pdf); 

	    // if you want ro run with fake data as data.
	    dataHandler( TString name, TH3F *h3pdf, int N); 

	    ~dataHandler();
	    
		void initialize();

	    TTree *DMdata;

	    //	    TH3F    *asimovData;

	    TFile   *file;

	    TString  Name;

		TString  TreePrefix;

	    TH3F    *grid;

	    TString x_name;
	    TString y_name;
		TString z_name;

		Float_t x;
	    Float_t y;
		Float_t z;
	    Float_t weight;
		TNtuple *weights;
	    Float_t sumOfWeights;

	    int dataType;

	    void  drawxy(TString opt="");
	    

	    vector<int> getSimulatedInfo(unsigned int size);
	    
	    Float_t getX(long N);
	    Float_t getY(long N);
		Float_t getZ(long N);
		Float_t getW(long N);
	    
	    TGraph getXY();
		
	    Long64_t getEntries();

	    Float_t getSumOfWeights();

	    void printSummary();

	    void     getEntry(Long64_t entry);

	    void     generateAsimov( TH3F *background);

	    void generateDataSet(TH3F *h3pdf, int N);
	
	    static Float_t   integrate(TH3F *histo, Float_t x_min, Float_t x_max, Float_t y_min, Float_t y_max, Float_t z_min, Float_t z_max);
	   
	    void setData(int Type) {dataType = Type; };

	    void fillDataHisto(TH3F *hist); 

	   Float_t getValFromPdf( TH3F &histo );	   
	   
	   //change the tree pointer to data, assumed to be in same file.
	   void   setDataTree(TString nameTree); 

		//! \brief change the tree pointer to data.
		void   setDataTree(TTree *tree); 

	   //change the  file.
	   void   setFile(TString PathtoFile); 

	   void setXname(TString name) {x_name = name;};
	   void setYname(TString name) {y_name = name;};
	   void setZname(TString name) {z_name = name;};

	   //! \brief set the prefix for the name of the tree collection to iterate on, usefull for neyman construction.
	   void setPrefixTree( TString prefix ) { TreePrefix = prefix ; } ;

		//! \brief load the tree in file according to the name convention: "TreePrefix"+"index"
		void setTreeIndex( int index );	

	   //change the  file and tree in one go.
	   void   setFileAndTree(TString PathtoFile, TString nameTree); 

	   void addToDataSet(TH3F *h3pdf, int N);
	  
	   //! \brief useful to get the truth generated value of parameter stored in tree
	   vector<string> getTrueParamsNames();

	   //! \brief useful to get the truth generated value of parameter stored in tree
	   vector<double> getTrueParams();


};





#endif
