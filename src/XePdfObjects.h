#ifndef XEPDFOBJECTS_H
#define XEPDFOBJECTS_H

#include "XeStat.h"
#include <TString.h>
#include "TGraph2D.h"
#include "TNtuple.h"
#include <cmath>
#include <csignal>
#include <iostream>
#include <vector>
#include <utility>      // std::pair, std::make_pair
#include "dataHandler.h"
#include "XeUtils.h"
#include "TColor.h"

using namespace std;


/**
 * \class shapeSys 
 * \brief it represents a Shape Systematic. 
 * It inherits from LKParameter. Here you can use the LKParameter functions to set initial val,
 * min, max and setp
*/
class shapeSys : public LKParameter {

   public:

	/** 
	 * Constructor for class shapeSys.
	 * @param name: Name of the sys as appears in histogram name, see pdfComponent for a detailed use case.
	 */
	shapeSys(TString name);


	//! returns the nearest HISTOGRAM value in the grid, lower edge
	double getNearestLow(); 

	//! returns the nearest HISTOGRAM value in the grid, upper edge
	double getNearestHigh();
	
};


class scaleSys : public LKParameter {

	public:
		// relativeUncertainty = 1 sigma relative uncertainty (example 10% uncertainty would be 0.1) 
		scaleSys(TString name, double relativeUncertainty);
///		~scaleSys();

		// gives the Normalization scale factor according to t-valued uncertainty  --> (1 + t_val * relativeUncertainty ) 
		// Uses the current t-value which is automatically linked to its pdfLikelihood 
		double getNormModifier();

		
		//this case return a sys that is centered in zero, half a gaussian, t-value strictly positive 
		//and for a tvalue=0 have zero events. Histo is supposed to be normalized to 1
		void setNull(){	isNull = true; setMinimum(0.);};

	private:
		double relUnc;
		bool isNull;
	

};


class pdfComponent :errorHandler{

   public:

         pdfComponent(TString name, TString filename);
	~pdfComponent();


	void addScaleSys(scaleSys *addMe) { myScaleUnc.push_back(addMe); };

	void addShapeSys(shapeSys *addMe) { myShapeUnc.push_back(addMe); };

	//! load histogram according to the current value of the parameters
	void loadHistos();

	//! load default histogram, no sys.
	void loadDefaultHisto();

	//! returns the interpolated pdf over shape sys computed in (s1,s2).
	/** 
	 * it is normalized to number of events and scale uncertainty are also taken into account
	 */
	double getNormalizedDensity(double s1, double s2);

	//! returns the (s1,s2) bin content of the default histo, no shape nor scale sys is applied
	double getDefaultDensity(double s1, double s2);

	//! Returns the total integrated number of events taking into account shape sys and scale sys.
	double getNormalizedEvents();
	
	//! returns total integral of the default histo, no shape nor scale sys is applied
	double getDefaultEvents();

	//!Returns a copy of the interpolated histogram according to the "currentValue" of the shapeSys associated with it.
	/** 
	 * The current value of each uncertainty can
	 * be set simply by shapeSys::setCurrentValue(val) method.
	 */
	TH2F  getInterpolatedHisto();

	//! Returns a copy of the default histogram, no shape nor scale sys is applied
	TH2F  getDefaultHisto();

	//! returns the grid points in histogram space that needs to be loaded. this method need to be modified for arbitrary number of shape sys.
	TString getNearestHistoName(vector<bool> setOfVal);

	//! returns the default name of the histogram
	TString getDefaultHistoName();

	//! returns a string with current values of parameters
	TString getParamValueString();

	TString getParamValueWritable();

	//! returns the integral of the square between points (not bins).
	/**
	 * linearly interpolates the area if desn't correspond to an integer number of bins.
	 * Uses the default histo, does not consider shape or scale uncertainty.
	 */
	double getDefaultPdfIntegral(double s1_min, double s1_max, double s2_min, double s2_max);

	//! This scans the parameter space given a number of steps for each parameter
	/** 
	 * produce a projection given min and max, produce a PDF file for comparison
	 * of all interpolation and saves also the TH2F to file.
	 */
	void plotInterpolatedSpace(bool doProjectionX, double min, double max, int Nsteps, bool legend_left = false, double y_max = 0.);

	//! scale the pdf by VAL, this happend for the methods:
	/* getDefaultHisto, getInterpolatedHisto, getNormalizedDensity, 
	 * getDefaultDensity, getNormalizedEvents, getDefaultEvents
	 */
	void setScaleFactor(double val) {scaleFactor = val;};

	//! \brief Normalize the pdf in a way that its integral is nevents
	void setEvents(double nevents);  	

	scaleSys* getScaleSys(TString name);
	shapeSys* getShapeSys(TString name);
	void      replaceUncertainty(TString name, scaleSys* newScale);
	void      replaceUncertainty(TString name, shapeSys* newShape);

	TString getName()  {return pdf_name;};

	vector<scaleSys*> 		myScaleUnc; //! container of scale sys
	vector<shapeSys*> 		myShapeUnc; //! container of shape sys

	TString                         suffix;
	bool                            doExtend;

   private:
	TFile	          		*file;
	TH2F                            *defaultDistro;
  	vector<TH2F*>			histos;	       /** contains the 2^N histo for the hyperplane interpolation of shapeSys */
  	vector<double>		InterpFactors;  /** contains the 2^N interpolation factors for the hyperplane interpolation of shapeSys	*/
	TString 			pdf_name;
	vector<double>			old_t_val;    /** contains the last value interpolated, the interpolation is lazy, doesn't ricompute it if is for the same set of values.*/
	double                          scaleFactor;


	void extendHisto(TH2F &h);

};



/**
 * \class histoCompare.
 * \brief This class is supposed to compare a "base" histo 
 * with a set of other histogram. One can choose for ratio plot,
 * stack plot of different component, one can set the projected axis
 * and the slice on which do the projection, one can also rebin the histo.
 * the histo are supposed to be TH2F.
 */
class histoCompare :errorHandler{

    public:

	    histoCompare();
	    ~histoCompare();

      //---------OPTIONS-----------// 
	    int      rebinX;
	    int      rebinY;
	    bool     projectionX;
	    double   projectionMin;    /** minimum value on the projected axis */
	    double   projectionMax;    /** maximum value on the projected axis */
	    int      binMin;		/** same but in bin nuber */
	    int      binMax;            
	    bool     doStack;   /** True if you want to sum up histo in compareList */
	    TString  titleX;
	    TString  titleY;
       //---------------------------//

	    //! set histo with wich one one to compare the rest
	    void setBaseHisto(TH2F b, TString n=""){base = b; names[0] = n;};
	    //! set the comparison histos
	    void addHistoToList(TH2F h, TString n ="") {compareList.push_back(h); names.push_back(n);};
	    //! set Name of component
	    void setNameofComponent(unsigned int i, TString n);
	    //! draw just a normal comparison plot
	    void compare();
		
		//! print model content in csv format
		void printModels();

	    //! add a ratio plot to compare()
	    void compareWithRatio();

	    void drawLegend(TH1D *baseH, vector <TH1D*> list );
	    
	    //return a string of info
	    TString projectionInfo();

        //------ histo holders --------//
	    TH2F               base;
	    vector <TString>   names;       
	    vector <TH2F>      compareList;

	    TH1D               *projectedBase;
	    vector <TH1D*>      projectedList;

     private:

	    void setOptions(TH1D *h, bool dataLike = true, bool isBottom = false);
	    void project(); /** push histo into the TH1D*/


	   


};






#endif
