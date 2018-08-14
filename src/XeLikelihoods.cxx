#include "XeLikelihoods.h"

pdfLikelihood::~pdfLikelihood(){

	bkg_components.clear();
	safeguarded_bkg_components.clear();
//	delete data;
	delete asimovData;
//	delete dmData;

}

pdfLikelihood::pdfLikelihood(TString nam, double wimpMass) : ProfileLikelihood(nam) {
	//setting to a dummy experiment value, this parameter is usefull if you want to combine
	//likelihoods, in that case you need to override it and make sure the two likelihood
	//have different "experiment" number.
	setExperiment(1);

	wimp_mass = wimpMass;


	signal_component = NULL ;

	data = NULL;

	calibrationData = NULL;

	asimovData = NULL;

	dmData = NULL;

	siganlDefaultNorm = NONE;

	XsecMultiplier    = NONE;

	withSafeGuard     = false;

	safeGuardPosDef = true;

	safeGuardParam    = NULL;

	safeguardAdditionalComponent  = NULL;

	safeguard_fixValue = -9;

}



void pdfLikelihood::checkInputs(){
 	//some checks
	if(bkg_components.size() == 0) {
		Error("checkInputs","no bkg components are set.");
	}

	if(signal_component == NULL) {
		Error("checkInputs"," no signal component is set. Quit. ");
	}

	if(data == NULL) {
		Warning("checkInputs","no data are set. BETTER YOU KNOW WHAT YOU ARE DOING!");
	}


	if(siganlDefaultNorm == NONE) {
		Error("checkInputs","no signal reference cross section is set, use setSignalDefaultNorm(Xsec_value).");
	}

	if(withSafeGuard &&  calibrationData == NULL ) {
		Error("checkInputs","no calibration sample set, cannot add safeguard.");
	}
}

int pdfLikelihood::numberOfSafeguarded(){

	int N = 0;
	for(unsigned int k=0; k < safeguarded_bkg_components.size(); k++){
		if(safeguarded_bkg_components[k]) N++;
	}

	return N;
}



void pdfLikelihood::initialize(){
       cout << "pdfLikelihood::initialize - INFO :  initialize..... " << endl;

       //some checks
       checkInputs();

	   // clear the parameters (in case one calls initialize more than once)
	   clearTheParameters();

       //parameter of interest
	   LKParameter *POI = new SigmaParameter();
       addParameter(POI);

       //nuissance parameters  from bkg
       for(unsigned int k = 0; k < bkg_components.size(); k++){
              cout << "INFO :  adding sys for BKG component "<< bkg_components[k]->getComponentName() << endl;

	       for(unsigned int j=0; j <  bkg_components[k]->myScaleUnc.size() ; j++){
	       addParameter(  bkg_components[k]->myScaleUnc[j], AUTO );
      		 }

	       for(unsigned int j=0; j <  bkg_components[k]->myShapeUnc.size() ; j++){
	       addParameter(  bkg_components[k]->myShapeUnc[j], AUTO );
      		 }
       }


       //nuissance parameters  from signal
        cout << "INFO :  adding sys for SIGNAL component "<< signal_component->getComponentName() << endl;
	for(unsigned int j=0; j <  signal_component->myScaleUnc.size() ; j++){
	    addParameter(  signal_component->myScaleUnc[j], AUTO );
      	}

	for(unsigned int j=0; j <  signal_component->myShapeUnc.size() ; j++){
	    addParameter(  signal_component->myShapeUnc[j], AUTO );
      	}



	if(withSafeGuard){
		safeGuardParam = new scaleSys(getName() + "_SG", 1.);

		// this parameter has no additional constraint
		safeGuardParam->setType(FREE_PARAMETER);
			if(safeguard_fixValue > 0. ) safeGuardParam->setType(FIXED_PARAMETER);
		safeGuardParam->setInitialValue(1.0);
			if(safeguard_fixValue > 0. ) safeGuardParam->setInitialValue( safeguard_fixValue );
		safeGuardParam->setCurrentValue(1.0);
			if(safeguard_fixValue > 0. ) safeGuardParam->setCurrentValue( safeguard_fixValue );

		if (safeGuardPosDef) safeGuardParam->setMinimum(0.00001);
		else safeGuardParam->setMinimum(-20);
		safeGuardParam->setMaximum(20);
		safeGuardParam->setStep(0.2);
		addParameter(safeGuardParam, AUTO);

		safeguard_scaling = 1000.;
	}
	else if(numberOfSafeguarded() > 0 )
		cout << "\n------ WARNING -------  Safeguard is turned OFF altough you have set components to be safeguarded this is ignored -----\n" << endl;



	//set the multiplier such that the total integral of signal is 1 event.
	setSignalMultiplier( 1. /  signal_component->getDefaultEvents() ) ;


	cout << "pdfLikelihood::initialize - INFO :  initialization SUCCESSFUL. " << endl;

}




void pdfLikelihood::setData(int dataType) {
   switch (dataType) {
	case DM_DATA:
	  useDMData();
	  break;

	case ASIMOV_DATA:
	  useAsimovData();
	  break;

	default :
	  useDMData();
	  cout<< "pdfLikelihood::setData - ERROR - Wrong data Type, set to obsData " << endl;

   }


}


int pdfLikelihood::useAsimovData() {
	  if  (!(asimovData==NULL))
	    {data = asimovData; return (data->getEntries()); }
	  else {printf("XeLikelihoods:: cannot switch to asimovData dataset - it is empty"); return (-1);}
	}

int pdfLikelihood::useDMData() {
	  if  (!(dmData==NULL))
	    {data = dmData; return (data->getEntries()); }
	  else {printf("XeLikelihoods:: cannot switch to dmData dataset - it is empty"); return (-1);}
	}




void pdfLikelihood::addBkgPdfComponent(pdfComponent *addMe , bool Safeguarded){

        // check if component exist ----------------------------------------//
	for(unsigned int k=0; k < bkg_components.size(); k++) {
		if(bkg_components[k] == addMe){
		       cout << " pdfLikelihood::addBkgPdfComponent - ERROR: bkg component " << bkg_components[k]->getComponentName() << " already exist. Quit" << endl;
		       exit(100);
		}
	}
        //-----------------------------------------------------------------//


	safeguarded_bkg_components.push_back(Safeguarded);

	bkg_components.push_back(addMe);

	cout << "pdfLikelihood - INFO: bkg component named " << addMe->getComponentName() << " added to " ;//<< getName();
	if    (Safeguarded) cout << "   SAFEGUARDED" << endl;
	else  cout << "   NOT SAFEGUARDED" << endl;

}

void pdfLikelihood::generateAsimov(double mu_prime) {

        delete asimovData;
 	double scaleFactorSignal =  mu_prime * getSignalMultiplier()  ;

	// get copy of Histo signal default
	TH2F temp_signal( signal_component->getDefaultHisto() )	;

	// get copy of bkg and sum over all, assumes they are normalized.
	// NOTE: safeguard is irrelevant for Asimov dataset
	TH2F temp_histo  (bkg_components[0]->getDefaultHisto());

	for(unsigned int k=1; k < bkg_components.size(); k++){
		//important to avoid returning pointers to local variables
		TH2F very_temp_hist(bkg_components[k]->getDefaultHisto());
		temp_histo.Add(&very_temp_hist);
	}

	temp_histo.Add(&temp_signal,scaleFactorSignal);

	//generate a TH2F as fake data

	//	data->generateAsimov(scaleFactorSignal, &temp_signal, &temp_bkg);

	asimovData =new dataHandler(Form("ASIMOV_DATA_%.2f", mu_prime), &temp_histo); //scaleFactorSignal, &temp_signal, &temp_bkg);


}



void pdfLikelihood::generateToyDataset(double seed, double mu_prime) { cout << "ERROR - pdfLikelihood::generateToyDataset is not supported" << endl;}



double pdfLikelihood::getCurrentNs(){
	return  (getPOI()->getCurrentValue() * getSignalMultiplier() *  signal_component->getNormalizedEvents());
}

double pdfLikelihood::getSafeguardValue(){
	return safeGuardParam->getCurrentValue() / safeguard_scaling;
}

double pdfLikelihood::computeTheLogLikelihood() {

   Debug("pdfLikelihood::computeTheLogLikelihood"," ENTER");
  //Retriving Parameter of Interest value
    double sigma = getPOI()->getCurrentValue();

  //------------------------- CHECK FIT STATUS -----------------------------//
  //Check fit  status: check all parameters value, if any is NaN means the
  //fit is bad. It can happen in case of crazy value of mu, not necessarily
  //means that the limit is bad.
    TRAVERSE_PARAMETERS(it) {
	double param_value = (it->second)->getCurrentValue();
     	if(std::isnan(param_value)){
    		cout<<"pdfLikelihood::computeTheLogLikelihood - WARNING : Fit is unstable"<<endl;
    		return VERY_SMALL;
	}
     }
   //---------------------------------------------------------------//

     // LL = log(likelihood)
     double LL = 0;

   //------------- LOAD PDF WITH SYS VARIATION ---------------------//
     TH2F signalPdf(signal_component->getInterpolatedHisto());

     // get copy of bkg and sum over all, assumes they are normalized.
     TH2F bkgPdf;


     if(withSafeGuard)  {
	     bkgPdf = getSafeguardedBkgPdf() ;

             // MOSHE check here the safeguard implementation
	     TH2F bkgPdftemp;
             bkgPdftemp = (bkg_components[0]->getInterpolatedHisto());

     	     for(unsigned int k=1; k < bkg_components.size(); k++){

       		TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());
       		bkgPdftemp.Add(&temp_bkgPdf);
     	    }

            if (fabs(bkgPdf.Integral()-bkgPdftemp.Integral())>1e-2) {
		    Error("ComputeTheLikelihood", Form("OOOHHHHHHHHHHHH background integral : %f %f \n",bkgPdf.Integral(), bkgPdftemp.Integral()));
		}
     }
     else{
		  Debug("computeTheLikelihood","Adding bkgs:");
	  	  //just sum up components otherwise
          bkgPdf = (bkg_components[0]->getInterpolatedHisto());
		  Debug("computeTheLikelihood", TString::Format("\t%s  = %f events",bkgPdf.GetName(), bkgPdf.Integral()));

          for(unsigned int k=1; k < bkg_components.size(); k++){

	            TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());
		    	bkgPdf.Add(&temp_bkgPdf);
				Debug("computeTheLikelihood", TString::Format("\t%s  = %f events", temp_bkgPdf.GetName(), temp_bkgPdf.Integral()));

	    }
	}
   //---------------------------------------------------------------//




   //----------------------- POISSON TERM --------------------------//
     double Nb = bkgPdf.Integral();
     double   Ns    = sigma * getSignalMultiplier() *  signalPdf.Integral();
     double   Nobs  = data->getSumOfWeights();    // this is == Nentry in case of data, but is not in case of asimov

     //protection against uphysical values of Ns, this is very common in binned case
     //for unbinned should be impossible.
     if( Ns + Nb <= 0.)  {
	      cout << "pdfLikelihood::computeTheLogLikelihood - WARNING : Ns + Nb <= 0." << endl;
	     return VERY_SMALL;
     }

     LL += Nobs * log(Ns + Nb ) -Ns - Nb ;
   //---------------------------------------------------------------//

    Debug("pdfLikelihood::computeTheLogLikelihood" , Form(" Ns %f    Nobs %f   Nb %f ", Ns , Nobs,  Nb));
    Debug("pdfLikelihood::computeTheLogLikelihood" , Form(" Sigma %f    sigmaMultiplier %f  SignalHistoIntegral %f", sigma, getSignalMultiplier() , signalPdf.Integral() ));
    Debug("pdfLikelihood::computeTheLogLikelihood" , Form(" PoissonTerm %f ",Nobs * log(Ns + Nb ) -Ns - Nb ));



   //------------------------- PDF TERM ----------------------------//
     double extended_term   = 0.;

    Long64_t Nentry  = data->getEntries();

    Debug("pdfLikelihood::computeTheLogLikelihood" , Form(" Nentry %lld ", Nentry ));
    //loop over all data
    for(Long64_t event = 0; event < Nentry; event++){
		double ts1=data->getS1(event);
		double ts2=data->getS2(event);
		double tweight=data->getW(event);
	     //loads data TNtuple entry
	     //data->getEntry(event);

	     //double extended_signal =  data->getValFromPdf(signalPdf) * sigma * getSignalMultiplier();
     	     //double extended_bkg    =  data->getValFromPdf( bkgPdf ) ;
		//	     double extended_signal =  signalPdf.GetBinContent(signalPdf.GetXaxis()->FindBin(ts1), signalPdf.GetYaxis()->FindBin(ts2)) * sigma * getSignalMultiplier();
		//    double extended_bkg    =   bkgPdf.GetBinContent(bkgPdf.GetXaxis()->FindBin(ts1), bkgPdf.GetYaxis()->FindBin(ts2));
		double extended_signal =  signalPdf.GetBinContent(signalPdf.FindBin(ts1,ts2)) * sigma * getSignalMultiplier();
		double extended_bkg    =   bkgPdf.GetBinContent(bkgPdf.FindBin(ts1,ts2));

	    Debug("computeTheLogLikelihood", TString::Format("S1 %f  --- S2 %f  ---- weight %f  ---- Fs %f  ----- Fb %f", ts1, ts2, tweight,extended_signal, extended_bkg ));

	    // check physical result
	    if(extended_signal + extended_bkg < 0) {
	      	Warning("pdfLikelihood::computeTheLogLikelihood" , "NsFs + NbFb < 0 ");
	     	return VERY_SMALL;
	    }
		else if( extended_signal + extended_bkg > 0 )  // skipping the case of zero that ahime happens even tough my reccomendations on templates.
	    	extended_term += tweight * log( (extended_signal + extended_bkg) / (Ns + Nb) ) ; //data weight is 1 for DM data and whatever for asimov

		Debug("computeTheLogLikelihood", TString::Format("Extended term  %f", extended_term ));

    }

    LL += extended_term ;
   //---------------------------------------------------------------//
      Debug("pdfLikelihood::computeTheLogLikelihood" , Form(" Extended term %f ", extended_term ));



   //-------------------- Adding the safeguard term-------------//
     if(withSafeGuard) {
	     // Check if Safeguard is applicable, Ns > Nb*epsilon otherwise you can eat up your signal
	     // in this case safeguard can bring problem.
		 double epsilon = safeGuardParam->getCurrentValue() / safeguard_scaling;

		Debug("pdfLikelihood::computeTheLogLikelihood" , Form(" PoissonTerm %f ",Nobs * log(Ns + Nb ) -Ns - Nb ));

	     if(  epsilon <= 0. && safeGuardPosDef ) {
				 Warning("computeTheLogLikelihood", "safeguard not safe");
			     return VERY_SMALL;
	     }
			 LL += LLsafeGuard();
		//Info("pdfLikelihood::computeTheLogLikelihood" , Form(" Safeguard %f ", safeGuardParam->getCurrentValue() ));
     }
   //---------------------------------------------------------------//


   //----------------------- ADDING NP CONSTRAINTS -----------------//
      // this is now moved at higher level due to combination
	  // (in combination one would otherwise consider this term twice)
	  // it is done in likelihood::LikelihoodEvaluate()
   //---------------------------------------------------------------//



	 Debug("computeTheLogLikelihood", Form("LogLike %f", LL));

  return LL;

}







void  pdfLikelihood::drawAllOnProjection(bool isS1Projection){

   histoCompare comp = getModelCompare();

   comp.rebinX = 1;
   comp.rebinY = 10;

   if(isS1Projection) comp.titleX = ("cS1  [PE]");
   else comp.titleX = ("cS2  [PE]");

   comp.doStack =  true;

   comp.projectionX = isS1Projection;

   comp.compareWithRatio();


   TLatex latex;
   latex.SetTextSize(0.05);
   TString t = "";
     TRAVERSE_PARAMETERS(it) {
		 t.Append((it->second)->getName());
		 t.Append(Form("=%.1f ",(it->second)->getCurrentValue()));
     }

    latex.DrawLatex(.2,.9,t);

}



histoCompare pdfLikelihood::getModelCompare(){



     printCurrentParameters();

     double sigma = getPOI()->getCurrentValue();
     double scaleFactorSignal =  sigma * getSignalMultiplier()  ;
     TH2F signalPdf(signal_component->getInterpolatedHisto());
     signalPdf.Scale(scaleFactorSignal);


     TH2F *dataHisto = (TH2F*)signalPdf.Clone("dataHisto");
     dataHisto->Reset();
     data->fillDataHisto(dataHisto);



     histoCompare c;
     c.setBaseHisto(*dataHisto,"Data");

     double bkgIntegral = 0.;

     for(unsigned int k=0; k < bkg_components.size(); k++){
	     TH2F temp_hist(bkg_components[k]->getInterpolatedHisto());

     	     c.addHistoToList(temp_hist,"Background");

	     bkgIntegral +=temp_hist.Integral();
     }

     c.addHistoToList(signalPdf,"WIMP");


    // Info("Model Comparison", ".....");
     cout <<"Data   " << dataHisto->Integral() << endl;;
     cout <<"Bkg    "  << bkgIntegral << endl;;
     if(scaleFactorSignal != 0.)
              cout <<"Signal "  << signalPdf.Integral()/scaleFactorSignal << " POI " << sigma << "  SignalMultiplier  " << getSignalMultiplier() << endl;
     else
	     cout <<"Signal "  << 0 << " POI " << sigma << "  SignalMultiplier  " << getSignalMultiplier() << endl;

    cout << "...... That signal scaled is " << signalPdf.Integral() << endl;


   return c;
}







TH2F pdfLikelihood::getSafeguardedBkgPdfOnly(){

   	if(numberOfSafeguarded() == 0) {
	     cout << "pdfLikelihood::getSafeguardedBkgPdf - ERROR: none of the bkg component is safeguarded " << endl;
	     cout << "\t==> use: addBkgPdfComponent(pdfComponent *addMe , bool Safeguarded = true) \n\n...Quit!\n" << endl;
     	     exit(100);
     }

     TH2F bkg_plusSafeguard = bkg_components[0]->getInterpolatedHisto();
     bkg_plusSafeguard.Reset(); // just to take the skeleton


     float standard_integral = 0. ;
     double Nb_safeguard      = 0. ;
     double epsilon           = safeGuardParam->getCurrentValue() / safeguard_scaling;

     //computing Nb(1-epsilon)Fb for the safegurded components
     for(unsigned int k=0; k < bkg_components.size(); k++){
	     //adding up only the safeguarded components
	     if(!safeguarded_bkg_components[k])  continue;

	     TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());

		 Debug("getSafeguardedBkgPdfOnly",TString::Format("component %s  n-events = %f",temp_bkgPdf.GetName(), temp_bkgPdf.Integral()));
	     standard_integral += temp_bkgPdf.Integral();

	     // Nb_k *(1 -epsilon) Fb_k(x,y)
	     bkg_plusSafeguard.Add(&temp_bkgPdf, 1. - epsilon);

             // this is the Nb that will multiply Fs
             Nb_safeguard += temp_bkgPdf.Integral() ;

      }

     // Adding Nb*epsilon*Fs
	 TH2F Fs (signal_component->getInterpolatedHisto());

	 Debug("getSafeguardedBkgPdfOnly", TString::Format("safeguard_value %f   corresponding to events = %f  and Nb= %f", epsilon, epsilon * Nb_safeguard, Nb_safeguard ));
	 //bkg_plusSafeguard.Add(&Fs, epsilon * Nb_safeguard / Fs.Integral() ) ;
	 plotHelpers::addHisto(&bkg_plusSafeguard, &Fs, epsilon * Nb_safeguard / Fs.Integral() );

	 //cross check:
     if( fabs(bkg_plusSafeguard.Integral() - standard_integral) > 0.001||
		     bkg_plusSafeguard.Integral() <= 0.) {
	  cout <<"pdfLikelihood::getSafeguardedBkgPdf - ERROR: probability is not conserved in safeguard." << endl;
	  cout << bkg_plusSafeguard.Integral() << " !=  " << standard_integral << "\nQuit!. " <<endl;
     	  exit(100);
	}



     return bkg_plusSafeguard;
}



TH2F pdfLikelihood::getSafeguardedBkgPdf(){

     TH2F safeguard_only = getSafeguardedBkgPdfOnly();

     //adding up all the other non safeguarded components
     for(unsigned int k=0; k < bkg_components.size(); k++){
	     if(safeguarded_bkg_components[k])  continue;

	     TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());
		 Debug("getSafeguardedBkgPdfOnly",TString::Format("component %s  n-events = %f",temp_bkgPdf.GetName(), temp_bkgPdf.Integral()));
	     safeguard_only.Add(&temp_bkgPdf);
      }

      return safeguard_only;
}

TH2F pdfLikelihood::getOverallBkgHisto(){

		TH2F bkgPdf = (bkg_components[0]->getInterpolatedHisto());
		for(unsigned int k=1; k < bkg_components.size(); k++){

	            TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());
		    	bkgPdf.Add(&temp_bkgPdf);
	    }

		return bkgPdf;

}

bool pdfLikelihood::isNegativeAnywhere(TH2F histo){


	for(int k=1; k <= histo.GetNcells(); k++){
		if(histo.GetBinContent(k) < 0. ) {
	     	Warning("computeTheLogLikelihood", Form(" safeGuard component <= 0 --> %f ", histo.GetBinContent(k) ));
			return true;
		}
	}

	return false;
}



double pdfLikelihood::LLsafeGuard(){

     double LL = 0;

     Long64_t Nentry = calibrationData->getEntries();

     TH2F safeguard_only = getSafeguardedBkgPdfOnly();

     //Check that the PDF is positive everywhere, gives trouble otherwise.
     /*if(isNegativeAnywhere(safeguard_only)){
	     	return VERY_SMALL;
	}
	*/

	Debug("pdfLikelihood::LLSafeguard" , Form("Calibration NEntry %lli", Nentry));

     //adding the "additional" component: meant to be for AC which is different
     if(safeguardAdditionalComponent) {
	     safeguard_only.Add(safeguardAdditionalComponent);
		 Debug("LLsafeGuard", TString::Format("Additional component integral %f", safeguardAdditionalComponent->Integral()));
	 }

     double safeguard_only_integral=safeguard_only.Integral();
     //if(printLevel > 4)  cout << "safeguard_only_integral=safeguard_only.Integral  = " << safeguard_only_integral << endl;
     //loop over all data
     for(Long64_t event = 0; event < Nentry; event++){
       double ts1=calibrationData->getS1(event);
       double ts2=calibrationData->getS2(event);
       double tweight=calibrationData->getW(event);

	   // Debug("pdfLikelihood::LLSafeguard" , Form("Calibration Event: S1 %f ; S2 %f ; weight %f",ts1, ts2, tweight));

	     //loads data TNtuple entry
	     //	     calibrationData->getEntry(event);
	     //Nb*Fb(1-epsilon) + epsilon*Nb*Fs
	     //	     double NbFb =  calibrationData->getValFromPdf(safeguard_only);
             //	     double NbFb =  tweight * safeguard_only.GetBinContent(safeguard_only.GetXaxis()->FindBin(ts1), safeguard_only.GetYaxis()->FindBin(ts2));
	     double NbFb =  tweight * safeguard_only.GetBinContent(safeguard_only.FindBin(ts1,ts2));

	     // check physical result
	     if(NbFb <=0) {
	     	cout << "pdfLikelihood::computeTheLogLikelihood - WARNING : safeGuard component <= 0. " << NbFb << endl;
	     	return VERY_SMALL;
	     }

	     //	     LL +=  log( NbFb / safeguard_only.Integral() ) ;
	     LL +=  log( NbFb / safeguard_only_integral ) ;

	     //if(printLevel > 4)    cout<<"....................... - INFO : Safeguard NB*Fb(x)" << NbFb << "  LL for data i " << LL << endl;
     }

	Debug("LLsafeGuard", TString::Format("LL safeguard term %f", LL));
    return LL ;

}


void pdfLikelihood::setTreeIndex(int index){

	data->setTreeIndex(index);
	if(withSafeGuard) calibrationData->setTreeIndex(index);
}

pdfComponent* pdfLikelihood::getBkgComponent(TString search_name) {

	for(unsigned int i=0; i < bkg_components.size(); i++){
		if(bkg_components[i]->getComponentName() == search_name) return bkg_components[i];
	}
	Error("getBkgComponent", "component " + search_name + " not found.");

	return NULL;
}

vector<pdfComponent*> pdfLikelihood::getBkgComponents() {

	return bkg_components;
}

void pdfLikelihood::printEventSummary(bool isForWiki){

	if(!isForWiki)
	{
		cout << "\n\n--------------Event Summary------------------\nPdfComponent Name \tEvents"<< endl;

		for(unsigned int i=0; i < bkg_components.size(); i++){
			cout << TString::Format("%s \t %1.4f",bkg_components[i]->getComponentName().Data(), bkg_components[i]->getNormalizedEvents()) << endl;
		}

    	cout << TString::Format("Signal \t %1.4f", signal_component->getNormalizedEvents()) << endl;
		if(data !=NULL) data->printSummary();
	}
	else {

		cout << "\n\n--------------Event Summary For Wiki------------------\nPdfComponent Name \tEvents"<< endl;

		cout << "^ ";

		for(unsigned int i=0; i < bkg_components.size(); i++){
			cout << TString::Format(" %s ^",bkg_components[i]->getComponentName().Data()) ;
		}
		cout << " Signal ^ Data ^" << endl;
		cout << "| ";

		for(unsigned int i=0; i < bkg_components.size(); i++){
			cout << TString::Format(" %1.4f |", bkg_components[i]->getNormalizedEvents()) ;
		}

		int entries = 0;
		if(data != NULL) entries = data->getEntries();
	    cout << TString::Format(" %1.4f | %d |", signal_component->getNormalizedEvents(), entries) << endl;
	}

}
