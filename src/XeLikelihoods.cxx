#include "XeLikelihoods.h"

pdfLikelihood::~pdfLikelihood(){ 

	bkg_components.clear();
	safeguarded_bkg_components.clear();
//	delete data;
	delete asimovData;
//	delete dmData;

}

pdfLikelihood::pdfLikelihood(string name, double wimpMass) : ProfileLikelihood(name) {  
	//setting to a dummy experiment value, this parameter is usefull if you want to combine
	//likelihoods, in that case you need to override it and make sure the two likelihood 
	//have different "experiment" number.
	setExperiment(1);  

	wimp_mass = wimpMass;

	POI = new SigmaParameter();
	
	signal_component = NULL ;

	data = NULL;

	calibrationData = NULL;

	asimovData = NULL;

	dmData = NULL;

	siganlDefaultNorm = NONE;

	XsecMultiplier    = NONE;

	withSafeGuard     = false;

	safeGuardParam    = NULL;

	safeguardAdditionalComponent  = NULL;

}



void pdfLikelihood::checkInputs(){
 	//some checks
	if(bkg_components.size() == 0) {
		cout << "pdfLikelihood::initialize - ERROR: no bkg components are set. Quit. " << endl;
		exit(100);
	}

	if(signal_component == NULL) {
		cout << "pdfLikelihood::initialize - ERROR: no signal component is set. Quit. " << endl;
		exit(100);
	}

	if(data == NULL) {
		cout << "pdfLikelihood::initialize - ERROR: no data are set. Quit. " << endl;
		exit(100);
	}


	if(siganlDefaultNorm == NONE) {
		cout << "pdfLikelihood::initialize - ERROR: no signal reference cross section is set, use setSignalDefaultNorm(Xsec_value).      Quit. " << endl;
		exit(100);
	}

	if(withSafeGuard &&  calibrationData == NULL ) {
		cout << "pdfLikelihood::initialize - ERROR: no calibration sample set, cannot add safeguard.   Quit. " << endl;
		exit(100);
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

       //parameter of interest
       addParameter(POI);

       //nuissance parameters  from bkg
       for(unsigned int k = 0; k < bkg_components.size(); k++){
              cout << "INFO :  adding sys for BKG component "<< bkg_components[k]->getName() << endl; 

	       for(unsigned int j=0; j <  bkg_components[k]->myScaleUnc.size() ; j++){
	       addParameter(  bkg_components[k]->myScaleUnc[j], AUTO );
      		 }

	       for(unsigned int j=0; j <  bkg_components[k]->myShapeUnc.size() ; j++){
	       addParameter(  bkg_components[k]->myShapeUnc[j], AUTO );
      		 }
       }


       //nuissance parameters  from signal
        cout << "INFO :  adding sys for SIGNAL component "<< signal_component->getName() << endl; 
	for(unsigned int j=0; j <  signal_component->myScaleUnc.size() ; j++){
	    addParameter(  signal_component->myScaleUnc[j], AUTO );
      	}

	for(unsigned int j=0; j <  signal_component->myShapeUnc.size() ; j++){
	    addParameter(  signal_component->myShapeUnc[j], AUTO );
      	}



	if(withSafeGuard){
		safeGuardParam = new scaleSys("SafeGuard", 0.);

		// this parameter has no additional constraint
		safeGuardParam->setType(FREE_PARAMETER);  
		addParameter(safeGuardParam, AUTO);
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
		       cout << " pdfLikelihood::addBkgPdfComponent - ERROR: bkg component " << bkg_components[k]->getName() << " already exist. Quit" << endl;
		       exit(100);
		}
	}
        //-----------------------------------------------------------------//


	safeguarded_bkg_components.push_back(Safeguarded);

	bkg_components.push_back(addMe);

	cout << "pdfLikelihood - INFO: bkg component named " << addMe->getName() << " added to " ;//<< getName();
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
	return  (POI->getCurrentValue() * getSignalMultiplier() *  signal_component->getNormalizedEvents());
}

double pdfLikelihood::computeTheLogLikelihood() {

   if(printLevel > 3)   cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : ENTER"<<endl;
  //Retriving Parameter of Interest value
    double sigma = POI->getCurrentValue();

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
   if(printLevel > 3)   cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : STAGE 1  sigma val " << sigma <<endl;

     // LL = log(likelihood) 
     double LL = 0;
    
   //------------- LOAD PDF WITH SYS VARIATION ---------------------//
     TH2F signalPdf(signal_component->getInterpolatedHisto());
//	signalPdf.Rebin2D(8,2);      //REMOVEME

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
		    printf (" OOOHHHHHHHHHHHH background integral : %f %f \n",bkgPdf.Integral(), bkgPdftemp.Integral());
		    exit(100);
		}
     }
     else{
	 //just sum up components otherwise
          bkgPdf = (bkg_components[0]->getInterpolatedHisto());

          for(unsigned int k=1; k < bkg_components.size(); k++){

	            TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());
		    bkgPdf.Add(&temp_bkgPdf);

	    }
//		bkgPdf.Rebin2D(8,2);      //REMOVEME
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

     if(printLevel > 0)    cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : STAGE 3  Ns "<< Ns << " Nobs " << Nobs  << "   Nb " << Nb <<endl;
     if(printLevel > 0)    cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : STAGE 3  sigma "<< sigma << " getSignalMultiplier() " <<getSignalMultiplier()  << "   signalPdf.Integral() " << signalPdf.Integral() <<endl;
     if(printLevel > 1)    cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : STAGE 3  Ns + Nb "<< Ns + Nb <<endl;
     if(printLevel > 1)    cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : STAGE PoissonTerm "<< Nobs * log(Ns + Nb ) -Ns - Nb <<endl;


  
   //------------------------- PDF TERM ----------------------------//
     double extended_term   = 0.;

     Long64_t Nentry  = data->getEntries();  
   
      if(printLevel > 3)    cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : Nentry "<< Nentry <<endl;
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

	//cout << "signalPdf " << signalPdf.GetBinContent(signalPdf.FindBin(ts1,ts2)) << "   sigma  " << sigma << "  getSignalMultiplier  "<<  getSignalMultiplier() << endl;
	     
	     // check physical result 
	     if(extended_signal + extended_bkg <=0) {	     
	     	if(printLevel > 0) cout << "pdfLikelihood::computeTheLogLikelihood - WARNING : NsFs + NbFb <= 0. " << extended_signal  << " <-- NsFs  NbFb--> " << extended_bkg << endl;
	     	return VERY_SMALL;
	     }

	     extended_term += tweight * log( (extended_signal + extended_bkg) / (Ns + Nb) ) ; //data weight is 1 for DM data and whatever for asimov

	     if(printLevel > 0)    cout<<"....................... - INFO : STAGE 3.1  ex_S "<< extended_signal << "   ex_B "<< extended_bkg << "  W " << tweight << " ex-T " << extended_term   <<endl;
     }
     
     LL += extended_term ;
   //---------------------------------------------------------------//
     if(printLevel > 1)  cout<<"pdfLikelihood::computeTheLogLikelihood - INFO : STAGE 4  extendedTerm "<< extended_term   <<endl;



   //-------------------- Adding the safeguard term-------------//
     if(withSafeGuard) {
	     // Check if Safeguard is applicable, Ns > Nb*epsilon otherwise you can eat up your signal
	     // in this case safeguard can bring problem.
	     double epsilon = safeGuardParam->getCurrentValue() ;
	     if(  epsilon <= 0. ) {
	     	    if(printLevel > 1) cout << "pdfLikelihood::computeTheLogLikelihood - WARNING : safeguard not safe " <<  endl;
		    return VERY_SMALL; 
	     }
    	     LL += LLsafeGuard();
     }     
   //---------------------------------------------------------------//


   //----------------------- ADDING NP CONSTRAINTS -----------------//
     TRAVERSE_PARAMETERS(it) {
	 if( (it->second)->getType() == NUISANCE_PARAMETER || (it->second)->getType() == FIXED_PARAMETER) 
		 LL +=   (it->second)->getLLGausConstraint();
     }
   //---------------------------------------------------------------//


  
     if(printLevel > 3)      cout << "computeTheLogLikelihood = " << LL << endl;

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

     double sigma = POI->getCurrentValue();
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


     Info("Model Comparison", ".....");
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
     double epsilon           = safeGuardParam->getCurrentValue();
    
     //computing Nb(1-epsilon)Fb for the safegurded components 
     for(unsigned int k=0; k < bkg_components.size(); k++){
	     //adding up only the safeguarded components
	     if(!safeguarded_bkg_components[k])  continue;

	     TH2F temp_bkgPdf (bkg_components[k]->getInterpolatedHisto());

	     standard_integral += temp_bkgPdf.Integral();

	     // Nb_k *(1 -epsilon) Fb_k(x,y)
	     bkg_plusSafeguard.Add(&temp_bkgPdf, 1. - epsilon);

             // this is the Nb that will multiply Fs
             Nb_safeguard += temp_bkgPdf.Integral() ;
	     
      }

     // Adding Nb*epsilon*Fs     
     TH2F Fs (signal_component->getInterpolatedHisto());	
     bkg_plusSafeguard.Add(&Fs, epsilon * Nb_safeguard / Fs.Integral() ) ;

  
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
	
	     safeguard_only.Add(&temp_bkgPdf);
      }
      
      return safeguard_only;
}



bool pdfLikelihood::isNegativeAnywhere(TH2F histo){


	for(int k=1; k <= histo.GetNcells(); k++){
		if(histo.GetBinContent(k) < 0. ) {
	     		if(printLevel > -1) cout << "pdfLikelihood::computeTheLogLikelihood - WARNING : safeGuard component <= 0 "<< histo.GetBinContent(k) << endl;
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



     //adding the "additional" component: meant to be for AC which is different
     if(safeguardAdditionalComponent) 
	     safeguard_only.Add(safeguardAdditionalComponent);

     double safeguard_only_integral=safeguard_only.Integral();
   if(printLevel > 4)  cout << "safeguard_only_integral=safeguard_only.Integral  = " << safeguard_only_integral << endl;
     //loop over all data
     for(Long64_t event = 0; event < Nentry; event++){
       double ts1=calibrationData->getS1(event);
       double ts2=calibrationData->getS2(event);
       double tweight=calibrationData->getW(event);
       

	     //loads data TNtuple entry
	     //	     calibrationData->getEntry(event);
	     //Nb*Fb(1-epsilon) + epsilon*Nb*Fs
	     //	     double NbFb =  calibrationData->getValFromPdf(safeguard_only);
             //	     double NbFb =  tweight * safeguard_only.GetBinContent(safeguard_only.GetXaxis()->FindBin(ts1), safeguard_only.GetYaxis()->FindBin(ts2));
	     double NbFb =  tweight * safeguard_only.GetBinContent(safeguard_only.FindBin(ts1,ts2));

	     // check physical result 
	     if(NbFb <=0) {	     
	     	if(printLevel > 1) cout << "pdfLikelihood::computeTheLogLikelihood - WARNING : safeGuard component <= 0. " << NbFb << endl;
	     	return VERY_SMALL;
	     }

	     //	     LL +=  log( NbFb / safeguard_only.Integral() ) ;
	     LL +=  log( NbFb / safeguard_only_integral ) ; 

	     if(printLevel > 4)    cout<<"....................... - INFO : Safeguard NB*Fb(x)" << NbFb << "  LL for data i " << LL << endl;
     }
     
     return LL ;

}

