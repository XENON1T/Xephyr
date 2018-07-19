{
	
	
	//===================================================//
	// 		SIMPLE script to run asymptotic limits		 //
	//===================================================//




	// loading the likelihood definition
	gROOT->ProcessLine(".L src/likelihoodDef.cxx");
        
	
	// this is to load a combined likelihood of 3 volumes, "other" is default to zero
	// and is not needed for SR1, but could be usefull as additional parameter for your signal model
	// See signalDef.h  	
	// CombinedProfileLikelihood* getDMCombinedLikelihood( double mass , double other = 0.) {
	
	CombinedProfileLikelihood *pl = getDMCombinedLikelihood(50. );
	
	// this is to load a likelihood only for one volume [0=egg, 1=U-volume, 2=wall]
	// other as above is dummy. 
	// pdfLikelihood* getDMLikelihood(double mass , int volume, double other=0 );
	
	// pdfLikelihood *pl = getDMLikelihood(50., 0);

	// print a summary of the event content of the histograms
	// if the argument is true, will print with wiki format
	pl->printEventSummary(false); 
	 			       

	// This creates the handler for Asymptotic limit production
	// This method uses a one-sided test statistic
	// 0.05 stands for 95% CL, we choose 95% so that we can compare with 90% CL 
	// produced with a two-sided test static.
	AsymptoticExclusion ae (pl, 0.05); 


	ae.computeLimits();

	ae.writeToFile("asymtpotic_limit_");


}

