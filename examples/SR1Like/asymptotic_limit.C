{



	// loading the likelihood definition
	 gROOT->ProcessLine(".L src/likelihoodDef.cxx");

	 pdfLikelihood *pl = getDMLikelihood(50., 0);



	 pl->printEventSummary(false); // print a summary of the event content of the histograms
	 			       // if the argument is true, will print with wiki format

	 pl->maximize(false);  // maximize the likelihood (remember we injected signal 30 events)
	 		       
		 	      



}

