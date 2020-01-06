{

// simple script for some likelihood games.
	
	// loading the likelihood definition
	 gROOT->ProcessLine(".L src/likelihoodDef.cxx");

	 pdfLikelihood *pl = getTheLikelihood(); // return a pointer to the builded likelihood

	 pl->initialize();  // initialize the likelihood, this is necessary. If you want to add
	 		    // more parameter you need to do it before than this.

	 pl->printEventSummary(false); // print a summary of the event content of the histograms
	 			       // if the argument is true, will print with wiki format

	 pl->maximize(false);  // maximize the likelihood (remember we injected signal 30 events)
	 		       // try to change that in "produceDataAndModels.C" and see what happens
		 	      // if argument is true will do the fit with fixed signal strength

	// pdfComponent *bkg = pl->getBkgComponent("bkg");  // return the pdfComponent called "bkg" from likelihood
	
	// // this produce a PDF file plot with an example of 10 interpolation points of the shape uncertainty
	// // see https://xenon1t.github.io/Xephyr/class_reference/classpdfComponent.html#a181072320088088fd8e84a5d970fc847
	// // for more details on this method 
    // bkg->plotInterpolatedSpace (true, 0, 100, 10); 



    // // class example to compare histograms after fit
    // histoCompare *p = pl->getModelCompare();
    // p.setPrintLevel(DEBUG);
    // p.setNameofComponent(1, "background");
    // p.setNameofComponent(2, "Signal");
    // p.rebinX = 2;                     // rebin the axis to make it look nicer ;)
    // p.doStack = true;		      // stack histo on top of each other	
    // p.titleY = "Entries/[PE]";	  
    // p.titleX = "cS1 [PE]";
    // p.projectionMin = 0.;            // minimum Y from which project 
    // p.projectionMax = 1.;	     // maximum Y from which project 
    // p.projectionX = true;            // which axis project X (true) or Y (false)
    // p.compareWithRatio();                     // Draw stuff

    // p.printModels();                 // print content of each bin

}
