{

// Build the directory
	
 	TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
	gSystem->Exec("mkdir -p " + xeDir + "xephyr_examples/likelihood1D/data");
 	TString inputDir = xeDir + "xephyr_examples/likelihood1D/data/";

	gStyle->SetOptTitle(kFALSE);
      	gStyle->SetOptStat(kFALSE);
	
	TFile *outFile =  new TFile(inputDir + "modelsAndData.root", "RECREATE");

// Production of the models
	
	// BACKGROUND histo definition
	      // NOTE: X axis with 100 bins (for example S1), Y axis with one dummy bin only
	TH2F *background =    new TH2F("bkg_sigma_0.00","best guess background", 100, 3., 100., 1, 0., 1.);
	TH2F *background_p1 = new TH2F("bkg_sigma_1.00","bkg +1 sigma", 100, 3., 100., 1, 0., 1.);
	TH2F *background_m1 = new TH2F("bkg_sigma_-1.00","bkg -1 sigma", 100, 3., 100., 1, 0., 1.);

	double bkg_expected = 600.; // events
	     
	// GENERATING the bkg model with shape uncertainty
	// the bkg model will be just a Gaussian and the shape uncertainty will be
	// defined as a change in variance of the Gaussian.
	TRandom3 rambo;
	double mean = 30.; 	    // Gaussian with mean at 30 PE and sigma 10.
	double sigma = 10.;
	double systematic = 0.20;   // 20% change in Gauss sigma
	double sigma_p1 = sigma * (1. + systematic);
	double sigma_m1 = sigma * (1. - systematic);
	int N = 100000;

	for(int i=0; i <N; i++){
		
		background->Fill( rambo.Gaus(mean, sigma), 0.5 );  // Y is filled at 0.5, just in the middle of the bin
		background_p1->Fill(rambo.Gaus(mean, sigma_p1), 0.5 );
		background_m1->Fill(rambo.Gaus(mean, sigma_m1), 0.5);
	}

	// SCALING histogram to expected number of events (here the same for each sys)
	background->Scale( bkg_expected / ((double)N) ); 
	background_p1->Scale( bkg_expected / ((double)N) ); 
	background_m1->Scale( bkg_expected / ((double)N) ); 

	// plotting
	background_m1->ProjectionX("p2")->Draw("PLC hist");
	background->ProjectionX("p")->Draw("same PLC hist");
	background_p1->ProjectionX("p1")->Draw("same PLC hist");

	// SIGNAL production (without shape uncertainty)
	TH2F *Signal =    new TH2F("Signal","signal", 100, 3., 100., 1, 0., 1.);
	double signal_expected = 4.; // events for a given cross section

	double mean_signal = 59.; // S1
	double sigma_signal = 3.;

	for(int i=0; i <N; i++) Signal->Fill( rambo.Gaus(mean_signal, sigma_signal), 0.5);
	
	Signal->Scale( signal_expected / ((double)N) );
	Signal->ProjectionX("s")->Draw("same PLC hist");



// DATA production
	
        TTree toyTree ("exampleTree", "generated toy sample");
        float cs1 = 0.; 
        float cs2 = 0.;
        string type = "DummyLabel";
        
        toyTree.Branch("cs1",&cs1,"cs1/F");
        toyTree.Branch("cs2",&cs2,"cs2/F");
        toyTree.Branch("type",&type);

	// Let's produce a dataset made of background with a small signal injection of 2 events
	int n_bkg = rambo.Poisson( bkg_expected ) ;
	int n_signal = 30. ; // injection of 30 events

	// LOOP over bkg events
	for(int j =0; j < n_bkg; j++) {
		double temp_cs1, temp_cs2;    // GetRandome returns double.
		background->GetRandom2( temp_cs1 , temp_cs2 );  // the extraction of cs2 here is irrelevant.
		cs1 = (float) temp_cs1;
		cs2 = (float) temp_cs2;
		type = "bkg";
		toyTree.Fill();
	}	

	// LOOP over signal events
	for(int j =0; j < n_signal; j++) {
		double temp_cs1, temp_cs2;    // GetRandome returns double.
		Signal->GetRandom2( temp_cs1 , temp_cs2 );
		cs1 = (float) temp_cs1;
		cs2 = (float) temp_cs2;
		type = "signal";
		toyTree.Fill();
	}	

	toyTree.Draw("cs1","","same PMC PLC E");	

	gPad->BuildLegend();

// SAVING into a file
	outFile->cd();
	toyTree.Write();
	Signal->Write();
	background_p1->Write();
	background->Write();
	background_m1->Write();



}
