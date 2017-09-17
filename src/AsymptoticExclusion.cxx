#include "AsymptoticExclusion.h"

AsymptoticExclusion::AsymptoticExclusion(ProfileLikelihood *pl, double confidenceLevel): XeObject(), 
	expectedMedian(1), expected2Sigma(1), observedLimit(1),observedLimitNoCLS(1),sigmaScan(1),qTestScan(1), 
	qTestScanData(1), pulls_uncond("pulls_uncond","",1,0.,20.), pulls_cond("pulls_cond","",1,0.,20.),
	Histo2D_mass_vs_x("Histo2D_mass_vs_x","",1001,-5.,10005.,61,-2.5,302.5)

{ 

	// a more serious error handling should be put in place... no time now... FIXME Ale

	if(pl == NULL)  cout <<" ProfileLikelihood is empty ---- this will crash!" << endl;

	name = "AsymptoticExclusion_" + pl->getName();

	plike=pl;
        cl = confidenceLevel;

	nScanPoints = 100;
 	
	stored3Sigma_mu = -999.;
	storedMedian_mu = -999.;

	scanMin = UNDEFINED;
	scanMax = UNDEFINED;
	
	useQtilde =  false;	

	AlternativeX = -999;

	Obslimit=-999;
	ObslimitnoCLS=-999;
	
}


AsymptoticExclusion::~AsymptoticExclusion() {}


void AsymptoticExclusion::generateAndSetAsimov(double mu_prime) {

 plike->generateAsimov(mu_prime); // trigger the run to generated data and store it under ASIMOV_DATA


 plike->setData(ASIMOV_DATA); 	  // set the data to use as observed data.

}


void AsymptoticExclusion::setRealData() {

	 plike->setData(DM_DATA);
}


void AsymptoticExclusion::setToyDataset(double seed, double mu_prime){

	plike->generateToyDataset(seed, mu_prime);  //trigger the run to generated data and store in DM_SIMULATED_DATA
	plike->setData(DM_SIMULATED_DATA);

}


double AsymptoticExclusion::computeQTestStat(double mu, bool useStoredFit){
 
  double logD   = UNDEFINED; 		 // Handler for the denominator of q
  double mu_hat = UNDEFINED;		 // Handler for MLE of signal strenght
  double qStat  = UNDEFINED;		 // Handler on the outcome

  // if unconditional fit has never been done, or you want to recompute it
  if(plike->getSigmaHat() == UNDEFINED || !useStoredFit) {
	logD   = plike->maximize(false);  // NOTE: this function automatially stores the maximum 
					  // likelihood and sigma_hat in plike.
	mu_hat = plike->getSigmaHat();
  }

  //If the unconditional fit has been done and you don't want to recompute it, usefull for limits.
  if(useStoredFit) {
	logD    = plike->getLogD();
	mu_hat  = plike->getSigmaHat();
  } 

 

  //in case one wants to use Q-tilde test statistic (formula 16 of http://arxiv.org/abs/1007.1727)
  if(useQtilde){
	if(mu_hat < 0.) {
	    cout << "AsymptoticExclusion::computeQTestStat  -- INFO :  mu_hat in negative!  mu_hat = " << mu_hat << endl;
	    cout << "AsymptoticExclusion::computeQTestStat  -- INFO :  QTilde test statistic was choosen. Setting mu_hat to zero..." << endl;

	    //substitute the denominator with conditional fit mu=0
	    plike->setParameterValue(PAR_SIGMA, 0.);
	    logD = plike->maximize(true);

	    //set mu hat to zero
	    mu_hat = 0;
	}
  }
 
  //Do the conditional fit
  plike->setParameterValue(PAR_SIGMA, mu);  // set signal strenght to mu
  double logN = plike->maximize(true);	    // conditional fit freezing PAR_SIGMA parameter	


  // Computing Q Test Statistic:

 if(mu_hat <= mu) {
	qStat = 2.* (logD - logN );
  }
  
  else if( mu_hat > mu ) qStat = 0.;
 
 //  cout << "logD " << logD << "  logN " << logN << "  qStat " << qStat << "  mu_hat " << mu_hat << endl;
  return qStat;

}


double AsymptoticExclusion::computeSigmaAsimov(double mu_prime){

  generateAndSetAsimov(mu_prime);	// tell plike to set reference data to Asimov



  //LOOP: scan over mu values and compute several times sigma_A, this depends sightly
  //on the value of mu, so we take an average. here we go from mu[0,10], but maybe better
  //to do as a function of step size on NP. FIXME
  double sigma_A =  0.;
  double sigma_temp = 0.;
  double counter= 0.;
  double sigma_min = 999.;
  double sigma_max = 0.;
  double qTest_temp = -999.;

  double XsecScale = plike->getSignalMultiplier() * plike->getSignalDefaultNorm();

  for(double mu_step= mu_prime +0.05; mu_step <= 4. + mu_prime; mu_step = mu_step +0.1){ // the step sise is only preliminary... FIXME

	if(mu_step == mu_prime)  continue;  // this case brings singularity 0/0
	
	qTest_temp = computeQTestStat(mu_step);
  	sigma_temp = fabs(mu_step - mu_prime)/sqrt(qTest_temp);

	if(sigma_min > sigma_temp) sigma_min = sigma_temp;
	if(sigma_max < sigma_temp) sigma_max = sigma_temp; // just to check the spans of sigma

	sigma_A += sigma_temp;

	counter = counter +1. ;
	//  	cout << "sigma_A " << sigma_A << "  counter "<< counter<< endl;
	sigmaScan.SetPoint(sigmaScan.GetN() , mu_step, sigma_temp * XsecScale);
	qTestScan.SetPoint(qTestScan.GetN() , mu_step,qTest_temp);

  }


  // throw a warning in case sigma values differ for more than 10%
  if(sigma_max / sigma_min > 1.1) cout << "AsymptoticExclusion::computeSigmaAsimov - WARNING : sigmas differs by more than 10% " << sigma_max << "  " << sigma_min << endl;
  if(sigma_max / sigma_min > 2)   { 
	  cout << "AsymptoticExclusion::computeSigmaAsimov - ERROR : sigmas differs by more than 100% " << sigma_max << "  " << sigma_min << "  Put smaller scale factor for signal ------> Quit "<< endl;
//	  exit(100);
  }
  if(counter > 0) sigma_A = sigma_A / counter;  // just average  

  return sigma_A ;
}


double AsymptoticExclusion::conditionalFit(double mu){
  //Do the conditional fit
  plike->setParameterValue(PAR_SIGMA, mu);  // set signal strenght to mu
  return plike->maximize(true);	    // conditional fit freezing PAR_SIGMA parameter	

}


double AsymptoticExclusion::unconditionalFit(){

	cout << " ------- START-------- \n AsymptoticExclusion::unconditionalFit() " << endl;
	plike->resetParameters();
	return plike->maximize(false);

}

double AsymptoticExclusion::getMu_hat(){

return plike->getSigmaHat();

}
double AsymptoticExclusion::computeSensitivityHagar(){

    //This computes the limit in mu, not in Xsection, you need to scale	
   
    generateAndSetAsimov(0.);	// tell plike to set reference data to Asimov

    // valid for qtilde and q	
    return plike->returnLimitHagar(cl) ; 

}


void AsymptoticExclusion::computeSensitivity() {

	//sigma asimov for Hypothesis mu=0
//	double sigma_0 = computeSigmaAsimov(0.);   // Old way UNCOMMENT TO GET BACK TO PREVIOUS


	//Signal Strenght multiplier to get back to the "right" cross section
	double XsecScale = plike->getSignalMultiplier() * plike->getSignalDefaultNorm();
	
	//According to the followi("formula the expected CLs limits with bands are:
	// N_sigma_band = +-N*sigma_0 + sigma_0 * Inverse_PHI(1- cl*PHI(+-N)) 
	// where PHI is the cumulative of Normal(0,1) from [-Inf : x], and cl is confidence level

	double ex_POI      =  computeSensitivityHagar();  //COMMENT
	double sigma_0     =  ex_POI / sqrt(ROOT::Math::normal_quantile(1-cl,1)); //COMMENT 
	// the new way
	double ex_hagar      = XsecScale * ex_POI; //COMMENT

	double ex_median      = XsecScale * computeExpectedLimit(sigma_0, 0., cl); //UNCOMMENT

	double ex_err_band1   = XsecScale * computeExpectedLimit(sigma_0, 1., cl) - ex_median;
	double ex_err_band_m1 = ex_median - computeExpectedLimit(sigma_0, -1., cl) * XsecScale;

	double ex_err_band2   = XsecScale * computeExpectedLimit(sigma_0, 2., cl) - ex_median;
	double ex_err_band_m2 = ex_median - computeExpectedLimit(sigma_0, -2., cl) * XsecScale ;
	
	cout << "ex_median " << ex_hagar << "  ex_err_band High " << ex_err_band1 << " ex_err_band Low " << ex_err_band_m1 << endl;

	double M = getMass();
	if(AlternativeX != -999) M = AlternativeX;  // to set an alternative value instead of mass

 	expectedMedian.SetPoint(0, M , ex_hagar);     //COMMENT
 	//expectedMedian.SetPoint(0, M , ex_median);  //UNCOMMENT
 	expected2Sigma.SetPoint(0, M , ex_hagar);    //COMMENT
 	//expected2Sigma.SetPoint(0, M , ex_median);  //UNCOMMENT
	expectedMedian.SetPointError(0, 0., 0., ex_err_band_m1, ex_err_band1);
	expected2Sigma.SetPointError(0, 0., 0., ex_err_band_m2, ex_err_band2);
	
	storedMedian_mu = ex_median / XsecScale;  // this is for limit computation range
	stored3Sigma_mu = 3. * ex_err_band1 / XsecScale;  // this is for limit computation range
}

double AsymptoticExclusion:: computeExpectedLimit(double sigma_0, double N, double CL){
 
   double cdf_of_sigma = ROOT::Math::normal_cdf(N, 1.);
   double limit = N*sigma_0 + sigma_0* ROOT::Math::normal_quantile( 1. -CL*cdf_of_sigma, 1.); 
   //   double limitnocls = sigma_0 * ( N + ROOT::Math::normal_quantile(1-CL,1));
   //cls place holder hyl
   return limit;
}

void AsymptoticExclusion::writeToFile(TString prefix) {

   TString nameFile = "";
   double mass = getMass();

   nameFile =  TString(Form("%1.F",mass));
/*   if(mass >= 10. )   nameFile =  formatI(mass , 2);
   if(mass >= 100. )  nameFile = formatI(mass , 3);
   if(mass >= 1000. ) nameFile = formatI(mass , 4);
   if(mass >= 10000. ) nameFile = formatI(mass , 5);
*/

   if(AlternativeX != -999)   
	   nameFile.Append("_"+(TString::Itoa((int) AlternativeX,10)));  // to set an alternative value instead of mass

   

   TFile limitFile(prefix +nameFile+".root", "RECREATE");

   expectedMedian.Write("expectedMedian");
   expected2Sigma.Write("expected2Sigma");
   observedLimit.Write("observedLimit");
   observedLimitNoCLS.Write("observedLimitNoCLs");
  
   sigmaScan.SetTitle("sigma scan for: " + name); 
   qTestScan.SetTitle("qTestScan scan for: " + name); 
   qTestScanData.SetTitle("qTestScan scan on Data for: " + name); 
   sigmaScan.GetXaxis()->SetTitle("#sigma_{test} #times 10^{-45} cm^{2}");
   sigmaScan.GetYaxis()->SetTitle("#sigma_{asimov} #times 10^{-45} cm^{2}");
   qTestScan.GetXaxis()->SetTitle("#sigma_{test} #times 10^{-45} cm^{2}");
   qTestScan.GetYaxis()->SetTitle("q");
   qTestScanData.GetXaxis()->SetTitle("#sigma_{test} #times 10^{-45} cm^{2}");
   qTestScanData.GetYaxis()->SetTitle("q");
   
   sigmaScan.Write("sigmaScan_"+nameFile);
   qTestScan.Write("qTestScan_"+nameFile);
   qTestScanData.Write("qTestScanData_"+nameFile);

   pulls_uncond.Write("pulls_uncond_"+nameFile);
   pulls_cond.Write("pulls_cond_"+nameFile);

   if(AlternativeX != -999)   {
	   Histo2D_mass_vs_x.GetXaxis()->SetTitle("M_{\\chi} [GeV]");
	   Histo2D_mass_vs_x.GetYaxis()->SetTitle("User Alternative Variable");
	   Histo2D_mass_vs_x.Write();
   }

  limitFile.Close();
}

void AsymptoticExclusion::computeSetAsimovSigma() {

   //this function generates a set of sigma Asimov values for the Hypotesis bkg only
   // computed for different mu value, exactly the same mu values that are used for the 
   // scan in the limits.

   asimovSigmaSet.clear();

   muStepsSet.clear();
  
   generateAndSetAsimov(0.);	// tell plike to set reference data to Asimov
    

   double range =  2.*stored3Sigma_mu;  // range is from ~ -3sigma to +3sigma
   double step = range / ((double) nScanPoints);

   double min = storedMedian_mu - stored3Sigma_mu; 
   double max = storedMedian_mu + stored3Sigma_mu;

   if(scanMin!=UNDEFINED && scanMax!=UNDEFINED){  //one can also define his own range
	min  = scanMin;
	max  = scanMax;
	step = (max - min) / ((double) nScanPoints);
   }
    

   if(min <0. ) { min = 0.; cout <<"AsymptoticExclusion::computeSetAsimovSigma() - WARNING min sigma scan < 0. " << endl; }

   for(double mu_point = min; mu_point <= max; mu_point = mu_point + step){

	muStepsSet.push_back(mu_point);

	if(mu_point == 0.) { asimovSigmaSet.push_back(UNDEFINED); continue; }
	
	
	double qTest_temp = computeQTestStat(mu_point);
	asimovSigmaSet.push_back( mu_point/sqrt(qTest_temp) );

	//	cout << "mu_val = " << mu_point << "  sigma_A = " << mu_point/sqrt(qTest_temp) << endl;
   }



}



void AsymptoticExclusion::computeLimits() {
  Obslimit=-999;
  ObslimitnoCLS=-999;
  
  // we need to define first what is the expected limit so it is easy to scan for mu and find best limit
  //computeSensitivity();

  // get the set of asimov sigma for bkg hipotesis for Pb computation. Best to do it here.
 // computeSetAsimovSigma();   // UNCOMMENT TO RESTORE

  // Set the likelihood data to be the DM data (not asimov)
  setRealData(); 

  // Do unconditional fit once
  double logL_uncond_fit = unconditionalFit();

  
  pulls_uncond = plike->getPullsHisto();    
  plike->printCurrentParameters();



   /// LOOP over values of mu within a defined range, the range and steps are defined in computeSetAsimovSigma() 
   double q_obs = -999.;
   double mu_limit = -999.;

   double q_obs_nocls = -999.;
   double mu_limit_nocls = -999.;

   qTestScanData.Set(0);


   mu_limit_nocls = plike->returnLimitHagar(0.1, true);  //COMMENT
/*
      //UNCOMMENT
   for(unsigned int i=0; i < muStepsSet.size(); i++){

	double mu_point =  muStepsSet[i];
	double sigmaA   =  asimovSigmaSet[i]; 

	//observed value of test stat. for defined hypotesis of mu
	//NOTE: in case of UseQtilde == true, then we are using wrong cumulative formulaes.  FIXME
	double temp_q_obs = computeQTestStat(mu_point, true); //NOTE: here we use stored unconditional fit just computed 
	double pval_S_B = compute_pval_s_plus_b(temp_q_obs);
	double pval_B   = compute_pval_b(temp_q_obs, mu_point, sigmaA);

	// While you do the scan push everything to qTestScanData Tgraph
	qTestScanData.SetPoint(qTestScanData.GetN() , mu_point,temp_q_obs);

	// Continue the scan until you reached two points:
	// 1).  no cls limit - that is P(s+b)=cl  (output saved into mu_limit_nocls, q_obs_nocls and later to (double) ObslimitNoCLS and (TGraph) observedLimitNoCLS
	if (pval_S_B < cl && mu_limit_nocls==-999)
	  {  mu_limit_nocls = mu_point;
		q_obs_nocls    = temp_q_obs;
	  }

	// 2).  cls limit - that is P(s+b)/Pb=cl  (output saved into mu_limit, q_obs and later to (double) Obslimit and (TGraph) observedLimit 
	if(pval_S_B / pval_B < cl && mu_limit==-999){
		mu_limit = mu_point;
		q_obs    = temp_q_obs;
	} 
	if (mu_limit!=-999 && mu_limit_nocls!=-999) break;
	
   }	

  */

   mu_limit  = mu_limit_nocls;

   Obslimit  = mu_limit *  plike->getSignalMultiplier() * plike->getSignalDefaultNorm();
   ObslimitnoCLS  = mu_limit_nocls *  plike->getSignalMultiplier() * plike->getSignalDefaultNorm();

   pulls_cond = plike->getPullsHisto(); // get pulls after conditional fit for mu value = mu_limit

   if(mu_limit < 0. ) cout <<"AsymptoticExclusion::computeLimits()  -- WARNING : Observed limit out of scanned range!! default range is within 3 sigma from expected limit. You may want to change this using SetScanMin and SetScanMax" << endl;
   else { cout << "OBSERVED limit not scaled mu " << mu_limit 
	  <<"  cross section no cLS " << ObslimitnoCLS
	  << "  test stat value no cls " << q_obs_nocls
          <<"  cross section  " << Obslimit 
     	  << "  test stat value " << q_obs << endl;

        double M = getMass();
	if(AlternativeX != -999) M = AlternativeX;  // to set an alternative value instead of mass

	//Store Observed limit in Graph
   	observedLimit.SetPoint(0, M, Obslimit);
	observedLimitNoCLS.SetPoint(0, M, ObslimitnoCLS);

	//Fill 2D histo
	if(AlternativeX != -999) Histo2D_mass_vs_x.Fill(getMass(),AlternativeX,Obslimit);
   }



}





double AsymptoticExclusion::compute_pval_s_plus_b(double q_obs) {
	return ROOT::Math::normal_cdf_c(sqrt(q_obs), 1.);  
}

double AsymptoticExclusion::compute_pval_b(double q_obs, double mu_val, double sigma_val) {

	if(sigma_val == UNDEFINED || sigma_val == 0.) return ROOT::Math::normal_cdf_c(sqrt(q_obs), 1.);

	double x = sqrt(q_obs) - mu_val / sigma_val;

	return ROOT::Math::normal_cdf_c( x , 1.);


}

void  AsymptoticExclusion::LikelihoodScan(){

  // Set the likelihood data to be the DM data (not asimov)
  setRealData(); 

  // Do unconditional fit once
  double logL_uncond_fit = unconditionalFit();



   /// LOOP over values of mu within 
   double q_obs = -999.;

   for(double mu_step = 0.; mu_step < 5.; mu_step = mu_step + 0.1 ){


	//observed value of test stat. for defined hypotesis of mu
	double temp_q_obs = computeQTestStat(mu_step, true); //NOTE: here we use stored unconditional fit just computed 
	
        double sigmaTemp = mu_step * plike->getSignalMultiplier() * plike->getSignalDefaultNorm();
	qTestScanData.SetPoint(qTestScanData.GetN() , mu_step,temp_q_obs);
	} 
	



}



void  AsymptoticExclusion::setAlternative2DHistoRange(int Nx, double xmin, double xmax, int Ny, double ymin, double ymax){

	Histo2D_mass_vs_x.SetBins(Nx, xmin, xmax, Ny, ymin, ymax);

}



