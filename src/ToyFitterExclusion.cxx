#include "ToyFitterExclusion.h"


ToyFitterExclusion::ToyFitterExclusion(TString CollectionName):errorHandler("ToyExclusion"){


    likeHood = NULL;
    dirPath  = "./";
    OutDir = "./";

    treeName = CollectionName;
    calTreeName = "";
    Suffix      = "";
    likelihood_uncond = 0.;
    likelihood_cond = 0.;
    limit_converged = false;
    testStat_limit  = 0.;
    name_params.clear();
    CurrentTreeIndex = 0;
    IndexHolder = -9;
    mu_fit = 0.;
    testStat =0.;
    q_tilde  = 0.;
    numberOfParams = 0.;
    randomizeMeasure = true;
    graph_of_quantiles = NULL;
    mu_limit = 0.;
    limit    = 0.;
    lower_limit = -1.;
    lower_mu_limit = -1.; 
    testStat_at0 = 0.;
    Gen = UNDEFINED_INT;
}

void ToyFitterExclusion::for_each_tree( double (ToyFitterExclusion::*p2method)(double), TTree *outTree,  double mu, int stopAt){
    

    // setting up branches on outTree
    mu_fit    = mu;
    testStat = 0.;
    q_tilde  = 0.;
    int inputTreeIndex = -9;
    numberOfParams =  likeHood->getParameters()->size();
    double mass    =  likeHood->getWimpMass();
    outTree->Branch("mu_fit", &mu_fit, "mu_fit/D");
    outTree->Branch("mass", &mass, "mass/D");
    outTree->Branch("q_mu", &testStat, "q_mu/D");
    outTree->Branch("q_tilde", &q_tilde, "q_tilde/D");
    outTree->Branch("n_params", &numberOfParams, "n_params/I");
    outTree->Branch("LL_cond", &likelihood_cond, "LL_cond/D");
    outTree->Branch("LL_uncond", &likelihood_uncond, "LL_uncond/D");
    outTree->Branch("true_params", true_params,"true_params[n_params]/D");
    outTree->Branch("measured_params", measured_params,"measured_params[n_params]/D");
    outTree->Branch("uncond_params", uncond_params,"uncond_params[n_params]/D");
    outTree->Branch("cond_params", cond_params,"cond_params[n_params]/D");
    outTree->Branch("name_params", &name_params);
    outTree->Branch("name_true_params", &name_true_params);
    outTree->Branch("inputTreeIndex", &inputTreeIndex, "inputTreeIndex/I");
    outTree->Branch("generation", &Gen, "generation/I");

    // reset CurrentTreeIndex
    CurrentTreeIndex = 0;

    while ( CurrentTreeIndex < stopAt ) {
        
        Info("fit", TString::Format("Fitting tree index, %d", CurrentTreeIndex));

        inputTreeIndex = CurrentTreeIndex;

        // fill true values from generated tree
        fillTrueParams(); 

        // set toy data for fit (except case is only one dataset, like you are fitting data)
        if(stopAt != 1) likeHood->setTreeIndex(CurrentTreeIndex);

        // reset the parameter to their nominal initial value (othrwise takes longer to fit)
        likeHood->resetParameters();

        // Dice the measured parameters (note that if a parameter is free the t-value exytacted here has no effect)
        // Note: this MUST be called after "fillTrueParams"
        if(randomizeMeasure) measureParameters();

        // Fancy coding isn't it? ;)  
        // This is a functional: using a pointer to a function of ToyFitterExclusion
        // so that we can run this same loop for different purposes
        testStat = (this->*p2method)(mu);
        
        outTree->Fill();

        CurrentTreeIndex++;
    }


}

void ToyFitterExclusion::fit(double mu, int stopAt){
    
    TFile f_out(OutDir + "post_fit_" + treeName + Suffix + ".root","RECREATE");

    // output tree, here intentionally all out tree will have the same name so we can hadd
    TTree *outTree = new TTree("post_fit_tree", "output tree for a given mu, hadd me");

    // read each tree in input file "f" nad applies the computeTS function to it 
    for_each_tree( &ToyFitterExclusion::computeTS, outTree, mu, stopAt );
        
    f_out.cd();
    outTree->Write();
    f_out.Close();

}



void ToyFitterExclusion::spitTheLimit(TGraphAsymmErrors *ninety_quantiles, int stopAt){
    
    TFile f_out(OutDir + "limits_" + treeName + ".root","RECREATE");

    graph_of_quantiles = ninety_quantiles;

    // output tree, here intentionally all out tree will have the same name so we can hadd
    TTree *outTree = new TTree("limit_tree", "tree containing limits, hadd me");
    
    // attach a few additional branch to out tree
    mu_limit = 0.;
    limit    = 0.;
    lower_limit = -1.; 
    lower_mu_limit = -1.; 
    testStat_at0 = 0.;
    outTree->Branch("mu_limit", &mu_limit, "mu_limit/D");
    outTree->Branch("lower_mu_limit", &lower_mu_limit, "lower_mu_limit/D");
    outTree->Branch("limit", &limit, "limit/D");
    outTree->Branch("lower_limit", &lower_limit, "lower_limit/D");
    outTree->Branch("testStat_limit", &testStat_limit, "testStat_limit/D");
    outTree->Branch("testStat_at0", &testStat_at0, "testStat_at0/D");
    outTree->Branch("limit_converged", &limit_converged, "limit_converged/O");

    // read each tree in input file "f" and applies the computeTS function to it 
    for_each_tree( &ToyFitterExclusion::limitLoop, outTree, -9., stopAt );
              
    f_out.cd();
    outTree->Write();
    f_out.Close();
    // f->Close();
    
}

double ToyFitterExclusion::limitLoop(double initial_mu){

    // PROCEDURE: 
    // 1- Find mu_hat and TS value for mu = 0.
    // 2- if TS(mu=0.5) > 90% quantile alternative hypo @ mu=0.5 ==> we need to produce an interval
    // 3- Find the lower limit if needed
    // 4- Find the upper limit

    // some checks
    if(graph_of_quantiles == NULL) Error("limitLoop", "graph of quantiles is not defined.");

    limit_converged = false;
    testStat_limit  = 0.;
    lower_limit     = -1.;
    lower_mu_limit = -1.; 
    bool do_interval = false;

    // Finding TS @ mu=0
    testStat_at0 = computeTS( 0. ) ;
    double mu_hat = likeHood->getSigmaHat();
    // check if we cross @ low mu
    if(testStat_at0 > graph_of_quantiles->Eval(0.) ) do_interval = true;

    /////////////////// FIND mu region using Asymptotics best guess  /////////////////////////
    //                                                                                      //
    double one_over_sigma = computeTS( mu_hat + 1. ) ;

    double mu_asym_best = getBestAsympoticGuessForMu( one_over_sigma, mu_hat);
    
    Info("limitLoop",TString::Format("Asymptotics best guess for mu = %1.2f",mu_asym_best));
    //                                                                                      //
    //////////////////////////////////////////////////////////////////////////////////////////

    // from ROOT docs:  https://root.cern.ch/root/html/ROOT__Math__BrentMinimizer1D.html
    //                  https://root.cern.ch/numerical-minimization
    //                  https://root.cern.ch/how-implement-mathematical-function-inside-framework
    ROOT::Math::Functor1D func(this,&ToyFitterExclusion::eval_testStatMinuit); 
    ROOT::Math::BrentMinimizer1D bm; // this guy does numerical minimization scanning a range, not the fastest but quite robust.
    bm.SetNpx(3);                    // divide in 3 intervals
    
    if(do_interval) {
        // CASE IN WHICH WE COMPUTE FC INTERVALS 
        Info("limitLoop","WE GOT INTERVALS THIS TIME! Computing lower limit.");
        bm.SetFunction(func, 0., mu_hat );                 // interval from [0, mu_hat] in mu 
        limit_converged = bm.Minimize(10, 0.05, 0.01);     // max 10 iteration absolute error on mu = 0.05 relative error 1% 
        lower_mu_limit = bm.XMinimum();
        lower_limit = lower_mu_limit * likeHood->getSignalMultiplier() * likeHood->getSignalDefaultNorm();
    }

    // NORMAL CASE UPPER LIMIT
    double mu_max =  mu_asym_best + mu_asym_best*0.5 + 3.; // usually mu_asym_best should not be wrong by more than 20% 
    Info("limitLoop",TString::Format("Limit computed in the range of mu [%1.2f, %1.2f]",mu_hat, mu_max ));

    bm.SetFunction(func, mu_hat , mu_max );   // interval from [mu_hat,15] in mu (by construction limit should be around 2.3)
    limit_converged = bm.Minimize(10, 0.05, 0.01);     // max 10 iteration absolute error on mu = 0.05 relative error 1% 

    double q_stat  =  bm.FValMinimum();  // test stat value at limit

    mu_limit = bm.XMinimum();
    limit    = mu_limit * likeHood->getSignalMultiplier() * likeHood->getSignalDefaultNorm();

    Info("limitLoop", TString::Format("%s with TS val= %1.3f", (limit_converged ? "CONVERGED" : "NOT - CONVERGED" ), q_stat ) );
    Info("limitLoop",TString::Format("computed for Tree index %d ---> mu_limit= %f (~events) xsec = %E cm^2", CurrentTreeIndex, mu_limit, limit));
    
    return q_stat;
}


double ToyFitterExclusion::eval_testStatMinuit( double mu )  {

    // computing the difference between H0 qstat for a given mu and H_mu qstat.
    double qstat = computeTS(mu);
    double delta  = graph_of_quantiles->Eval(mu) - qstat;
    double last_mu = (graph_of_quantiles->GetX())[graph_of_quantiles->GetN() - 1];
     
    if(mu > last_mu) delta  = graph_of_quantiles->Eval(last_mu) - qstat; // for large value where the graph is not defined  no need to extrapolate
    
    // the plus 50 is to make it asymmetric with respect to zero, 
    // it had a lot of difficulties in finding the right minima otherwise for some cases.
    // The problem might be solved (or at least minimized)  with more than 3 point scan.
    // The initial scan values might be such that going to lower mu get lower delta, but the crossing, teh zero is actually
    // for higher mu.
    // This is a work around, a better solution might be neded FIXME.
    if( delta > 0.7 ) delta = delta + 50.;   // 0.7 so that it will be a bit faster and will be symmetric around 0.
    
    Debug("limitLoop", TString::Format("current_mu = %f   ;  current_qstat = %f  ;  delta = %f" ,mu, qstat, delta));        
    
    testStat_limit  =  qstat; // assuming the last call will set the qstat_limit properly
    
    return fabs(delta);
}


double ToyFitterExclusion::getBestAsympoticGuessForMu(double one_over_sigma_squared, double mu_hat){
    
    // according to asymptotic formulae:
    // mu_lim = mu_hat + sigma * sqrt(90_quantile)

    double quantile_90 = ROOT::Math::chisquared_quantile(0.9, 1);
    double sigma = sqrt( 1. /  one_over_sigma_squared );

    double asym_limit =  mu_hat + sigma * sqrt(quantile_90) ;

    return asym_limit;
}

double ToyFitterExclusion::computeTS(double mu) {

    // MAYBE make sense to brake this in two functions conditional and unconditional fit
    // instead of this messy thing here, see fore example above limit_loop (that sucks) FIXME.
    
    // we use a twosided TS t_tilde from equation 11 of: https://arxiv.org/abs/1007.1727
    // furthermore we report also q_tilde from eq. 16.
    // both are suitable for upper limit in which the parameter of interest must be >0 
    
    double qstat = 0.;
    
    // this is for limit case, where we are running many fit with different mu_test
    // on the same data tree. In those cases mu_hat is always the same, you don't want
    // to do unconditional fit again.
    bool DoMaximize  = ( IndexHolder != CurrentTreeIndex );
    Debug("computeTS", TString::Format("maximize Tree index %d ", CurrentTreeIndex));
    
    if(DoMaximize) IndexHolder = CurrentTreeIndex;

    // performing unconditional fit and saving it in unconditional likelihood  for outTree
    if(DoMaximize)    likelihood_uncond = likeHood->maximize(false) ;
    double mu_hat  = likeHood->getSigmaHat();   // the likelihood store safely mu_hat and overrides it only when you call maximize(false)

    double LL_denominator = likelihood_uncond;

    // printing 
    if(getPrintLevel() < WARNING) 
        likeHood->printCurrentParameters();
    
    // save parameters of unconditional fit 
    if(DoMaximize)  saveParameters(uncond_params);
    
    // override denominator in case mu_hat < 0.
    if( mu_hat < 0.){
        likeHood->getParameter(PAR_SIGMA)->setCurrentValue(0.);    
        LL_denominator = likeHood->maximize(true) ;
    }

    // perform conditional fit
    likeHood->getParameter(PAR_SIGMA)->setCurrentValue(mu);
    double LL_numerator = likeHood->maximize(true) ;
    
    // saving the conditional likelihood for outTree
    likelihood_cond = LL_numerator;
    
    // printing 
    if(getPrintLevel() < WARNING) 
        likeHood->printCurrentParameters();

    // save parameters of conditional fit 
    saveParameters(cond_params);

    qstat = 2. * ( LL_denominator - LL_numerator );  // -2*logLikelihoodRatio

    // following the q_tilde construction, this could have been put above with just return 0.
    // I just wanted to do the conditional fit anyway to return conditional parameter for study
    q_tilde = qstat;
    if( mu_hat > mu )  q_tilde = 0. ;
    return qstat;
}


 
void ToyFitterExclusion::saveParameters(double *outParams){

    map <int, LKParameter*> *params = likeHood->getParameters();
    
    int itr = 0;
    name_params.clear();
    
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
        
        LKParameter *param = ip->second;
        outParams[itr] = param->getCurrentValue();
        name_params.push_back(param->getName().Data()); // this is a bit stupid to do here since it does it twice (for cond and uncond)        
        itr++;
    }    
}


void ToyFitterExclusion::saveNames(string *names){
    
        map <int, LKParameter*> *params = likeHood->getParameters();
        
        int itr = 0;
        for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
            
            LKParameter *param = ip->second;
            names[itr] = param->getName();
            itr++;
        }    
}
    


void ToyFitterExclusion::fillTrueParams(){

    // retrive the previously saved list of parameters (done in ToyGenerator)
    name_true_params.clear();
    name_true_params = likeHood->getTrueParamsNames();

    vector<double> par = likeHood->getTrueParams();

    for(unsigned int i=0; i < par.size(); i++){
        true_params[i] = par[i];
    }
}



void ToyFitterExclusion::measureParameters(){

    // randomize measure of parameters around the true value.
    
    map <int, LKParameter*> *params = likeHood->getParameters();
    
    Info("measureParameters", "Randomizing parameters:");

    int parItr = -1;
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
            
        LKParameter *param = ip->second;
        parItr++;   // start from zero 
        
        if(param->getType() == FIXED_PARAMETER || param->isOfInterest() 
            || param->getType() == FREE_PARAMETER ) {

             measured_params[parItr] = param->getCurrentValue();
             Info("---->","Skipping paramater: " + param->getName()); 
             continue; 
        }
            
        // getting range
        double min = param->getMinimum();
        double max = param->getMaximum();
            
        // sample gauss tvalue for parameter
        double random_tvalue = 0.;
        random_tvalue = rambo.Gaus(true_params[parItr],1.);
        // extract again if out of range
        while(random_tvalue > max || random_tvalue < min)
                    random_tvalue = rambo.Gaus(true_params[parItr],1.);
    
    
        param->setT0value(random_tvalue);
        param->setCurrentValue(random_tvalue);  // you wanna start the fit from here
        measured_params[parItr] = random_tvalue;
        Info("---->",TString::Format("new T0-Value: %s = %1.2f",param->getName().Data(),random_tvalue)); 
    }

}


TGraphAsymmErrors ToyFitterExclusion::computeTSDistros(TString fileName, double *mu_list, int mu_size){


    TFile *f = new TFile(fileName);
    
    TTree *tree = (TTree*) f->Get(treeName);
 
    return computeTSDistros(tree, mu_list, mu_size );
    
        
}


TGraphAsymmErrors ToyFitterExclusion::computeTSDistros(TTree *tree, double *mu_list, int mu_size){

    double wimpMass = likeHood->getWimpMass();
    return plotHelpers::giveTSquantiles(tree,mu_list, mu_size, OutDir, wimpMass);

}


