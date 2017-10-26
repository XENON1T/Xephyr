#include "ToyFitterExclusion.h"


ToyFitterExclusion::ToyFitterExclusion(TString fileName):errorHandler("ToyExclusion"){


    likeHood = NULL;
    dirPath  = fileName;
    OutDir = "./";
    treeName = "";
    likelihood_uncond = 0.;
    likelihood_cond = 0.;
    name_params.clear();


}


void ToyFitterExclusion::fit(double mu, TString nameTree){
    
    // open the file and loop over all trees with given name
    TFile *f = TFile::Open(dirPath);
    if(f == NULL) Error("fit", TString::Format("file %s does not exist",dirPath.Data()));
    
    TFile f_out(OutDir + "post_fit_" + nameTree + "_mu"+ TString::Format("%2.2f",mu)+ ".root","RECREATE");

    // output tree, here intentionally all out tree will have the same name so we can hadd
    TTree *outTree = new TTree("out_" + nameTree, "output tree for a given mu, hadd me");

    // setting up branches on outTree
    double mu_fit    = mu;
    double testStat = 0.;
    int numberOfParams =  likeHood->getParameters()->size();
    outTree->Branch("mu_fit", &mu_fit, "mu_fit/D");
    outTree->Branch("q_mu", &testStat, "testStat/D");
    outTree->Branch("n_params", &numberOfParams, "n_params/I");
    outTree->Branch("LL_cond", &likelihood_cond, "LL_cond/D");
    outTree->Branch("LL_uncond", &likelihood_uncond, "LL_uncond/D");
    outTree->Branch("true_params", true_params,"true_params[n_params]/D");
    outTree->Branch("measured_params", measured_params,"measured_params[n_params]/D");
    outTree->Branch("uncond_params", uncond_params,"uncond_params[n_params]/D");
    outTree->Branch("cond_params", cond_params,"cond_params[n_params]/D");
    outTree->Branch("name_params", &name_params);

    // input tree
    TTree *readTree = NULL;

    // Data handler
    dataHandler data("toyDMData");
    data.setData(DM_DATA);

    // setting random seed for measured param generation
    // TRandom likes a integer seed, the multiplication just makes sure that even if mu is small you get a new seed 
    rambo.SetSeed(((ULong_t) mu * 1000.)); 

    // looping over all tree that are prefixed with nameTree
    TIter next(f->GetListOfKeys());
    TObject *treeKey = NULL;
    while ((treeKey = next())) {
     
        if(TString(treeKey->GetName()).Contains(nameTree))
            readTree = (TTree*)f->Get(treeKey->GetName());
        else { 
            Debug("fit", TString::Format("skipping, %s",treeKey->GetName()));
            continue; 
        }

        Debug("fit", TString::Format("Matched tree name, %s", treeKey->GetName()));

        // fill true values from generated tree
        fillTrueParams(readTree);

        // set toy data for fit
        data.setDataTree(readTree);
        likeHood->setDataHandler(&data);

        // Dice the measured parameters (note that if a parameter is free the t-value exytacted here has no effect)
        measureParameters();

        testStat = computeTS(mu);
        
        outTree->Fill();

    }
    
    f_out.cd();
    outTree->Write();
    f_out.Close();

}


double ToyFitterExclusion::computeTS(double mu){
    
    // we use q_tilde from equation 16 of: https://arxiv.org/abs/1007.1727
    // this is suitable for upper limit in which the parameter of interest must be >0.
    
    double qstat = 0.;
    
    // performing unconditional fit
    double LL_denominator = likeHood->maximize(false) ;

    // saving the unconditional likelihood  for outTree
    likelihood_uncond = LL_denominator;

    // printing 
    if(getPrintLevel() < 2) 
        likeHood->printCurrentParameters();
    
    // save parameters of unconditional fit 
    saveParameters(uncond_params);
    
    // override denominator in case mu_hat < 0.
    if(likeHood->POI->getCurrentValue() < 0.){
        likeHood->POI->setCurrentValue(0.);    
        LL_denominator = likeHood->maximize(true) ;
    }

    // perform conditional fit
    likeHood->POI->setCurrentValue(mu);
    double LL_numerator = likeHood->maximize(true) ;
    
    // saving the conditional likelihood for outTree
    likelihood_cond = LL_numerator;
    
    // printing 
    if(getPrintLevel() < 2) 
    likeHood->printCurrentParameters();

    // save parameters of conditional fit 
    saveParameters(cond_params);

    qstat = 2. * ( LL_denominator - LL_numerator );  // -2*logLikelihoodRatio

    // following the q_tilde construction, this could have been put above with just return 0.
    // I just wanted to do the conditional fit anyway to return conditional parameter for study
    if( likeHood->POI->getCurrentValue() > mu )  qstat = 0. ;

    

    return qstat;
}


 
void ToyFitterExclusion::saveParameters(double *outParams){

    map <int, LKParameter*> *params = likeHood->getParameters();
    
    int itr = 0;
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
        
        LKParameter *param = ip->second;
        outParams[itr] = param->getCurrentValue();
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
    


void ToyFitterExclusion::fillTrueParams(TTree *inputTree){

    // retrive the previously saved TList of parameters (done in ToyGenerator)
    TIter iterateMe(inputTree->GetUserInfo());

    TParameter<double> *parameter = NULL;
    name_params.clear();

    int i = 0;
    // saving the parameters of the new tree
    while ((parameter = (TParameter<double>*)iterateMe())) {
        true_params[i] = parameter->GetVal();
        name_params.push_back(parameter->GetName());
        i++;
    }
}



void ToyFitterExclusion::measureParameters(){

    map <int, LKParameter*> *params = likeHood->getParameters();
    
    Info("measureParameters", "Randomizing parameters:");

    int parItr = -1;
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
            
        LKParameter *param = ip->second;
        parItr++;   // start from zero 
        
        if(param->getType() == FIXED_PARAMETER || param->isOfInterest() 
            || param->getType() == FREE_PARAMETER ) {

             Info("---->","Skipping paramater: " + param->getName()); continue; 
        }
            
        // getting range
        double min = param->getMinimum();
        double max = param->getMaximum();
            
        // sample gauss tvalue for parameter
        double random_tvalue = 0.;
        random_tvalue = rambo.Gaus(0.,1.);
        // extract again if out of range
        while(random_tvalue > max || random_tvalue < min)
                    random_tvalue = rambo.Gaus(0.,1.);
    
    
        param->setT0value(random_tvalue);
        measured_params[parItr] = random_tvalue;
        Info("---->",TString::Format("new T0-Value: %s = %1.2f",param->getName().Data(),random_tvalue)); 
    }

}