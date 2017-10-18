#include "ToyGenerator.h"

ToyGenerator::ToyGenerator(TString sampleName, TString outDir):errorHandler("ToyGenerator"){

    treeName = sampleName;
    dir = outDir;
    averageCalEvnt  = -9;
    averageDataEvnt = -9;
    likeHood = NULL;

}


void ToyGenerator::setSeed(int seed){

    randomizeMyass.SetSeed(seed);
}


void ToyGenerator::saveParameters(TTree *tree){
    
    TList *config = tree->GetUserInfo();

    map <int, LKParameter*> *params = likeHood->getParameters();
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
            
        LKParameter *param = ip->second;
        TParameter<double> *temp_p = new TParameter<double>(param->getName(), param->getCurrentValue());

        config->Add(temp_p);
    }
    
}

void ToyGenerator::generateCalibration(int N){

}


void ToyGenerator::generateData(double mu, int N){

    // the idea is to generate N toys with the current seed
    // Each dataset will contain as average averageDataEvnt events
    // which will be Poisson random distributed.

    // set the parameter of interest to the specified value
    likeHood->getParameter(PAR_SIGMA)->setCurrentValue(mu);

    TFile f(dir+treeName+".root","RECREATE");

    
    
    TTree toyTree (treeName, "generated toy data");

    saveParameters(&toyTree);
       
    for(int evnt =0; evnt < N ; evnt++){

    }

    toyTree.Write();
    f.Close();

}

void ToyGenerator::randomizeNuissanceParameter(){
    map <int, LKParameter*> *params = likeHood->getParameters();

    Info("randomizeNuissanceParameter", "Randomizing parameters:");
    for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
        
        LKParameter *param = ip->second;
        if(param->getType() == FIXED_PARAMETER || param->isOfInterest() ) 
            { Info("","Skipping paramater: " + param->getName()); continue; }
        
        // getting range
        double min = param->getMinimum();
        double max = param->getMaximum();
        
        // if param is Free then sample uniform otherwise gauss
        double random_tvalue = 0.;
        if(param->getType() == FREE_PARAMETER ){
            Warning("","Following parameter is FREE and will be extracted from uniform distro. " + param->getName());
            random_tvalue = randomizeMyass.Uniform(min,max);
        }
        else{ 
            random_tvalue = randomizeMyass.Gaus(0.,1.);
            // extract again if out of range
            while(random_tvalue > max || random_tvalue < min)
                random_tvalue = randomizeMyass.Gaus(0.,1.);
        }

        param->setCurrentValue(random_tvalue);

    }


    likeHood->printCurrentParameters();


}
