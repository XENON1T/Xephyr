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


void ToyGenerator::generateCalibration(int N){

}


void ToyGenerator::generateData(double mu, int N){

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
