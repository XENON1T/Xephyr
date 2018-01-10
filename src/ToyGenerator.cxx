#include "ToyGenerator.h"

ToyGenerator::ToyGenerator(TString sampleName, TString outDir):errorHandler("ToyGenerator"){

    treeName = sampleName;
    dir = outDir;
    averageCalEvnt  = -9;
    averageDataEvnt = -9;
    likeHood = NULL;

}


void ToyGenerator::setSeed(int seed){

    rambo.SetSeed(seed);
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

void ToyGenerator::generateCalibration(int N, bool randomizeNP ){

    // generate calibration like toys, where averageCalEvnt is the poisson median
    // of the number of generated event per toy. N is the number of toy datasets.
    // The pdfs are taken only from pdfcomponent that are "safeguarded" in the 
    // likelihood, their relative rates are respected. Also the additional 
    // contribution from "AdditionalSafeGuardComponent" is included.


    if(!(averageCalEvnt>0)) Error("generateCalibration", "you MUST set averageCalEvnt.");

    // randomize initial 'true' values of NP
    if(randomizeNP)  randomizeNuissanceParameter();
        
    // rescaling to defined Calibration events
    double default_evnt = getModelIntegralSafeguarded();
    if(!(default_evnt> 0.)) Error("generateCalibration", "you MUST set Safeguarded components.");
    
    Debug("generateCalibration: default_event =", TString::Itoa(default_evnt,10));
    
    double scaleFactor  =  averageCalEvnt / default_evnt ;

    Debug("generateCalibration: scaleFactor =", TString::Itoa(scaleFactor,10));
    
    // retrive a vector of TH2F of interpolated bkg components
    vector <TH2F> backgrounds = getTH2OfBkg();
    
    TFile f(dir+treeName+"_Cal.root","RECREATE");
        
    // necessary because TH2F::GetRandom uses ROOT::gRandom
    gRandom = &rambo;
    
    // actual generation of N toys with poisson fluctuating events.
    for(int toyItr =0; toyItr < N ; toyItr++){
        
        TString name = treeName + "_Cal_" + TString::Itoa(toyItr,10); 
        TTree toyTree (name, "generated toy Calibration");
        float cs1 = 0.; 
        float cs2 = 0.;
        string type = "DummyLabel";
        
        toyTree.Branch("cs1",&cs1,"cs1/F");
        toyTree.Branch("cs2",&cs2,"cs2/F");
        toyTree.Branch("type",&type);
        
        saveParameters(&toyTree);
                       
        // loop over each component extract N events and dice s1-s2
        for(unsigned int bkgItr=0; bkgItr < backgrounds.size(); bkgItr++){
            
            // only safeguarded conponent
            if(!(likeHood->safeguarded_bkg_components[bkgItr])) continue;
            
            Debug("generateCalibration","");

            int N_events   = rambo.Poisson(scaleFactor * backgrounds[bkgItr].Integral());
                    
            type = (likeHood->bkg_components[bkgItr])->getName();

            Debug("generateCalibration", TString::Format("Generating %d events for %s, with median %f",N_events, (likeHood->bkg_components[bkgItr])->getName().Data(), scaleFactor * backgrounds[bkgItr].Integral()));
            for(int evt =0; evt < N_events; evt++){
        
                double temp_cs1 = 0., temp_cs2 = 0.;
                backgrounds[bkgItr].GetRandom2(temp_cs1,temp_cs2);
                cs1 = (float) temp_cs1;// TODO FIXME this is just for back compatibility cs1 and cs2 must be double
                cs2 = (float) temp_cs2;
                toyTree.Fill();
            }
        }

        // adding the "additional component" events
        int N_additional   = 
            (likeHood->safeguardAdditionalComponent !=NULL) ? rambo.Poisson(scaleFactor * likeHood->safeguardAdditionalComponent->Integral()) : 0;

        type = "additional";
        Debug("generateCalibration", TString::Format("Generating %d events for additional component",N_additional));
        for(int evt =0; evt < N_additional; evt++){
            double temp_cs1 = 0., temp_cs2 = 0.; // TODO FIXME this is just for back compatibility cs1 and cs2 must be double
            likeHood->safeguardAdditionalComponent->GetRandom2(temp_cs1,temp_cs2);
            cs1 = (float) temp_cs1;// TODO FIXME this is just for back compatibility cs1 and cs2 must be double
            cs2 = (float) temp_cs2;
            toyTree.Fill();
        }

        toyTree.Write();
    }

    f.Close();

}


void ToyGenerator::generateData(double mu, int N, bool randomizeNP){

    // the idea is to generate N toys with the current seed
    // Each dataset will contain as average averageDataEvnt events
    // which will be Poisson random distributed.

    // randomize initial 'true' values of NP
    if(randomizeNP) randomizeNuissanceParameter();
    // set the parameter of interest to the specified value
    likeHood->getParameter(PAR_SIGMA)->setCurrentValue(mu);

    // rescaling to defined data events
    double default_evnt = getModelIntegral();
    double scaleFactor  = (averageDataEvnt>0.) ? averageDataEvnt / default_evnt : 1.;

    // retrive a vector of TH2F of interpolated bkg components
    vector <TH2F> backgrounds = getTH2OfBkg();

    TFile f(dir+treeName+".root","RECREATE");
    
    // necessary because TH2F::GetRandom uses ROOT::gRandom
    gRandom = &rambo;

    // actual generation of N toys with poisson fluctuating events.
    for(int toyItr =0; toyItr < N ; toyItr++){

        TString name = treeName + "_" + TString::Itoa(toyItr,10); 
        TTree toyTree (name, "generated toy data");
        float cs1 = 0.; 
        float cs2 = 0.;
        string type = "DummyLabel";

        toyTree.Branch("cs1",&cs1,"cs1/F");
        toyTree.Branch("cs2",&cs2,"cs2/F");
        toyTree.Branch("type",&type);

        saveParameters(&toyTree);
    
        
        // loop over each bkg extract N events and dice s1-s2
        for(unsigned int bkgItr=0; bkgItr < backgrounds.size(); bkgItr++){
            int N_events   = rambo.Poisson(scaleFactor * backgrounds[bkgItr].Integral());
            
            type = (likeHood->bkg_components[bkgItr])->getName();
            for(int evt =0; evt < N_events; evt++){
                double temp_cs1 = 0., temp_cs2 = 0.;
                backgrounds[bkgItr].GetRandom2(temp_cs1,temp_cs2);
                cs1 = (float) temp_cs1;// TODO FIXME this is just for back compatibility cs1 and cs2 must be double
                cs2 = (float) temp_cs2;
                toyTree.Fill();
            }
        }

        if(mu > 0) {
          // filling signal if any
          type = likeHood->signal_component->getName();
          int N_signal  = rambo.Poisson(likeHood->getCurrentNs());
          TH2F signal   = likeHood->signal_component->getInterpolatedHisto();
          
          for(int evt =0; evt < N_signal; evt++){
            double temp_cs1 = 0., temp_cs2 = 0.;
            signal.GetRandom2(temp_cs1,temp_cs2);
            cs1 = (float) temp_cs1;// TODO FIXME this is just for back compatibility cs1 and cs2 must be double
            cs2 = (float) temp_cs2;
            toyTree.Fill();
          }
        
        }

        toyTree.Write();

    }

    likeHood->printCurrentParameters();
    
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
            random_tvalue = rambo.Uniform(min,max);
        }
        else{ 
            random_tvalue = rambo.Gaus(0.,1.);
            // extract again if out of range
            while(random_tvalue > max || random_tvalue < min)
                random_tvalue = rambo.Gaus(0.,1.);
        }

        param->setCurrentValue(random_tvalue);

    }


    likeHood->printCurrentParameters();


}


double ToyGenerator::getModelIntegral(){
    
    double total_integral = 0.;
    unsigned int n = likeHood->bkg_components.size();

    for (unsigned int i=0; i < n; i++ ){

        total_integral += (likeHood->bkg_components[i])->getDefaultEvents();
    }

    return total_integral;
}


double ToyGenerator::getModelIntegralSafeguarded(){
    
    double total_integral = 0.;
    unsigned int n = likeHood->bkg_components.size();

    Info("getModelIntegralSafeguarded","Components used for calibration generation:");

    for (unsigned int i=0; i < n; i++ ){

        if(!(likeHood->safeguarded_bkg_components[i]))  continue;
        Info("---->",(likeHood->bkg_components[i])->getName() );
        total_integral += (likeHood->bkg_components[i])->getDefaultEvents();
    }

    // adding the integral of additional component
    total_integral += (likeHood->safeguardAdditionalComponent != NULL) ? likeHood->safeguardAdditionalComponent->Integral() : 0.;

    return total_integral;
}



vector<TH2F> ToyGenerator::getTH2OfBkg(){

    vector<TH2F> temp_v;
    unsigned int n = likeHood->bkg_components.size();
    
        for (unsigned int i=0; i < n; i++ ){
    
            temp_v.push_back((likeHood->bkg_components[i])->getInterpolatedHisto());
        }

    return temp_v;

}