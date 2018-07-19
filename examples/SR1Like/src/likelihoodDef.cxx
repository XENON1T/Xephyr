#include "XeLikelihoods.h"
#include "XePdfObjects.h"
#include "dataHandler.h"
#include "TSystem.h"
#include "signalDef.h"
#include <vector>

///---------------     PARAMETERS ------------         ///

//============ SR1 SPECIFIC ============//
const double livedays_SR1 = 246.74 ;  // days
const double xe_total_mass_SR1  = 1.363 ;  // tonne --- THIS IS TOTAL MASS TO WHICH HISTO ARE EXTRUDED

const double xe_inner_mass_SR1_old  = 1.0064;  // tonne    USED FOR SCALING ER
const double xe_egg_mass_SR1  =   0.646157 ;  // tonne    
const double xe_inner_mass_SR1  = 0.360478 ;  // tonne    
const double xe_outer_mass_SR1  = 0.238908;  // tonne

const double inner_r_SR1  = 36.4822 ; // cm
const double outer_r_SR1  =  41.2396 ;   // cm
const double year         = 365. ;   // days
const double exposure_SR1 = livedays_SR1 * xe_total_mass_SR1 / year ; // tonne * year

const double er_events_inner_volume_SR1_old = 408 ;
const double er_events_outer_volume_SR1 =  er_events_inner_volume_SR1_old * xe_outer_mass_SR1 / xe_inner_mass_SR1_old ;
const double er_events_inner_volume_SR1 =  er_events_inner_volume_SR1_old * xe_inner_mass_SR1 / xe_inner_mass_SR1_old ;
const double er_events_egg_volume_SR1 =    er_events_inner_volume_SR1_old * xe_egg_mass_SR1 / xe_inner_mass_SR1_old ;

const double calibration_events_egg_SR1   = 3394. ; // events
const double calibration_events_inner_SR1 = 1751. ; // events
const double calibration_events_outer_SR1 = 1203. ; // events


//============== GLOBAL for TYPE LIKELIHOOD ==================//
          ///    {  SR1_EGG  ;  SR1_INNER ;   SR1_OUTER  }
const int n_types = 3;
const double xe_mass[n_types]  = { xe_egg_mass_SR1, xe_inner_mass_SR1, xe_outer_mass_SR1 } ;  // tonne  
const double exposure[n_types] = { exposure_SR1 , exposure_SR1 , exposure_SR1 } ;

const double exposure_signal[n_types] = { livedays_SR1 *xe_egg_mass_SR1  / year , 
                                          livedays_SR1 * xe_inner_mass_SR1 / year , 
                                          livedays_SR1 * xe_outer_mass_SR1 / year } ;
                                      
const double livedays[n_types] = { livedays_SR1 , livedays_SR1 , livedays_SR1 } ;

const double er_exp_events[n_types] = { er_events_egg_volume_SR1 , er_events_inner_volume_SR1 , er_events_outer_volume_SR1 };
const double calibration_events[n_types] = {calibration_events_egg_SR1,  calibration_events_inner_SR1, calibration_events_outer_SR1  } ;

  // LAX Version
const TString LAX = "lax_1.5.1_egg";

  // INPUT FILES
const TString ER_file =  "ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180420.root" ;

const TString AC_file =  "Background_wall_ac_templates_v6_inegg_bin_SR1_2018-04-30.root";

const TString Radio_file = "RadiogenicNRBackgroundModel_SR1_RunSourceCombined180420.root";

const TString Radio_NX_file =  "RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombined180420.root";

const TString CNNS_file = "CNNSBackgroundModel_SR1_RunSourceCombinedFit180420.root";



// HISTO NAMES: each volume region has its own set of histograms
// The ER is special, the names are actually suffixes
const TString ER_name[n_types] =       { "_yx_r0_f" ,  "_yx_r1_f" , "_yx_r2_f" };
const TString AC_name[n_types] =       { "acbg_yx_r0_f"      ,  "acbg_yx_r1_f", "acbg_yx_r2_f" };
const TString AC_rn_name[n_types] =    { "acrn220_yx_r0_f", "acrn220_yx_r1_f", "acrn220_yx_r2_f" };
const TString Wall_name[n_types] =     { "wallbg_0.00_yx_r0_f" , "wallbg_0.00_yx_r1_f" , "wallbg_0.00_yx_r2_f" };
const TString Wall_rn_name[n_types] =  { "wallrn_0.00_yx_r0_f" , "wallrn_0.00_yx_r1_f" , "wallrn_0.00_yx_r2_f" };
const TString Radio_name[n_types] =    { "hnrbkg_yx_r0_f", "hnrbkg_yx_r1_f" , "hnrbkg_yx_r2_f" };
const TString Radio_NX_name[n_types] = { "hnrbkg_yx_r0_f", "hnrbkg_yx_r1_f" , "hnrbkg_yx_r2_f" };
const TString CNNS_name[n_types] =     { "hmc_extruded_yx_r0_f", "hmc_extruded_yx_r1_f" , "hmc_extruded_yx_r2_f" };


// Uncertainties, in case the uncertainty would depend on the volume region (egg, U, wall)  
const double wall_unc[n_types] =  { 1. , 1., 0.5 };

/// --------- ///


// This function returns a skeleton likelihood definition, without data.
pdfLikelihood* getTheLikelihood_SR1( double mass, unsigned int type, double other = 0.){

 TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
 TString inputDir = xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/"+ LAX + "/sliced_2d_templates/";


 if(type != 1 && type != 0 && type !=2 )  exit(100);

 TString type_str = TString::Itoa(type, 10);
 
 pdfComponent *ER = new pdfComponent("hbkg", inputDir + ER_file );

    shapeSys *PY = new shapeSys("_py0_");     
    PY->setStep(1.);                     
    PY->setMinimum(-1.);   
    PY->setMaximum(1.);    
    PY->setType(FREE_PARAMETER);
    shapeSys *RF = new shapeSys("_rf0_");
    RF->setStep(1.);                     
    RF->setMinimum(-2.);   
    RF->setMaximum(2.);    
    RF->setType(FREE_PARAMETER);
    scaleSys *s1 = new scaleSys("ERscale"+type_str,1.);
    s1->setMinimum(-0.5);
    s1->setMaximum(0.5);
    s1->setType(FREE_PARAMETER);
    ER->addScaleSys(s1);   

    ER->addShapeSys(PY);   
    ER->addShapeSys(RF);

    ER->suffix = ER_name[type];
    ER->setEvents( er_exp_events[type] ); // Set the normalization factor for ER in events.


  
  pdfComponent *AC = new pdfComponent(AC_name[type], inputDir + AC_file );

    scaleSys *sAC = new scaleSys("ACscale"+type_str,0.2);
    sAC->setMinimum(-2.);
    sAC->setMaximum(2.);
    //sAC->setType(FIXED_PARAMETER);
    AC->addScaleSys(sAC);


    TFile *f1 = TFile::Open(inputDir + AC_file );
    TH2F* AC_for_safeguard = (TH2F*) f1->Get( AC_rn_name[type] );
    TH2F* wall_for_safeguard = (TH2F*) f1->Get( Wall_rn_name[type] );

    double AC_scale_safeguard = er_exp_events[type]  / calibration_events[type];
    
    AC_for_safeguard->Add(wall_for_safeguard);
    AC_for_safeguard->Scale( AC_scale_safeguard );  // must be scaled to ER

  
    pdfComponent *Wall = new pdfComponent( Wall_name[type] , inputDir + AC_file );

    scaleSys *sWall = new scaleSys("Wallscale"+type_str, wall_unc[type] );
    sWall->setMinimum( -1. / wall_unc[type] + 0.001);
    sWall->setMaximum(3.);
    //sWall->setType(FIXED_PARAMETER);
    Wall->addScaleSys(sWall);

  
    pdfComponent *Radio = new pdfComponent( Radio_name[type] , inputDir + Radio_file );

    scaleSys *sRadio = new scaleSys("Radioscale"+type_str,0.5 );
    sRadio->setMinimum(-2.);
    sRadio->setMaximum(2.);
    //sRadio->setType(FIXED_PARAMETER);
    Radio->addScaleSys(sRadio);
    
    Radio->setScaleFactor( livedays[type] / year );


    pdfComponent *RadioNX = new pdfComponent( Radio_NX_name[type] , inputDir + Radio_NX_file );

    scaleSys *sRadioNX = new scaleSys("RadioscaleNX"+type_str,0.5 );
    sRadioNX->setMinimum(-2.);
    sRadioNX->setMaximum(2.);
    //sRadio->setType(FIXED_PARAMETER);
    RadioNX->addScaleSys(sRadioNX);
    
    RadioNX->setScaleFactor( livedays[type] / year );


  
    pdfComponent *CNNS = new pdfComponent( CNNS_name[type] , inputDir + CNNS_file );

    CNNS->setScaleFactor( exposure[type] );  // has been cutted in R after extrusion.
    scaleSys *sCNNS = new scaleSys("CNNSscale"+type_str,0.3);
    //sCNNS->setType(FIXED_PARAMETER);
    sCNNS->setMinimum(-2.);
    sCNNS->setMaximum(2.);
    CNNS->addScaleSys(sCNNS);


    TString Signal_name = getSignalHistoName(mass, type, other);
    TString signal_file_with_mass = getSignalFile(mass, type, other);
    double  signal_unc = getSignalUncertainty(mass, type, other);

    pdfComponent *Signal = new pdfComponent( Signal_name , signal_file_with_mass );
    Signal->setScaleFactor( exposure_signal[type] );

    scaleSys *s2 = new scaleSys("SignalScale"+type_str, signal_unc) ;
    Signal->addScaleSys(s2);
    

  pdfLikelihood *pl = new pdfLikelihood("L_SR1_V"+type_str, mass);
    
    pl->setExperiment(1);
    pl->addBkgPdfComponent(CNNS, false);
    pl->addBkgPdfComponent(RadioNX, false);
    pl->addBkgPdfComponent(Radio, false);
    pl->addBkgPdfComponent(AC, false);
    pl->addBkgPdfComponent(Wall, false);
    pl->addBkgPdfComponent(ER, true);
    pl->setSignalPdf(Signal);       
    pl->setSignalDefaultNorm(1.E-45);
    //pl.setPrintLevel(WARNING);
    pl->setWithSafeGuard(false);
    //pl->setFixedValueForSafeguard(0.00001);
    pl->setAdditionalSafeGuardComponent(AC_for_safeguard);
      
    return pl;

}



dataHandler* getDataToFit( TString collection, int mu, double mass, int type, int Gen ){ 
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
  TString inputDir = xeDir + "build/RESULTS/GENtrees/";
  
  TString data_filename = collection + TString::Format("_M%1.0f_mu%d_G%d_V%d.root",mass, mu, Gen, type); 
  TString data_treeName = collection +  TString::Format("_M%1.0f_mu%d_G%d_V%d",mass, mu, Gen, type) ;
   
  dataHandler *data = new dataHandler("dmData",inputDir + data_filename, data_treeName+ "_0");
  data->setPrefixTree(data_treeName);
  data->setPrintLevel(DEBUG);

  return data;
}
   

dataHandler* getCalibToFit( TString collection,  double mass, int type, int Gen ){ 
  
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
  TString inputDir = xeDir + "build/RESULTS/GENtrees/";
  
  TString data_filename = collection + TString::Format("_M%1.0f_mu0_G%d_V%d_Cal.root", mass, Gen, type); 
  TString data_treeName = collection +  TString::Format("_M%1.0f_mu0_G%d_V%d_Cal", mass, Gen, type) ;
   
  dataHandler *data = new dataHandler("calibration",inputDir + data_filename, data_treeName+ "_0");
  data->setPrefixTree(data_treeName);
  data->setPrintLevel(DEBUG);

  return data;

}


pdfLikelihood* getTheLikelihoodToFit(TString collection, int mu, double mass, double other, int type, int Gen){
  
  dataHandler *data  = getDataToFit(collection, mu, mass, type, Gen ) ;
  dataHandler *calib = getCalibToFit(collection,mass, type, Gen );

  pdfLikelihood *pl = getTheLikelihood_SR1(mass, type, other);
  pl->setDataHandler(data);
  pl->setCalibrationData(calib);
  
  pl->setWithSafeGuard(true);
  pl->initialize();
  
  return pl;
}

dataHandler* getDMdata( int type){ 
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
  TString inputDir = xeDir + "xephyr_examples/SR1Like/data/";

  TString data_filename = "xephyr_none_SR1_pax6.8.0_hax2.4.0_lax1.5.1_cs1LT200_fv1_cuts1_pruned.root" ;
  TString data_treeName = "tree_"+TString::Itoa(type,10);

  dataHandler *data = new dataHandler("data_sr1_"+TString::Itoa(type,10),inputDir + data_filename, data_treeName);

  return data;

}

dataHandler* getDMCalibration( int type){ 
  
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
  TString inputDir = xeDir + "xephyr_examples/SR1Like/data/";

  TString data_filename = "xephyr_Rn220_SR1_pax6.8.0_hax2.4.0_lax1.5.1_cs1LT200_fv1_cuts1_pruned.root" ;
  TString data_treeName = "tree_"+TString::Itoa(type,10);

  dataHandler *data = new dataHandler("calibration_sr1_"+TString::Itoa(type,10),inputDir + data_filename, data_treeName);

  return data;

}

pdfLikelihood* getDMLikelihood( double mass, int type, double other = 0.){
  
  dataHandler *data  = getDMdata(type) ;
  dataHandler *calib = getDMCalibration(type);

  pdfLikelihood *pl = getTheLikelihood_SR1(mass, type, other );
  pl->setDataHandler(data);
  pl->setCalibrationData(calib);
  pl->setWithSafeGuard(true);
  pl->initialize();
  
  return pl;
}


CombinedProfileLikelihood* combine(vector <pdfLikelihood*> p) {
  
  CombinedParameter *PY_sr1 = new CombinedParameter("CORR_PY_SR1");
  CombinedParameter *RF_sr1 = new CombinedParameter("CORR_RF_SR1");
  CombinedParameter *ER_sr1 = new CombinedParameter("CORR_ER_SR1");
  CombinedParameter *cnns_sr1 = new CombinedParameter("CORR_CNNS_SR1");
  CombinedParameter *signal_sr1 = new CombinedParameter("CORR_SIGNALScale_SR1");

  // loop over likelihood to combine and set common parameters
  for(unsigned int i=0 ; i < p.size(); i++) {

    TString type = TString::Itoa(i, 10);

    p[i]->setExperiment(i + 1);
    p[i]->setPrintLevel(INFO);
    
    PY_sr1->correlateParameter( p[i]->getBkgComponent("hbkg")->getShapeSys("_py0_") );
    RF_sr1->correlateParameter( p[i]->getBkgComponent("hbkg")->getShapeSys("_rf0_") );
    ER_sr1->correlateParameter( p[i]->getBkgComponent("hbkg")->getScaleSys("ERscale"+type) );
    cnns_sr1->correlateParameter( p[i]->getBkgComponent("hmc_extruded_yx_r"+type+"_f")->getScaleSys("CNNSscale"+type) );
    signal_sr1->correlateParameter( p[i]->signal_component->getScaleSys("SignalScale"+type) );

  }

  // set limits of parameters 
  PY_sr1->setType(FREE_PARAMETER); 
  PY_sr1->setMinimum(-1.); PY_sr1->setMaximum(1.);

  RF_sr1->setType(FREE_PARAMETER); 
  RF_sr1->setMinimum(-2.); RF_sr1->setMaximum(2.);

  ER_sr1->setType(FREE_PARAMETER);
  ER_sr1->setMinimum(-0.5); ER_sr1->setMaximum(0.5);

  cnns_sr1->setMinimum(-2); cnns_sr1->setMaximum(2);
  

  signal_sr1->setMinimum(-2); signal_sr1->setMaximum(2);


  CombinedProfileLikelihood *cpl = new CombinedProfileLikelihood("SR1_combo");
  cpl->addParameter(RF_sr1, AUTO);
  cpl->addParameter(PY_sr1, AUTO);
  cpl->addParameter(ER_sr1,AUTO);
  cpl->addParameter(cnns_sr1,AUTO);
  cpl->addParameter(signal_sr1,AUTO);

  // combine! 
  for(unsigned int i=0 ; i < p.size(); i++) cpl->combine(p[i]);

  cpl->initialize();
  cpl->setPrintLevel(INFO);

  return cpl;
  
}

CombinedProfileLikelihood* getTheCombinedLikelihood(TString collection, int mu, double mass, int Gen, double other = 0.) {

  vector <pdfLikelihood*> like;
  like.push_back( getTheLikelihoodToFit(collection, mu, mass, other, 0, Gen)  );
  like.push_back( getTheLikelihoodToFit(collection, mu, mass, other, 1, Gen)  );
  like.push_back( getTheLikelihoodToFit(collection, mu, mass, other, 2, Gen)  );

  CombinedProfileLikelihood* cpl = combine( like );

  return cpl;

}


CombinedProfileLikelihood* getDMCombinedLikelihood( double mass , double other = 0.) {

  vector <pdfLikelihood*> like;
  like.push_back( getDMLikelihood(mass, 0, other) );
  like.push_back( getDMLikelihood(mass, 1, other) );
  like.push_back( getDMLikelihood(mass, 2, other) );

  CombinedProfileLikelihood* cpl = combine( like );

  return cpl;
}


