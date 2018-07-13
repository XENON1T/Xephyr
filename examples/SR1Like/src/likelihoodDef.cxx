#include "XeLikelihoods.h"
#include "XePdfObjects.h"
#include "dataHandler.h"
#include "TSystem.h"
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
const double livedays[n_types] = { livedays_SR1 , livedays_SR1 , livedays_SR1 } ;

const double er_exp_events[n_types] = { er_events_egg_volume_SR1 , er_events_inner_volume_SR1 , er_events_outer_volume_SR1 };
const double calibration_events[n_types] = {calibration_events_egg_SR1,  calibration_events_inner_SR1, calibration_events_outer_SR1  } ;

  // LAX Version
const TString LAX = "lax_1.5.1_egg";

  // INPUT FILES
const TString ER_file[n_types] = { "ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180420.root" , 
                             "ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180420.root",
                             "ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180420.root"  };

//const TString ER_file[n_types] = { "v1_ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180417.root" , 
//                             "v1_ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180417.root",
//                             "v1_ERBackgroundModel_Stitched_SR1_RunSourceCombinedFit180417.root"  };


const TString AC_file[n_types] = { "Background_wall_ac_templates_v6_inegg_bin_SR1_2018-04-30.root", 
                                   "Background_wall_ac_templates_v6_inegg_bin_SR1_2018-04-30.root",
                                   "Background_wall_ac_templates_v6_inegg_bin_SR1_2018-04-30.root" };

const TString Radio_file[n_types] = { "RadiogenicNRBackgroundModel_SR1_RunSourceCombined180420.root",
                                      "RadiogenicNRBackgroundModel_SR1_RunSourceCombined180420.root",
                                      "RadiogenicNRBackgroundModel_SR1_RunSourceCombined180420.root"};

//const TString Radio_file[n_types] = { "v1_RadiogenicNRBackgroundModel_SR1_RunSourceCombinedFit180417.root",
//                                      "v1_RadiogenicNRBackgroundModel_SR1_RunSourceCombinedFit180417.root",
//                                      "v1_RadiogenicNRBackgroundModel_SR1_RunSourceCombinedFit180417.root"};

//const TString Radio_NX_file[n_types] = { "v1_RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombinedFit180417.root",
//                                         "v1_RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombinedFit180417.root",
//                                         "v1_RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombinedFit180417.root"};

const TString Radio_NX_file[n_types] = { "RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombined180420.root",
                                        "RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombined180420.root",
                                         "RadiogenicNeutronXBackgroundModel_SR1_RunSourceCombined180420.root"};

const TString CNNS_file[n_types] = { "CNNSBackgroundModel_SR1_RunSourceCombinedFit180420.root",
                                     "CNNSBackgroundModel_SR1_RunSourceCombinedFit180420.root",  
                                     "CNNSBackgroundModel_SR1_RunSourceCombinedFit180420.root" };

const TString Signal_file[n_types] = { "SignalModel_SR1_GeV_FromRunSourceCombinedFit180420.root", 
                                       "SignalModel_SR1_GeV_FromRunSourceCombinedFit180420.root", 
                                       "SignalModel_SR1_GeV_FromRunSourceCombinedFit180420.root" } ;
                                 
//const TString Signal_file[n_types] = { "v1_SignalModel_SR1_GeV_FromRunSourceCombinedFit180417.root", 
//                                       "v1_SignalModel_SR1_GeV_FromRunSourceCombinedFit180417.root", 
//                                       "v1_SignalModel_SR1_GeV_FromRunSourceCombinedFit180417.root" } ;

  // HISTO NAMES
const TString ER_name[n_types] =       { "_yx_r0_f" ,  "_yx_r1_f" , "_yx_r2_f" };
const TString AC_name[n_types] =       { "acbg_yx_r0_f"      ,  "acbg_yx_r1_f", "acbg_yx_r2_f" };
const TString AC_rn_name[n_types] =    { "acrn220_yx_r0_f", "acrn220_yx_r1_f", "acrn220_yx_r2_f" };
const TString Wall_name[n_types] =     { "wallbg_0.00_yx_r0_f" , "wallbg_0.00_yx_r1_f" , "wallbg_0.00_yx_r2_f" };
const TString Wall_rn_name[n_types] =  { "wallrn_0.00_yx_r0_f" , "wallrn_0.00_yx_r1_f" , "wallrn_0.00_yx_r2_f" };
const TString Radio_name[n_types] =    { "hnrbkg_yx_r0_f", "hnrbkg_yx_r1_f" , "hnrbkg_yx_r2_f" };
const TString Radio_NX_name[n_types] = { "hnrbkg_yx_r0_f", "hnrbkg_yx_r1_f" , "hnrbkg_yx_r2_f" };
const TString CNNS_name[n_types] =     { "hmc_extruded_yx_r0_f", "hmc_extruded_yx_r1_f" , "hmc_extruded_yx_r2_f" };
const TString Signal_name[n_types] =   { "hmc_extruded_yx_r0_f", "hmc_extruded_yx_r1_f" , "hmc_extruded_yx_r2_f" } ;


  // Uncertainties 
const double wall_unc[n_types] =  { 1. , 1., 0.5 };


/// --------- ///


pdfLikelihood* getTheLikelihood_SR1( double mass, unsigned int type){

 TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
 TString inputDir = xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/"+ LAX + "/sliced_2d_templates/";

 if(type != 1 && type != 0 && type !=2 )  exit(100);

 TString type_str = TString::Itoa(type, 10);
 
 pdfComponent *ER = new pdfComponent("hbkg", inputDir + ER_file[type] );

    //Now some settings:
    /*shapeSys *Eth = new shapeSys("_Eth_");     
    Eth->setStep(1.);                      // FIXME
    Eth->setMinimum(-2.);   
    Eth->setMaximum(2.);    
    Eth->setType(FIXED_PARAMETER);
    */

    shapeSys *PY = new shapeSys("_py0_");     
    PY->setStep(1.);                      // FIXME
    PY->setMinimum(-1.);   
    PY->setMaximum(1.);    
    PY->setType(FREE_PARAMETER);
    shapeSys *RF = new shapeSys("_rf0_");
    RF->setStep(1.);                     // FIXME
    RF->setMinimum(-2.);   
    RF->setMaximum(2.);    
    RF->setType(FREE_PARAMETER);
    scaleSys *s1 = new scaleSys("ERscale"+type_str,1.);
    s1->setMinimum(-0.5);
    s1->setMaximum(0.5);
    s1->setType(FREE_PARAMETER);
    ER->addScaleSys(s1);   

    //ER->addShapeSys(Eth);   
    ER->addShapeSys(PY);   
    ER->addShapeSys(RF);

    ER->suffix = ER_name[type];
    ER->setEvents( er_exp_events[type] ); // Set the normalization factor for ER in events.


  
  pdfComponent *AC = new pdfComponent(AC_name[type], inputDir + AC_file[type] );

    scaleSys *sAC = new scaleSys("ACscale"+type_str,0.2);
    sAC->setMinimum(-2.);
    sAC->setMaximum(2.);
    //sAC->setType(FIXED_PARAMETER);
    AC->addScaleSys(sAC);


    TFile *f1 = TFile::Open(inputDir + AC_file[type] );
    TH2F* AC_for_safeguard = (TH2F*) f1->Get( AC_rn_name[type] );
    TH2F* wall_for_safeguard = (TH2F*) f1->Get( Wall_rn_name[type] );

    double AC_scale_safeguard = er_exp_events[type]  / calibration_events[type];
    
    AC_for_safeguard->Add(wall_for_safeguard);
    AC_for_safeguard->Scale( AC_scale_safeguard );  // must be scaled to ER

  
    pdfComponent *Wall = new pdfComponent( Wall_name[type] , inputDir + AC_file[type] );

    scaleSys *sWall = new scaleSys("Wallscale"+type_str, wall_unc[type] );
    sWall->setMinimum( -1. / wall_unc[type] + 0.001);
    sWall->setMaximum(3.);
    //sWall->setType(FIXED_PARAMETER);
    Wall->addScaleSys(sWall);

  
    pdfComponent *Radio = new pdfComponent( Radio_name[type] , inputDir + Radio_file[type] );

    scaleSys *sRadio = new scaleSys("Radioscale"+type_str,0.5 );
    sRadio->setMinimum(-2.);
    sRadio->setMaximum(2.);
    //sRadio->setType(FIXED_PARAMETER);
    Radio->addScaleSys(sRadio);
    
    Radio->setScaleFactor( livedays[type] / year );


    pdfComponent *RadioNX = new pdfComponent( Radio_NX_name[type] , inputDir + Radio_NX_file[type] );

    scaleSys *sRadioNX = new scaleSys("RadioscaleNX"+type_str,0.5 );
    sRadioNX->setMinimum(-2.);
    sRadioNX->setMaximum(2.);
    //sRadio->setType(FIXED_PARAMETER);
    RadioNX->addScaleSys(sRadioNX);
    
    RadioNX->setScaleFactor( livedays[type] / year );


  
    pdfComponent *CNNS = new pdfComponent( CNNS_name[type] , inputDir + CNNS_file[type] );

    CNNS->setScaleFactor( exposure[type] );  // has been cutted in R after extrusion.
    scaleSys *sCNNS = new scaleSys("CNNSscale"+type_str,0.3);
    //sCNNS->setType(FIXED_PARAMETER);
    sCNNS->setMinimum(-2.);
    sCNNS->setMaximum(2.);
    CNNS->addScaleSys(sCNNS);


    
    //TString Signal_file = "sliced_2d_templates/InelasticWIMPSignalModel_1000GeV_Delta100GeV_FromCombinedFit180212.root" ;
    //TString Signal_file = "sliced_2d_templates/InelasticWIMPSignalModel_50GeV_Delta100keV_FromCombinedFit180212.root" ;
    TString signal_file_with_mass = Signal_file[type];
    TString m = TString::Format("%1.0f", mass);
    signal_file_with_mass.Insert( signal_file_with_mass.Index("GeV"), m );  // adding the mass to the file name


    pdfComponent *Signal = new pdfComponent( Signal_name[type] , inputDir + signal_file_with_mass );
    Signal->setScaleFactor( exposure[type] );

    TFile f_uncertainty(xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/"+LAX+"/models_extended/Acceptance_RunSourceCombined180420_SR1_WIMP.root");
    //TFile f_uncertainty(xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/"+LAX+"/models_extended/v1_Acceptance_RunSourceCombined180417_SR1_WIMP.root");
    TGraph *signal_unc = (TGraph*)f_uncertainty.Get("Graph");
    cout << "scale Signal uncertainty: " << signal_unc->Eval(mass) << endl;
    
    scaleSys *s2 = new scaleSys("SignalScale"+type_str, signal_unc->Eval(mass));
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
    //pl->setSignalDefaultNorm(2.9627e-38);
    //pl.setPrintLevel(WARNING);
    pl->setWithSafeGuard(false);
    //pl->setFixedValueForSafeguard(0.00001);
    pl->setAdditionalSafeGuardComponent(AC_for_safeguard);
      
    return pl;

}



pdfLikelihood* getTheLikelihoodByType(double mass, unsigned int type){
  
  TString tipo = TString::Itoa(type, 10);

  if(type == 0 || type == 1 || type == 2) {
    cout << "Likelihood Type: SR1 Volume "+ tipo + ", mass  "<< mass << endl;
    return getTheLikelihood_SR1(mass, type);
  }
  else cout << "ERROR no such a Likelihood Type " << endl;
  
  return new pdfLikelihood("DUMMY",999);

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


pdfLikelihood* getTheLikelihoodToFit(TString collection, int mu, double mass, int type, int Gen){
  
  dataHandler *data  = getDataToFit(collection, mu, mass, type, Gen ) ;
  dataHandler *calib = getCalibToFit(collection,mass, type, Gen );

  pdfLikelihood *pl =  getTheLikelihoodByType(mass, type);
  pl->setDataHandler(data);
  pl->setCalibrationData(calib);
  
  // bool add_safeguard = ( type ==0 ) ;  // not really sure I wanna do this THINK
  pl->setWithSafeGuard(true);
  pl->initialize();
  
  return pl;
}

dataHandler* getDMdata( int type){ 
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
  TString inputDir = xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/" + LAX + "/data/";

  //TString data_filename = "xephyr_none_SR1_pax6.8.0_hax2.4.0_lax1.5.1_cs1LT200_fv1_cuts1.root" ;
  TString data_filename = "xephyr_none_SR1_pax6.8.0_hax2.4.0_lax1.5.1_cs1LT200_fv1_cuts1_pruned.root" ;
  TString data_treeName = "tree_"+TString::Itoa(type,10);

  dataHandler *data = new dataHandler("data_sr1_"+TString::Itoa(type,10),inputDir + data_filename, data_treeName);

  return data;

}

dataHandler* getDMCalibration( int type){ 
  
  TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
  TString inputDir = xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/" + LAX +"/data/";

  TString data_filename = "xephyr_Rn220_SR1_pax6.8.0_hax2.4.0_lax1.5.1_cs1LT200_fv1_cuts1_pruned.root" ;
  TString data_treeName = "tree_"+TString::Itoa(type,10);

  dataHandler *data = new dataHandler("calibration_sr1_"+TString::Itoa(type,10),inputDir + data_filename, data_treeName);

  return data;

}

pdfLikelihood* getDMLikelihood( double mass, int type){
  
  dataHandler *data  = getDMdata(type) ;
  dataHandler *calib = getDMCalibration(type);

  pdfLikelihood *pl =  getTheLikelihoodByType(mass, type);
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

CombinedProfileLikelihood* getTheCombinedLikelihood(TString collection, int mu, double mass, int Gen) {

  vector <pdfLikelihood*> like;
  like.push_back( getTheLikelihoodToFit(collection, mu, mass, 0, Gen)  );
  like.push_back( getTheLikelihoodToFit(collection, mu, mass, 1, Gen)  );
  like.push_back( getTheLikelihoodToFit(collection, mu, mass, 2, Gen)  );

  CombinedProfileLikelihood* cpl = combine( like );

  return cpl;

}


CombinedProfileLikelihood* getDMCombinedLikelihood( double mass ) {

  vector <pdfLikelihood*> like;
  like.push_back( getDMLikelihood(mass, 0) );
  like.push_back( getDMLikelihood(mass, 1) );
  like.push_back( getDMLikelihood(mass, 2) );

  CombinedProfileLikelihood* cpl = combine( like );

  return cpl;
}



pdfLikelihood* getTheLikelihood_SR1_calibration( double mass, TString type){
 
 TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
 TString inputDir = xeDir + "SR1/StatisticalAnalyses/inputs_for_likelihood/";

 if(type != "1" && type != "0" )  exit(100);

 int type_int = type.Atoi()  ;

 TString ER_file = "lax_1.5.1/sliced_2d_templates/ERBackgroundModel_SR1_FromSourceCombinedFit180212.root";

 pdfComponent *ER = new pdfComponent("hbkg_Eth_0.00", inputDir + ER_file);
    
    shapeSys *PY = new shapeSys("_py0_");     
    PY->setStep(1.);                      // FIXME
    PY->setMinimum(-2.);   
    PY->setMaximum(2.);    
    PY->setType(FREE_PARAMETER);
    shapeSys *RF = new shapeSys("_rf0_");
    RF->setStep(1.);                     // FIXME
    RF->setMinimum(-2.);   
    RF->setMaximum(2.);    
    RF->setType(FREE_PARAMETER);

    scaleSys *s1 = new scaleSys("ERscale"+type,1.);
    s1->setMinimum(-0.5);
    s1->setMaximum(0.5);
    s1->setType(FREE_PARAMETER);
    ER->addScaleSys(s1);   

    //ER->addShapeSys(Eth);   
    ER->addShapeSys(PY);   
    ER->addShapeSys(RF);

    ER->suffix = "_extruded_yx_r"+type+"_f";
    ER->setEvents( calibration_events[type_int] ); // Set the normalization factor for ER in events.


  TString AC_file = "lax_1.5.1/sliced_2d_templates/Background_wall_ac_templates_v5_SR1_2018-02-28.root";
  pdfComponent *AC = new pdfComponent("acrn220_yx_r"+type+"_f", inputDir +AC_file );

  TString Wall_file = "lax_1.5.1/sliced_2d_templates/Background_wall_ac_templates_v5_SR1_2018-02-28.root";
    pdfComponent *Wall = new pdfComponent("wallrn_1.00_yx_r"+type+"_f", inputDir + Wall_file );

    TString Signal_file = TString::Format("lax_1.5.1/sliced_2d_templates/SignalModel_SR1_%1.0fGeV_FromSourceCombinedFit180212.root",mass) ;

    pdfComponent *Signal = new pdfComponent( "hmc_extruded_yx_r"+type+"_f", inputDir + Signal_file);
    Signal->setScaleFactor( exposure_SR1 );

    TString data_filename = "lax_1.5.1/data/xephyr_Rn220_SR1_pax6.8.0_hax2.4.0_lax1.5.1_cs1LT200_fv1_cuts1_pruned.root";
    TString data_treeName = (type == "0") ? "tree_inner" : "tree_outer";
    dataHandler *data = new dataHandler("calibration",inputDir + data_filename, data_treeName);

  pdfLikelihood *pl = new pdfLikelihood("L_SR1_V"+type, mass);
    
      pl->setExperiment(1);
    pl->addBkgPdfComponent(AC, false);
        pl->addBkgPdfComponent(Wall, false);
    pl->addBkgPdfComponent(ER, true);
    pl->setSignalPdf(Signal);       
    pl->setSignalDefaultNorm(1.E-45);
    //pl.setPrintLevel(WARNING);
    pl->setWithSafeGuard(false);    
    pl->setDataHandler(data);

    pl->initialize();

    pl->printEventSummary();

  return pl;

}
