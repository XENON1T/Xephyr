#include "XeLikelihoods.h"
#include "XePdfObjects.h"
#include "dataHandler.h"
#include "TSystem.h"
#include <vector>

///---------------     PARAMETERS ------------         ///



pdfLikelihood* getTheLikelihood(){

 TString xeDir(gSystem->Getenv("XEPHYR_DIR"));
 TString inputDir = xeDir + "Xephyr/examples/likelihood1D/data/";

 
 pdfComponent *ER = new pdfComponent("hbkg", inputDir + "bkg.root");


    shapeSys *PY = new shapeSys("parA");     
    PY->setStep(1.);                      // FIXME
    PY->setMinimum(-1.);   
    PY->setMaximum(1.);    
    PY->setType(FREE_PARAMETER);
    
    scaleSys *s1 = new scaleSys("ERscale"+type_str,1.);
    s1->setMinimum(-0.5);
    s1->setMaximum(0.5);
    s1->setType(FREE_PARAMETER);
    ER->addScaleSys(s1);   

    //ER->addShapeSys(Eth);   
    ER->addShapeSys(PY);   
    ER->addShapeSys(RF);


  
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
