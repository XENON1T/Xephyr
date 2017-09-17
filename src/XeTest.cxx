#include "XeTest.h"

XeTest t;

XeTest::XeTest() : XeStat("Test package") { 
  cout<<endl
      <<"--------------  Welcome to the XeTest package ------------"<<endl
      <<" An object named \"t\" has been created for interactive      "<<endl
      <<" It is an instance of XeTest, XeMath, XeStat, XeGraphics   "<<endl 
      <<" - Methods can be invoked by e.g. \"t.plotSigmaSI()\"        "<<endl
      <<" - Members can be accessed by e.g. \"t.run->print()\"        "<<endl
      <<"----------------------------------------------------------"<<endl;

  wimp           = NULL ; 
  target         = NULL ;
  interSI        = NULL ; 
  interSD        = NULL ;
  run            = NULL ; 
  simData        = NULL ;
  bands          = NULL ; 
  cuts           = NULL ;
  exclusion      = NULL ; 
  gr             = NULL ;
  mg             = NULL ;
  pv             = NULL ;
  pl             = NULL ;
}

XeTest::~XeTest() {}

/* ----------------------------------------------------------------------*/

void XeTest::plotSigmaSI(double mass, double sigma, int A){

  target=new Target(XENON,A);
  wimp=new Wimp(mass);
  target->setWimp(wimp);
  ErRange Er;
  string title="Rate for "+target->getName()
              +" M_{w}="+formatF(mass,6,1)+" GeV/c^{2}, "
              +" \\sigma_{0} = "+ formatG(sigma,6,1)+" cm^{2}"
              ;
  mg=new XeMultiGraph(title);
  for(int p=0;p<N_COMPARE_SI;p++) {
    target->setInteraction(new SIInteraction(sigma,p)); 
    XeGraph* gRate=target->newGraphOfRate(&Er);
    mg->add(gRate);
  }
  newCanvas("Differential_rate");
  mg->setDefaultYScale(LOG);
  mg->drawWithFrame();
  saveCanvas();
}


/* --------------------------------------------------------------------*/

void XeTest::plotLimitsSI(int rn){

  STATIC_CONST  bool   plotNEvents   =  true;
  STATIC_CONST  bool   plotExclusion =  false;
  STATIC_CONST  bool   plotByLEff    =  false ;
  STATIC_CONST  bool   plotByVEsc    =  false ;

  showTheReferences(false);
  setDebugLevel(1);
  setAnalysisMode(PL_ANALYSIS);
  run=new XeRun(rn,NO_DATA);
  if(!run->printIt(2)) return;
  run->getInteraction()->setSigmaNucleon(1.);
  PRLMassRange mr;

  if(plotNEvents){
    newCanvas("NEvents");
    gr=run->newGraphOfExpectedSignal(&mr);
    gr->drawFrame(AUTOMATIC,AUTOMATIC,.001,5.);
    gr->setLineColor(kRed);      
    gr->draw("L");

    XeGraph* gs1=run->newGraphOfExpectedSignal(&mr);
    gs1->setLineColor(kBlue);      
    gs1->draw("L");
  }

  if(plotByLEff) {
    newCanvas("NEvents by LEff");
    XeMultiGraph* mleff=run->newMultiGraphOfExpectedSignalByLEff(&mr);
    mleff->drawFrame(AUTOMATIC,AUTOMATIC,.001,5.); 
    mleff->draw();
  }

  if(plotByVEsc) {
    newCanvas("NEvents by VEsc");
    XeMultiGraph* mvesc=run->newMultiGraphOfExpectedSignalByVEsc(&mr);
    mvesc->drawFrame(AUTOMATIC,AUTOMATIC,.001,5.); 
    mvesc->draw();
  }


  if(plotExclusion){
    newCanvas("Exclusion");
    XeGraph *gsigma=run->newGraphOfTheoreticalUpperSigma(&mr);
    gsigma->drawFrame(AUTOMATIC,AUTOMATIC,1.E-46,1.E-40);
    gsigma->draw("L");
  }

}


//---------------------------------------------------------------------

void XeTest::plotS1SI(double mass,int rn){

  setDebugLevel(1);

  run=new XeRun(rn,NO_DATA);
  wimp=new Wimp(mass);
  run->setWimp(wimp);
  interSI=new SIInteraction(2.e-45);
  run->setInteraction(interSI);
  if(!run->printIt(1)) return;
  gr=run->newGraphOfSignalERecoilSpectrum();
  gr->drawWithCanvasAndFrame("L");

  mg=run->newMultiGraphOfSignalS1Spectrum(5);
  mg->drawWithCanvasAndFrame();
  

}

//------------------------------------------------------------------

void XeTest::plotSigmaSD(double mass,double sig) {

  ErRange *er=new ErRange(50,0.,50.); 
  for(int x=0;x<N_XENON_MODES;x++){
    target=XeRun::getTarget(XENON_MODES[x]);
    target->setWimpMass(mass);
    for(int n=NEUTRON;n<=PROTON;n++){
      string suffix=" "+SDPureInteraction::getModeName(n)
              +", interaction M_{w} = "
              +formatF(mass,6,1)+" GeV/c^{2}, "
              +" \\sigma_{0}="+ formatLatex(sig,1)+" cm^{2}"
              ;
      string title=target->getFrameName()+suffix;
      string lTitle=target->getLaTeXName()+suffix;
      mg=new XeMultiGraph(lTitle);
      for(int p=0;p<N_SD_INTERACTIONS;p++) {
        target->setInteraction(new SDPureInteraction(sig,n,p)); 
        gr=target->newGraphOfRate(er);
        mg->add(gr);
      }
      newCanvas(makeItAFileName(title));
      mg->drawWithFrame();
    }
  }
}

//------------------------------------------------------------------

void XeTest::plotS1SD(double mass,int runNumber) {

  run=new XeRun(runNumber,NO_DATA);
  run->setWimpMass(mass);
  for(int n=NEUTRON;n<=PROTON;n++){
    string title=run->getName()+" "+SDPureInteraction::getModeName(n)
                +" interaction MWimp= "+formatF(mass,6,1);
    mg=new XeMultiGraph(title);
    for(int p=0;p<N_SD_INTERACTIONS;p++) {
      run->setInteraction(new SDPureInteraction(1.E-40*cm2,n,p));
      XeGraph* rg=run->newGraphOfSignalS1Spectrum();
      mg->add(rg);
    }
    newCanvas(makeItAFileName(title));
    mg->drawWithFrame();
  }
}

//------------------------------------------------------------------------

void XeTest::plotLimitsSD(int runNumber) {

  STATIC_CONST  int    N_COMPARE=N_SD_INTERACTIONS+2;
  STATIC_CONST  int    REFERENCE=SF_SCHWENK_2BC+BAND_CENTER;
               
  XeGraph*            sigma[N_NUCLEONS][N_COMPARE];

  MassRange mr(20,10.,1000.);
  run=new XeRun(runNumber,NO_DATA);
  for(int n=NEUTRON;n<=PROTON;n++){
    string title=run->getName()+" theoretical limit for "
                 +SDPureInteraction::getModeName(n)
                 +" interaction";
    mg=new XeMultiGraph(title);
    for(int p=0;p<N_COMPARE;p++) {
      int w=p;
      if(p>=SF_SCHWENK_2BC) w=SF_SCHWENK_2BC;
      interSD=new SDPureInteraction(1.E-45*cm2,n,w);
      if(p>=SF_SCHWENK_2BC) interSD->setBandEdge(p-SF_SCHWENK_2BC);
      run->setInteraction(interSD);
      sigma[n][p]=run->newGraphOfTheoreticalUpperSigma(&mr,.95,false);
      mg->add(sigma[n][p]);
    }
    newCanvas(title);
    mg->drawWithFrame();
    mg->clear(); 
    saveCanvas();

    XeMultiGraph ma(title+" normalized to Schwenk");
    for(int p=0;p<N_COMPARE;p++) {
      XeGraph* ratio=sigma[n][p]->divide(sigma[n][REFERENCE]);
      ma.add(ratio);
    }   

    ma.newCanvas();
    ma.drawFrame(AUTOMATIC, AUTOMATIC, .0,2.,LINEAR,LOG);
    ma.draw();
    ma.clear();
    saveCanvas();

  }
}

//------------------------------------------------------------------------

void XeTest::plotFFSD(){
  
  STATIC_CONST  int     N_FF=N_SD_POLFF;
  STATIC_CONST  double  YMAX=1.;
  bool                 toPrint=false;
  bool                 plotSA=true;
  bool                 plotFormFactors=false;   
  bool                 plotStructureFactors=false;  
  bool                 plotOriginalBands=false;   
  bool                 printDetailedComputation=false;
  
  SchwenkYRange        syr;
  YRange               yr(51,0.,YMAX);
  SDFormFactor        *Form[N_FF];
  SDStructureFunction *sf[N_SD_SF];

  for(int f=0;f<N_FF;f++)    Form[f]=SDFormFactor::getFormFactor(f);
  for(int p=0;p<N_SD_SF;p++) {
    sf[p]=SDStructureFunction::getStructureFunction(p+FIRST_SD_SF);
  }  

  for(int xe=129;xe<=131;xe+=2) {
    Nucleus *nuc=new Nucleus(XENON,xe);
    string header=nuc->getFrameName();
    string htex=nuc->getLaTeXName();
    double xNorm=1./nuc->ErtoY(1.);
  
    if(plotFormFactors){
      for(int n=NEUTRON;n<=PROTON;n++){
        string title=" Form Factor for "+SDPureInteraction::getModeName(n);
        mg=new XeMultiGraph(htex+title);
        for(int f=0;f<N_FF;f++){
          XeGraph *gf2=Form[f]->newGraphOfFF(nuc,n,&yr);
          gf2->setName(header+title+" "+Form[f]->getName());
          mg->add(gf2);
          if(toPrint) gf2->printGraph(0,false,xNorm);
        }
        newCanvas(header+title);
        mg->drawWithFrame();
        mg->clear(); 
        saveCanvas();
      }
    }

    if(plotSA){

      XeGraph* gsa[N_NUCLEONS][N_FF];
      XeGraph* rsa[N_NUCLEONS][N_FF];

      XeGraph* gsw[N_NUCLEONS][N_SCHWENK_DELTA_A1];
      XeGraph* rsw[N_NUCLEONS][N_SCHWENK_DELTA_A1];

      XeGraph* gse[N_NUCLEONS][N_EDGES];
      XeGraph* rse[N_NUCLEONS][N_EDGES];

      XeGraph* gssf[N_SD_SF][N_NUCLEONS][N_EDGES];
      XeGraph* rssf[N_SD_SF][N_NUCLEONS][N_EDGES];

      XeMultiGraph msa(htex);
      for(int n=NEUTRON;n<=PROTON;n++){
        for(int f=0;f<N_FF;f++){
          gsa[n][f]=Form[f]->newGraphOfSA(nuc,n,&yr);
          msa.add(gsa[n][f]);
        }
        string title=" Structure function for "
                    +SDPureInteraction::getModeName(n);
        msa.setName(htex+title);
        for(int y=0;y<N_SCHWENK_DELTA_A1;y++) {
          gsw[n][y]=SchwenkDeltaA1::newGraphOfSA(nuc,n,y);
          msa.add(gsw[n][y]);
        }
    
        for(int e=0;e<N_EDGES;e++) {
          Form[FF_SCHWENK_2BC]->setBandEdge(e);
          if(plotOriginalBands){
            gse[n][e]=Form[FF_SCHWENK_2BC]->newGraphOfSA(nuc,n,&yr);
            if(toPrint) gse[n][e]->printIt(false);
            msa.add(gse[n][e]);
          }
          for(int p=0;p<N_SD_SF;p++){
            sf[p]->setBandEdge(e);
            gssf[p][n][e]=sf[p]->newGraphOfSA(nuc,n,&yr);
            if(toPrint) gssf[p][n][e]->printIt(false);
            msa.add(gssf[p][n][e]);
          }
        }
  
        newCanvas(header+title);
        msa.drawFrame(0.,YMAX);
        msa.draw();
        saveCanvas();

        msa.clear(); 
        msa.setName(htex+title+", normalized to Schwenk");
        Form[FF_SCHWENK_2BC]->resetCoefficients();
        XeGraph* ref1=gssf[FF_SCHWENK_2BC][n][BAND_CENTER];
        for(int f=0;f<N_FF;f++){
          rsa[n][f]=gsa[n][f]->divide(ref1);
          msa.add(rsa[n][f]);
        }
        for(int y=0;y<N_SCHWENK_DELTA_A1-1;y++) {
          rsw[n][y]=SchwenkDeltaA1::newGraphOfSA(nuc,n,y,SF_SCHWENK_2BC);
          msa.add(rsw[n][y]);
        }
        for(int e=0;e<N_EDGES;e++) {
          if(plotOriginalBands){
            rse[n][e]=gse[n][e]->divide(ref1);
            msa.add(rse[n][e]);
          }
          for(int p=0;p<N_SD_SF;p++){
            rssf[p][n][e]=gssf[p][n][e]->divide(ref1);
            msa.add(rssf[p][n][e]);
          }
        }
        newCanvas(header+title+" normalized to Schwenk");
        msa.drawFrame(0.,YMAX,0.,2.);
        msa.draw();
        msa.clear(); 
        saveCanvas();
      }
    }
 
    if(plotStructureFactors){
      XeGraph* gfc[N_SFACT][N_FF];
      for(int c=0;c<N_SFACT;c++){
        for(int f=0;f<N_FF;f++){
          gfc[c][f]=Form[f]->newGraphOfStructureFactor(nuc,c,&yr);
        }
        string title=SDFormFactor::getStructureFactorName(c);
        string ltitle=SDFormFactor::getStructureFactorLaTeXName(c);
        XeMultiGraph mff(htex+" "+ltitle,N_FF,&gfc[c][0]);
        newCanvas(header+"_"+title);
        mff.drawWithFrame();
        mff.clear(); 
        saveCanvas();
      }
    }
 
    if(printDetailedComputation){
      for(int f=0;f<N_FF;f++){
        Form[f]->printDetailedStructureFactor(nuc);
      }
    }

  }

}
//-------------------------------------------------`------------------

ebLKTest::~ebLKTest(){}

ebLKTest::ebLKTest(DataSet* da):LikelihoodFromDataSet(da){
  bMax=10.;
  addParameter(X0,PARAMETER_OF_INTEREST,"x0",1.,.1,0.,10.);
  addParameter(FB,NUISANCE_PARAMETER   ,"fb",.1,.1,0., 1.);
}

double ebLKTest::computeIndividualLog(double *val){
  double fb=getParameterValue(FB);
  double x0=getParameterValue(X0);
  double lk= (1.-fb)*exp(-val[0]/x0)/x0 + fb/bMax;
  return log(lk);
}

void XeTest::testLK(){
  STATIC_CONST  double fb           =     .2; 
  STATIC_CONST  double x0           =    0.9; 
  STATIC_CONST  int    nSignalEvents= 100000;

  setPrintLevel(3);
  SimulatedExponentAndBackground eb(x0);
  eb.generate(nSignalEvents,(int)(nSignalEvents*fb));

  ebLKTest eblk(&eb);
  eblk.maximize(false);
}


//-------------------------------------------------`------------------
void XeTest::testExclusionPL(int sMin,int sMax){

  STATIC_CONST  double bkg   = 5.;

  setPrintLevel(1);
  setAnalysisMode(PL_ANALYSIS);
  PLcountingSB plsb(1.,bkg);
  exclusion=new Exclusion(&plsb);
  exclusion->initialize();
  EventRange sr(201,0.,20.);

  XeMultiGraph gp("p values");
  XeMultiGraph gc("Corrected p values");
  XeMultiGraph gq("Q for exclusion`");
  for(int s=sMin;s<=sMax;s++) {
    plsb.setNObserved(s);
    exclusion->applyCLs(false);
    gp.add(exclusion->newGraphOfPValue(&sr),s);
    exclusion->applyCLs(true);
    gc.add(exclusion->newGraphOfPValue(&sr));
    gq.add(exclusion->newGraphOfQExclusion(&sr),s);
  }

  newCanvas("P values for Simple Counting PL");
  gp.drawFrame(AUTOMATIC,AUTOMATIC,1.E-6,10.);
  gc.setLineColor(kRed);
  gc.draw();
  gp.setLineColor(kBlack);
  gp.draw();
  gq.setLineColor(kBlue);
  gq.draw();
}

//-------------------------------------------------------------------

void XeTest::studyPoissonCoverage(){
  double CL=.9;
  double bkg=.3;
  PoissonCI* ci[N_CI_MODES];
  for(int mod=0; mod<N_CI_MODES;mod++) {
    ci[mod]=new PoissonCI(1.,0.,mod);
    if(ci[mod]->expectsBackground()) ci[mod]->setBackground(bkg);
  }
  for(double mu=0.01;mu<=10.;mu+=.01) {
    cout<<setw(5)<<mu;
    for(int mod=0; mod<N_CI_MODES;mod++){
      cout<<setw(10)<<ci[mod]->coverage(mu,CL);
    }
    cout<<endl;
  }
}

void XeTest::tableCombinedPL() {

  STATIC_CONST  int N_EXP=2;
  PLcountingSB* pls[N_EXP];

  setPrintLevel(0);

  CombinedProfileLikelihood cpl("Combining two experiments");
  for(int e=0;e<N_EXP;e++){
    double b=0.;
    pls[e]=new PLcountingSB(1.,b);
    pls[e]->setExperiment(e);
    cpl.combine(pls[e]);
  }
  
  cpl.initialize();
  exclusion=new Exclusion(&cpl);

  for(int n1=0;n1<11;n1++) {
    cout<<setw(2)<<n1;
    pls[0]->setNObserved(n1);
    for(int n2=0;n2<11;n2++){
      pls[1]->setNObserved(n2);
      double l=exclusion->computeUpperSigma(0.9);
      cout<<formatF(l,6,2);
    }
    cout<<endl;    
  }
}

void XeTest::tableCombinedCLs(bool withCLs) {

  STATIC_CONST  int N_EXP=2;
  PVcountingSB* pvsb[N_EXP];

  for(int e=0;e<N_EXP;e++){
    double b=2.;
    pvsb[e]=new PVcountingSB(1.,b);
  }
  
  exclusion=new Exclusion((PValue**)pvsb,N_EXP);
  exclusion->applyCLs(withCLs);

  for(int n1=0;n1<11;n1++) {
    cout<<setw(2)<<n1;
    pvsb[0]->setNObserved(n1);
    for(int n2=0;n2<11;n2++){
      pvsb[1]->setNObserved(n2);
      double l=exclusion->computeUpperSigma(0.9);
      cout<<formatF(l,6,2);
    }
    cout<<endl;    
  }
}

void XeTest::tableCLs(bool withCLs){

  double CL=.9;
  PVcountingSB psb(1.);
  exclusion=new Exclusion(&psb);
  exclusion->applyCLs(withCLs);

  for(double bkg=.0;bkg<20.;bkg+=.1){
    psb.setBackground(bkg);
    cout<<setw(10)<<bkg;
    for(int n=0;n<=10;n++){
      psb.setNObserved(n);
      cout<<setw(10)<<exclusion->computeUpperSigma(CL);
    }
    cout<<endl;
  }
}

//  

void XeTest::PLexclusionComputedBackground(bool withCLs, double sigToE) {
 
  STATIC_CONST  int    OBS_MAX =   10;
  STATIC_CONST  double BKG_MAX = 20.0;
  STATIC_CONST  double BKG_STEP=   .1;
  STATIC_CONST  double CL      =   .9;

  setPrintLevel(0);
  PLcountingSB plsb(sigToE);

  exclusion=new Exclusion(&plsb);
  exclusion->applyCLs(withCLs);
 
  for(double bkg=0.;bkg<=BKG_MAX;bkg+=BKG_STEP){
    plsb.setBackground(bkg);
    cout<<"+++"<<setw(5)<<bkg;
    for(int n=0;n<=OBS_MAX;n++){
      plsb.setNObserved(n);
      double l=exclusion->computeUpperSigma(CL);
      cout<<formatF(l,6,2);
    }
    cout<<endl;
  }
}

//-------------------------------------------------------------------

void XeTest::testYellin(int s, int b, int n) {

  double c=Exclusion::typicalUpperLimit(s+b);
  int    nc=c+1;
  TH1F * h=new TH1F("limit","limit",nc,0,nc);

  setPrintLevel(0);
  simData=new SimulatedGaussianAndBackground(1.,.1,2.);
  YellinPValue yel(simData); 
  exclusion=new Exclusion(&yel);

  for(int i=0;i<n;i++){
    simData->generate(s,b);
    double l=exclusion->computeUpperSigma(0.9);
    cout<<"===> "<<simData->getNEvents()<<" events, limit is "<<l<<endl;
    h->Fill(l,1.); 
  }

  newCanvas("Yellin's limits");
  h->Draw();
}


//------------------------------------------------------------------

void XeTest::plotAcceptance(){
  S1Range s1;
  XeGraph* ga[2];
  newCanvas("Acceptances",1,2);

  for(int i=0;i<2;i++){
    bool smeared=i==1;
    string sm=smeared? "smeared S1" : "unsmeared S1";
    newFrame("Cuts on "+sm,0,25,0.,1.1,"S1 (pe)","Acceptance");
    ga[i]=run->newGraphOfSelectionCutsAcceptance(&s1,smeared);
    ga[i]->Draw("L");
    nextSubCanvas();
  }
}

//------------------------------------------------------------------
void XeTest::setupTheRun(bool publishedBackground, double mass ,int runNumber){
  run= new XeRun(runNumber,REAL_DATA);
  if(mass>0.) run->setWimpMass(mass);
  if(publishedBackground)  run->setPublishedElectronBackgroundModel();
}

bool XeTest::printTheSetup(){
  if(run==NULL) cout<<"You didn't initialize the run!"<<endl;
  else if(run->printIt()) return run->printCandidates(5);
  return false;
}

//------------------------------------------------------------

void XeTest::showRunData(bool publishedER, double mass,int rn
                        ,int analysisM, int displayM) {
  setDebugLevel(0);
  setPrintLevel(0);
  setAnalysisMode(analysisM);
  setupTheRun(publishedER,mass,rn);

  if(!printTheSetup()) return;
  S1S2Display::setMode(displayM);
  cuts=run->getAllCuts();
  newCanvas("Data "+run->getName(),2,2);
  for(int type=0;type<N_DATA_TYPES;type++) {
    S1S2Data* s1s2=run->getS1S2Data(type);
    s1s2->drawS1S2(CURRENT_DISPLAY,false);
    cuts->draw();
    nextSubCanvas();
  }

  if(analysisMode!=PL_ANALYSIS) return;

  S1S2Bands* ref=run->getReferenceBands();
  const int nTypes=3;
  const int types[nTypes]={AM_BE_DATA,E_GAMMA_DATA,DM_DATA};

  for(int ty=0;ty<nTypes;ty++){
    int before=types[ty];
    int after=S1S2Data::typeAfterCut(before);
    newCanvas(S1S2Data::getTypeName(before),2,2);

    S1S2Data* s1s2=run->getS1S2Data(before);
    s1s2->drawS1S2();
    cuts->draw();
    ref->drawBands();

    bands=run->getS1S2Bands(after);
    nextSubCanvas();
    bands->drawCellContentWithFrame(LOG);
    nextSubCanvas();
    bands->drawBandContentWithFrame(LOG);

    nextSubCanvas();
    bands->newHistogramOfS1Spectrum(ALL,LINEAR);
  }

  newCanvas("Background bands",2,2);
  S1S2Display::setMode(BAND_VS_SLICE);
  for(int ty=FIRST_BACKGROUND;ty<=LAST_BACKGROUND;ty++){
    run->getS1S2Bands(ty)->drawS1S2();
    nextSubCanvas();
  }

  newCanvas("Background S1",2,2);
  for(int ty=FIRST_BACKGROUND;ty<LAST_BACKGROUND;ty++){
    run->newGraphOfBackgroundS1Distribution(ty,LINEAR);
    nextSubCanvas();
  }
  mg=run->newMultiGraphOfAllBackgroundsInBandsS1Distribution(LINEAR);

}

//-----------------------------------------------------------------

void XeTest::setupTheAnalysis(int mode
                             , bool doFitSystBkTvalue
                             , bool doFitLeffTValue
                             , bool withS1Shape
                             ){
  setAnalysisMode(mode);
  if(analysisMode==PL_ANALYSIS){
    pl=new S1S2PL(run);
    pl->fitLEffTValue(doFitLeffTValue);
    pl->fitSystBkgTValue(doFitSystBkTvalue);
    pl->withS1Likelihood(withS1Shape);
    if(pl->isInError()) return;
    pv=pl;
  }
  else {
    pv=new S1S2CLs(run);
  }
  exclusion=new Exclusion(pv);
}

// ---------------------------------------------------------------

void XeTest::plotLogLL(bool s1,bool publishedBkg,bool doFitSystBkTvalue 
                     ,double specialMass,int nMasses,int rn) {

  bool doFitLeffTValue      =  true;
  analysisMode              =  PL_ANALYSIS;
  setupTheRun(publishedBkg,specialMass,rn);
  setupTheAnalysis(analysisMode,doFitSystBkTvalue,doFitLeffTValue,s1);
  if(!printTheSetup()) return;
  mg=new XeMultiGraph("Log LL",LABEL_EVT,LABEL_LOG_LIKELIHOOD);

  MassRange *mr;
  if(specialMass>0.) mr=new MassRange(1,specialMass,specialMass);
  else               mr=XeRange::getDefaultMassRange();
  if(nMasses==0) nMasses=mr->getNPoints();
  vector<XeGraph*> grs(nMasses);
  for(int g=0;g<nMasses;g++){
    double mass=mr->getValue(g);
    cout<<"Setting mass to "<<mass<<endl;
    pl->setWimpMass(mass);
    exclusion->prepareForLimit();
    grs[g]=pl->newGraphOfLogLikelihood();
    grs[g]->setMassLegend(mass);
  }
  string title="Log LL";
  mg=new XeMultiGraph(title,grs,mr->getTheValues(),LABEL_EVT
                     ,LABEL_LOG_LIKELIHOOD);
  newCanvas(title);
  mg->setRainbowColors();
  mg->drawWithFrame("L");

}

//-----------------------------------------------------------------


void XeTest::runAnalysis(bool s1,bool publishedBkg ,bool doFitSystBkTvalue 
                        ,double specialMass,int rn,bool LimitSI,bool LimitSD) {

  bool doFitLeffTValue      =  true;
  analysisMode              =  PL_ANALYSIS;
  bool showPGraph           =  false;

  bool debugOneCaseInDetail =  specialMass>0.;
  
  setDebugLevel(0);
  setPrintLevel(1);

  MassRange mr(25);

  setupTheRun(publishedBkg,specialMass,rn);
  setupTheAnalysis(analysisMode,doFitSystBkTvalue,doFitLeffTValue,s1);
  if(!printTheSetup()) return;
  MethodCounter::reset();

  if(debugOneCaseInDetail) {
    setDebugLevel(1);
    setPrintLevel(2);
    exclusion->setWimpMass(specialMass);
    double s=exclusion->computeCrossSection()->getEvents();
    cout<<" ==========> "<<s<<endl;
    MethodCounter::list();
    return;
  }

  if(showPGraph){
    EventRange sr(21,0.,20.);
    pv->setWimpMass(specialMass);
    exclusion->newGraphOfPValue(&sr,AUTO);
  }

  if(LimitSI){
    XeLimits* limits=exclusion->newUpperLimits(&mr);
    XeGraph* gsl=limits->newGraph(UPPER_LIMIT,SIGMA_UNIT);
    XeGraph* gse=limits->newGraph(ESTIMATED,SIGMA_UNIT);
    XeMultiGraph* mgs=new XeMultiGraph("Sigma");
    mgs->add(gsl,kBlue,"Upper limit");
    mgs->add(gse,kBlack,"Estimated");
    mgs->drawWithCanvasAndFrame("L");

    XeGraph* gel=limits->newGraph(UPPER_LIMIT,EVENT_UNIT);
    XeGraph* gee=limits->newGraph(ESTIMATED,EVENT_UNIT);
    XeMultiGraph* mge=new XeMultiGraph("Number of events");
    mge->add(gel,kBlue,"Upper limit");
    mge->add(gee,kBlack,"Estimated");
    mge->drawWithCanvasAndFrame("L");

    limits->printIt(1); 
  }

  if(LimitSD){ 
    for(int n=NEUTRON;n<=PROTON;n++){
      string title=run->getName()+" Sigma limit for "
                   +SDPureInteraction::getModeName(n) +" interaction";
      mg=new XeMultiGraph(title);
      for(int p=0;p<N_SD_INTERACTIONS+2;p++) {
        int w=p;
        if(p>=SF_SCHWENK_2BC) w=SF_SCHWENK_2BC;
        interSD=new SDPureInteraction(1.E-45*cm2,n,w);
        if(p>=SF_SCHWENK_2BC) interSD->setBandEdge(p-SF_SCHWENK_2BC);
        run->setInteraction(interSD);
        XeGraph* gd=exclusion->newGraphOfUpperSigma(&mr);
        gd->setStyle(interSD->getStyleFromFF());
        gd->setLegend(interSD->getLegendFromFF());
        mg->add(gd);
      }
      newCanvas(title);
      mg->drawWithFrame();
      mg->clear();
    }
  } 
  MethodCounter::list();
}

//----------------------------------------------------------------------
void XeTest::combineRuns(int nr) {

  STATIC_CONST  int N_EXP=3;
  S1S2PL* pls[N_EXP];
  static const int runs[N_EXP]={8,10,12};

  setDebugLevel(0);
  setPrintLevel(1);

  string sr="Combination of runs";
  for(int r=0;r<nr;r++) sr=sr+formatI(runs[r],3);
  CombinedProfileLikelihood cpl(sr);
  
  for(int r=0;r<nr;r++){
    setupTheRun(true,0.,runs[r]);
    pls[r]=new S1S2PL(run);
    pls[r]->fitLEffTValue(true);
    pls[r]->fitSystBkgTValue(false);
    pls[r]->withS1Likelihood(true);
    if(pls[r]->isInError()) return;
    cpl.combine(pls[r]);
  }

  MassRange mr(25);
  cpl.initialize();
  exclusion=new Exclusion(&cpl);
  XeLimits* limits=exclusion->newUpperLimits(&mr);
  XeGraph* gsl=limits->newGraph(UPPER_LIMIT,SIGMA_UNIT);
  XeGraph* gse=limits->newGraph(ESTIMATED,SIGMA_UNIT);
  XeMultiGraph* mgs=new XeMultiGraph("Sigma");
  mgs->add(gsl,kBlue,"Upper limit");
  mgs->add(gse,kBlack,"Estimated");
  mgs->drawWithCanvasAndFrame("L");

  XeGraph* gel=limits->newGraph(UPPER_LIMIT,EVENT_UNIT);
  XeGraph* gee=limits->newGraph(ESTIMATED,EVENT_UNIT);
  XeMultiGraph* mge=new XeMultiGraph("Number of events");
  mge->add(gel,kBlue,"Upper limit");
  mge->add(gee,kBlack,"Estimated");
  mge->drawWithCanvasAndFrame("L");

  limits->printIt(1); 
 
}


//-------------------------------------------------------------------------
void XeTest::computeSensitivityBands( bool s1, bool publishedBkg
                                   , int nSimulations, int runNumber) {
  bool doFitSystBkTvalue=true;
  bool doFitLeffTValue=true;

  setupTheRun(publishedBkg,0.,runNumber);
  setupTheAnalysis(PL_ANALYSIS,doFitSystBkTvalue,doFitLeffTValue,s1);
  if(!printTheSetup()) return;

  exclusion=new Exclusion(pl);
  SensitivityBands* sBands=exclusion->newSensitivityBands(nSimulations);
  sBands->printIt(2);
}

//-------------------------------------------------------------------------
void XeTest::studyQValuesForBackground( int nm, int nSimulations, bool s1
                              , bool publishedBkg, int runNumber) {

  STATIC_CONST  int nMasses=5;
  double masses[nMasses]={6.,10.,20.,50.,1000.};
  bool doFitSystBkTvalue=false;
  bool doFitLeffTValue=true;

  if(nm==0) nm=nMasses;

  setPrintLevel(0);
  setupTheRun(publishedBkg,0.,runNumber);
  setupTheAnalysis(PL_ANALYSIS,doFitSystBkTvalue,doFitLeffTValue,s1);
  if(!printTheSetup()) return;

  XeValues* qValues=new XeValues("Q values for background");
  XeValues* sigmaValues=new XeValues("Estimated limits for background");
  XeValues* eventValues=new XeValues("Estimated N Events for background");
  mg=new XeMultiGraph("Q values for background");
  XeMultiGraph *mc=new XeMultiGraph("Cumulated Q values for background");
  XeMultiGraph *me=new XeMultiGraph("Estimated nEvents for background");

  Chi2Dist chi2(2);
  STATIC_CONST  int    Q_BINS=100;
  STATIC_CONST  double Q_MAX =5.;
  TabulatedDist chi2d("chi2",Q_BINS,0.,Q_MAX);
  chi2d.importDistribution(&chi2);
  chi2d.addDelta(0.,1.);  // yes, weight is 1.
  XeGraph* gc=chi2d.newGraph("Chi2","",NULL,"\\chi^2");
  gc->setLineColor(kBlack);
  XeGraph* gcc=chi2d.newCumulatedGraph("Chi2","",NULL,"\\chi^2");
  gcc->setLineColor(kBlack);

  for(int m=0;m<nm;m++){
    double mass=masses[m];
    cout<<"Now "<<mass<<endl;
    run->setWimpMass(mass);
    if(!run->checkAll()) return;
    exclusion-> simulateNQValues(nSimulations,qValues,sigmaValues,eventValues);
    string leg="Mass = "+formatF(mass,7,1);

    TabulatedDist tdist("q",Q_BINS,0,Q_MAX,qValues);
    XeGraph* g=tdist.newGraph("Q","density",NULL,leg);
    mg->add(g);
    XeGraph* c=tdist.newCumulatedGraph("Q","density",NULL,leg);
    mc->add(c);

    TabulatedDist edist("nEvents",100,-20,20, eventValues);
    XeGraph* h=edist.newGraph("nEvents","Density",NULL,leg);
    me->add(h);
   }
  newCanvas("Background Sigma and Q values");

  subCanvas(2,2);
  mg->setRainbowColors();
  mg->add(gc);
  mg->setDefaultYScale(LOG);
  mg->drawWithFrame("L");

  nextSubCanvas();
  mc->setRainbowColors();
  mc->add(gcc);
  mc->drawWithFrame("L");

  nextSubCanvas();
  me->setRainbowColors();
  me->drawWithFrame("L");
}
