#include "XePdfObjects.h"



scaleSys::scaleSys(TString name, double relativeUncertainty) : LKParameter(PAR_NOT_ASSIGNED, NUISANCE_PARAMETER, name.Data(), 0, 0.01, -5.,5.) {


	relUnc = relativeUncertainty;
	isNull = false;
}

//scaleSys::~scaleSys(){

//}

double scaleSys::getNormModifier(){
	//this case return a sys that is centered in zero 
	//and for a tvalue=0 have zero events. Histo is supposed to be normalized to 1
	if(isNull) return(getCurrentValue() * relUnc);

	return ( 1. +  getCurrentValue() * relUnc ) ; 
}



shapeSys::shapeSys(TString name) : LKParameter(PAR_NOT_ASSIGNED, NUISANCE_PARAMETER, name.Data(), 0, 1., -1.,1.) {


}

//shapeSys::~shapeSys(){

//}

double shapeSys::getNearestLow(){

	if( getCurrentValue() > getMaximum() || getCurrentValue() < getMinimum() ) {
		// this should never happen, minuit is limited in range
		cout << "shapeSys::getNearestLow() : ERROR - sys value outside range " << endl;
		exit(100);
	}


	double nearestLow = ( (int) ( getCurrentValue() / getStep() ) -1 ) * getStep() ;

	if(getCurrentValue() > 0 )
		nearestLow = ( (int) ( getCurrentValue() / getStep() ) ) * getStep() ;


	if( getCurrentValue() == getMinimum() ) 
		nearestLow = getMinimum() ; 

	if( getCurrentValue() == getMaximum() ) 
		nearestLow = getMaximum() - getStep() ; 

	return nearestLow ;
}

double shapeSys::getNearestHigh(){

	if( getCurrentValue() > getMaximum() || getCurrentValue() < getMinimum() ) {
		// this should never happen, minuit is limited in range
		cout << "shapeSys::getNearestHigh() : ERROR - sys value outside range " << endl;
		exit(100);
	}

	return ( getNearestLow()  + getStep() ) ;


}

pdfComponent::pdfComponent(TString name, TString filename) : errorHandler("pdfComponent"), pdf_name(name) {

  	file = TFile::Open(filename);
	if(file == NULL) {
		cout << "pdfComponent::pdfComponent:: ERROR - can't access file " << filename << endl;
		exit(100);
	}

	histos.push_back(NULL);

	InterpFactors.push_back(0.);

	old_t_val.push_back(-999.);

	defaultDistro = NULL;

	scaleFactor   = -999.;

	suffix = "";

	doExtend  = false;
}

pdfComponent::~pdfComponent(){
	for(unsigned int k=0 ; k < myScaleUnc.size(); k++) delete myScaleUnc[k];
	for(unsigned int k=0 ; k < myShapeUnc.size(); k++) delete myShapeUnc[k];

	myScaleUnc.clear();
	myShapeUnc.clear();
	
	delete defaultDistro;

	InterpFactors.clear();
 	
	old_t_val.clear();	

	for(unsigned int k=0 ; k < histos.size(); k++) {
			delete histos[k];
	}
	
	histos.clear();

	file->Close();
	delete file;

}

void pdfComponent::loadHistos() {
	
	//clear vector of pointers, this does not delete the histo
	//from memory, they remain attached to the TFile, this is a wanted
	//feature, we don't hit the disk each time, we put in memory all the
	//histo that have been read for later access. Histo are not duplicated. 
	histos.clear();  

	InterpFactors.clear();

	loadDefaultHisto();

	
	//end here if no shape sys
	if(myShapeUnc.size() ==0) return;
	
	// total number of combination of point needed for 
	// interpolation, 2^N_parameter. Each parameter can be loaded
	// with his closest high/low end on the grid of that parameter.
        // example: to interpolate Leff=0.5 we need to load Leff=0 and Leff=1
	// histograms.	
	int N = pow(2, myShapeUnc.size());
	
	//the idea is to compute the volume of the hypercube in parameter space
	//corresponding to that grid point and divide by the total volume,
	//this is the interpolation factor.
	double total_vol = 1.;

	//compute total volume
	for(unsigned int k =0; k< myShapeUnc.size(); k++){
		total_vol *= fabs( myShapeUnc[k]->getNearestHigh() - 
				myShapeUnc[k]->getNearestLow() );
	}

	//loop on the realization of parameter, we create a vector of bool
	//that say for each parameter what "step" (point in grid) hast to be loaded,
	//the low end or high end wrt the value asked.
	for(int i=0; i < N; i++){
		
		vector <bool> grid_point;

		double grid_point_vol = 1.;

		//loop on the parameter
		for(unsigned int k =0; k< myShapeUnc.size(); k++){
			
			int param_power = pow(2,k);
			//if false load the histo corresponding to the low end
			//of that parameter, the high end otherwise
			int lowOrHigh = ( i / param_power) % 2 ;


			if(lowOrHigh == 0)       grid_point.push_back(false);
			else if(lowOrHigh == 1)  grid_point.push_back(true);

			else Error("loadHistos","something very wierd happened");

			//compute the numerator of iterpolation factor 
			//for this grid point, area of the opposite
			if(lowOrHigh == 1) 
			    grid_point_vol *= fabs( myShapeUnc[k]->getCurrentValue() - 
					          myShapeUnc[k]->getNearestLow() );
			else 
			    grid_point_vol *= fabs( myShapeUnc[k]->getCurrentValue() - 
					          myShapeUnc[k]->getNearestHigh() );

		}

		TString histName = getNearestHistoName(grid_point);
//		Info("loadHistos","Loading histo: "+ histName + " interp. factor for  " + getParamValueString() + " : " +
//				+ TString(printTools::formatF(grid_point_vol / total_vol)) ) ;
		
		//check if name exist
		if( file->FindKey(histName) == NULL) 
			Error("loadHistos","Histogram does not exist in file: "+histName);
		//store histo pointer
		histos.push_back((TH2F*)file->Get(histName));


		//store the interpolation factor
		InterpFactors.push_back( grid_point_vol / total_vol );

	}	


}

void pdfComponent::loadDefaultHisto(){

  if(defaultDistro == NULL) {	

       if( file->FindKey(getDefaultHistoName()) == NULL) 
	Error("loadDefaultHisto","histo name " +getDefaultHistoName() + "  not found");

	defaultDistro    = (TH2F*)file->Get(getDefaultHistoName());

  }
	
	
}


TString pdfComponent::getNearestHistoName(vector<bool> setOfVal){

  // This part is going to be revisited once add multidimensional interpolation capability
  TString name_histo(pdf_name);

  if(setOfVal.size() != myShapeUnc.size())
	Error("getNearestHistoName","SetOfVal not compatible with number of shape sys");

  for(unsigned int j=0; j < myShapeUnc.size(); j++){

	  char value_temp[20];

	  TString sysName(myShapeUnc[j]->getName());
	  name_histo.Append( sysName );

	  //name_histo.Append("_" + sysName );
	  
	  //get nearest low or high value to the current one
	  if(!(setOfVal[j])) sprintf(value_temp,"%.2f%s", myShapeUnc[j]->getNearestLow(),suffix.Data());
	  else      sprintf(value_temp,"%.2f%s", myShapeUnc[j]->getNearestHigh(), suffix.Data());

	  name_histo.Append(value_temp);
   }

   return name_histo;
}

TString pdfComponent::getDefaultHistoName(){

  TString name_histo(pdf_name);

  for(unsigned int j=0; j < myShapeUnc.size(); j++){

	  char value_temp[20];

	  TString sysName(myShapeUnc[j]->getName());

	  name_histo.Append( sysName );

	  //name_histo.Append("_" + sysName );

	  // default value of all sys is zero.
	  sprintf(value_temp,"%.2f%s", 0.,suffix.Data());

	  name_histo.Append(value_temp);
   }

   return name_histo;


}



double pdfComponent::getNormalizedDensity(double s1, double s2) {
	
	//load histogram according to the current value of the parameters
	loadHistos();

	int s1_bin = defaultDistro->GetXaxis()->FindBin(s1);
	int s2_bin = defaultDistro->GetYaxis()->FindBin(s2);

	//use single histo if no shape uncertainties
	double interpolated_content = defaultDistro->GetBinContent(s1_bin, s2_bin);

	//use single histo if no shape uncertainties
	if(myShapeUnc.size() > 0) {

            interpolated_content = 0.;
	    for(unsigned int k=0; k< histos.size(); k++){	
		interpolated_content += histos[k]->GetBinContent(s1_bin, s2_bin) * InterpFactors[k]; 
	    }

	 }

	//scale uncertainty part
	for(unsigned int k=0; k < myScaleUnc.size() ; k++){

		interpolated_content *= myScaleUnc[k]->getNormModifier(); 
	}


	if(scaleFactor > 0.) interpolated_content *= scaleFactor;

	return interpolated_content;
	
}

double pdfComponent::getDefaultDensity(double s1, double s2){

	loadDefaultHisto();

	int s1_bin = defaultDistro->GetXaxis()->FindBin(s1);
	int s2_bin = defaultDistro->GetYaxis()->FindBin(s2);

	double content =  defaultDistro->GetBinContent(s1_bin,s2_bin);

	if(scaleFactor >0. ) content *= scaleFactor;

	return content;

}




double  pdfComponent::getNormalizedEvents() {

	//load histogram according to the current value of the parameters
	loadHistos();

	//use default histo if no shape uncertainties
	double all_content = defaultDistro->Integral();

	//use single histo if no shape uncertainties
	if(myShapeUnc.size() > 0 ) { 
            all_content = 0.;
	    for(unsigned int k=0; k< histos.size(); k++){	
		all_content += histos[k]->Integral() * InterpFactors[k]; 
	    }
	}

	//scale uncertainty part
	for(unsigned int k=0; k < myScaleUnc.size() ; k++){

		all_content *= myScaleUnc[k]->getNormModifier(); 
	}

	if(scaleFactor > 0.) all_content *= scaleFactor;

	return all_content;
}


double  pdfComponent::getDefaultEvents(){

	loadDefaultHisto();

	double integral = defaultDistro->Integral();

	if(scaleFactor > 0.) integral *= scaleFactor;

	return integral;

}



TH2F   pdfComponent::getInterpolatedHisto(){

	//load histogram according to the current value of the parameters
	loadHistos();

	//clone default
	TH2F h_temp = getDefaultHisto() ;
	h_temp.SetName("Interp_" + getParamValueString() );


	if(myShapeUnc.size() > 0) {
	    h_temp = *histos[0];
            h_temp.Reset();
	    for(unsigned int k=0; k< histos.size(); k++)
		h_temp.Add(histos[k], InterpFactors[k]);

	    //getDefault is scaled, but reset function bring back content to zero
	    if(scaleFactor > 0.) h_temp.Scale(scaleFactor);
	}

	//scale uncertainty part
	for(unsigned int k=0; k < myScaleUnc.size() ; k++){

		h_temp.Scale(myScaleUnc[k]->getNormModifier()); 
	}

	if(h_temp.GetNbinsX() == 63 || doExtend) 
		extendHisto(h_temp);

	return h_temp;	
}

TH2F   pdfComponent::getDefaultHisto(){

	loadDefaultHisto();

	TH2F h_temp(*defaultDistro);   //we don't want to return a pointer to local variable
	h_temp.SetName("DefaultHisto_" + pdf_name);

	if(scaleFactor > 0.) h_temp.Scale(scaleFactor);

	if(h_temp.GetNbinsX() == 63 || doExtend) 
		extendHisto(h_temp);

	return h_temp;
}





double pdfComponent::getDefaultPdfIntegral(double s1_min, double s1_max, double s2_min, double s2_max){
	
	loadDefaultHisto();

	return dataHandler::integrate(defaultDistro,s1_min,s1_max,s2_min,s2_max);

}



TString pdfComponent::getParamValueString(){
	
	TString gridPointName = "";
	if(myShapeUnc.size() ==0 && myScaleUnc.size() ==0) return TString("NONE");

	for(unsigned int k =0; k< myShapeUnc.size(); k++){
		if(k!=0) gridPointName.Append("_");
		gridPointName.Append(myShapeUnc[k]->getName());
		gridPointName.Append(printTools::formatF(myShapeUnc[k]->getCurrentValue(),1,2));
	}

	for(unsigned int k =0; k< myScaleUnc.size(); k++){
		if(gridPointName!="") gridPointName.Append("_");
		gridPointName.Append(myScaleUnc[k]->getName());
		gridPointName.Append(printTools::formatF(myScaleUnc[k]->getCurrentValue(),1,2));
	}

	return gridPointName;

}

TString pdfComponent::getParamValueWritable(){
	
	TString gridPointName = "";
	if(myShapeUnc.size() ==0 && myScaleUnc.size() ==0) return TString("NONE");

	for(unsigned int k =0; k< myShapeUnc.size(); k++){
		if(k!=0) gridPointName.Append("_");
		gridPointName.Append(myShapeUnc[k]->getName());
		gridPointName.Append(printTools::formatR(myShapeUnc[k]->getCurrentValue(),2));
	}

	for(unsigned int k =0; k< myScaleUnc.size(); k++){
		if(gridPointName!="") gridPointName.Append("_");
		gridPointName.Append(myScaleUnc[k]->getName());
		gridPointName.Append(printTools::formatR(myScaleUnc[k]->getCurrentValue(),2));
	}

	return gridPointName;

}


void pdfComponent::plotInterpolatedSpace(bool doProjectionX, double min, double max, int Nsteps){
	
	if(myShapeUnc.size() == 0 )
		Error("plotInterpolatedSpace","no shape uncertainties, no interpolation, no scan possible... Go home");

	

	//the +1 is to scan also the min and max values
	int nGridPoints = pow(Nsteps + 1,myShapeUnc.size());

	vector <TCanvas*>    canvases;
	vector <TH1D*>       projections;
	vector <TLegend*>    legends;

	unsigned int canvas_counter 	= 0;
	unsigned int projection_counter = 0;
	

	for(int i =0; i < nGridPoints; i++){

		if( ((int)projection_counter %  (Nsteps+1)) == 0) {


			canvas_counter++;
			canvases.push_back(new TCanvas());

			legends.push_back(new TLegend(0.7,0.9,0.99,0.7));
		}
		
		//set the value of all shape sys to one grid point of hyperspace
		for(unsigned int k=0; k< myShapeUnc.size(); k++){

			double step_size = (myShapeUnc[k]->getMaximum() - myShapeUnc[k]->getMinimum() ) / ((double)Nsteps); 
			int sys_power = pow(Nsteps+1,k);

			int step_number = ( i / sys_power ) % (Nsteps+1);

			double value_sys = (double)step_number * step_size + myShapeUnc[k]->getMinimum();

			myShapeUnc[k]->setCurrentValue(value_sys);	
		}

		TH2F h_temp(getInterpolatedHisto());
		
		TH1D *project_temp =  NULL;

		int bin_min=0;
		int bin_max =0 ;

		if(doProjectionX){
			bin_min      = h_temp.GetYaxis()->FindBin(min);
			bin_max      = h_temp.GetYaxis()->FindBin(max);
			project_temp = h_temp.ProjectionX("prX_"+getParamValueString(),bin_min, bin_max);
		}

		else{ 
			bin_min      = h_temp.GetXaxis()->FindBin(min);
			bin_max      = h_temp.GetXaxis()->FindBin(max);
			project_temp = h_temp.ProjectionY("prY_"+getParamValueString(),bin_min, bin_max);
		}


		//set the canvas
		canvases[canvas_counter -1]->cd();

		project_temp->SetLineColor(projection_counter+1);

		if( ((int)projection_counter %  (Nsteps+1)) == 0) project_temp->Draw();
		else     project_temp->Draw("same");

		legends[canvas_counter-1]->AddEntry(project_temp, getParamValueString());
		legends[canvas_counter-1]->Draw("same");


		//storing the projections
		projections.push_back(project_temp);

		projection_counter++;
			
	}

	//clean up and store in plots
	for(unsigned int j=0; j < canvases.size(); j++){

		if(j == canvases.size() -1) canvases[j]->Print("sys_shapes_iterpolations.pdf)");
		else canvases[j]->Print("sys_shapes_iterpolations.pdf(");

		delete canvases[j];

		delete legends[j];

	}


	canvases.clear();
	legends.clear();

	TFile proj_file("projection_sys_interpolated.root","RECREATE");

	for(unsigned int j=0; j < projections.size(); j++){
		projections[j]->Write();
		delete projections[j];

	}
	
	proj_file.Close();


}

void pdfComponent::extendHisto(TH2F &h){

	TH2F h_right(TString(h.GetName())+"_yep","",67,3,70,70,log10(50),log10(8000));

	for(int x=1; x <= h.GetNbinsX(); x++){
		for(int y=1; y <= h_right.GetNbinsY(); y++){
			double content = h.GetBinContent(x,y);
			if (content <=0 ) content = 1.E-8;
			h_right.SetBinContent(x,y, h.GetBinContent(x,y));
		}
	}

	h = h_right;	
}



scaleSys* pdfComponent::getScaleSys(TString search_name){
	
	for(unsigned int i=0; i< myScaleUnc.size(); i++){
		if(myScaleUnc[i]->getName() == search_name ) return myScaleUnc[i];
	}

	Error("getScaleSys", search_name + " not found.");

	return NULL;
}

shapeSys* pdfComponent::getShapeSys(TString search_name){
		
	for(unsigned int i=0; i< myShapeUnc.size(); i++){
		if(myShapeUnc[i]->getName() == search_name ) return myShapeUnc[i];
	}

	Error("getShapeSys", search_name + " not found.");

	return NULL;
}

void pdfComponent::replaceUncertainty(TString name, scaleSys* newScale){
	
	bool found = false;

	for(unsigned int i=0; i< myScaleUnc.size(); i++){
		if(myScaleUnc[i]->getName() == name ) {	
			myScaleUnc[i] = newScale;
			found = true;
		}
	}

	if(found == false) Error("replaceUncertainty", name + " not found.");

}


void pdfComponent::replaceUncertainty(TString name, shapeSys* newShape){

	bool found = false;

	for(unsigned int i=0; i< myShapeUnc.size(); i++){
		if(myShapeUnc[i]->getName() == name ) {	
			myShapeUnc[i] = newShape;
			found = true;
		}
	}

	if(found == false) Error("replaceUncertainty", name + " not found.");


}



/////-------------------  class for comparison   --------//

histoCompare::histoCompare():errorHandler("histoCompare"){
	
	  
	    rebinX = 1;
	    rebinY = 1;
	    projectionX = true;
	    doStack      = true;
	    projectionMin =  UNDEFINED;
	    projectionMax =  UNDEFINED;

	    binMin = 1;
	    binMax = -1;

	    names.push_back("");	    
	    projectedBase = NULL;

}

histoCompare::~histoCompare(){}




void histoCompare::compare(){

    if(base.GetEntries() == 0) 
	    Error("draw()","Base histo is not set or empty, use  setBaseHisto(TH2F )");

    if(compareList.size() == 0)
	    Error("draw()","Compare Histo list  is not set, use addHistoToList(TH2F )");

    // do the projection and fill projectedBase and projectedList
    project(); 
     
    if(doStack){
	//Stacked Plot
		gStyle->SetPalette(91);
		//gStyle->SetPalette(57);

        projectedBase->Draw("PE");
	
        for(unsigned int i=0; i < projectedList.size(); i++) {

	     	//stacking the histos
	     	if(i != 0) projectedList[i]->Add(projectedList[i-1]);
	    
   		projectedList[i]->SetLineColor(i+2); 
   		setOptions(projectedList[i], false);
	}
		// draw them in inverse order so can put colors on top of each other
        for(unsigned int i= 1; i <= projectedList.size(); i++) {
			unsigned int n = projectedList.size();
			projectedList[n - i]->Draw("PLC PFC hist Same");
		}

        projectedBase->SetLineColor(1);
		setOptions(projectedBase, true);
        projectedBase->Draw("samePE");
        drawLegend(projectedBase, projectedList);
	
    } 
    else{
	//just comparison
	projectedBase->Draw("hist");
        for(unsigned int i=0; i < compareList.size(); i++) {
   		 setOptions(projectedList[i], false);
   		 projectedList[i]->SetLineColor(i+2); 
		 projectedList[i]->Draw("sameHIST");
	}

	projectedBase->SetLineColor(1);
	setOptions(projectedBase, false);
	projectedBase->Draw("sameHIST");

        drawLegend(projectedBase, projectedList);

    }

    //check overflow, problem: due to re-binning, some bins may go to overflow
    double base_overflow = projectedBase->GetBinContent(projectedBase->GetNbinsX() +1);
    if(base_overflow)
	    cout << "WARNING:: You got " << base_overflow << " data events in overflows -- Check your rebinning! "<< endl;

	gPad->RedrawAxis();	
}


void histoCompare::compareWithRatio(){
	

  TCanvas *c1 = new TCanvas("c1","ratio plot",800,800);
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35 );
  pad1->SetTopMargin(0.05);
  pad1->SetRightMargin(0.05);
  pad2->SetTopMargin(0.1);
  pad2->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.30);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();

  //draw the top part
  compare();

  // fix the label offset 
  projectedBase->GetXaxis()->SetTitleOffset(1.);
  projectedBase->GetXaxis()->SetLabelOffset(1.);

  pad1->Update();

  //setting up the ratios and draw
  pad2->cd();
  vector <TH1D*>  ratios;

  if(doStack){
	  // if stack plot then only ratio between base and top histo
	  unsigned int n = projectedList.size();

  	  ratios.push_back((TH1D*)projectedList[n-1]->Clone("ratio"));
	  ratios[0]->Divide(projectedBase);
	  setOptions(ratios[0], true, true);
	  ratios[0]->Draw("PE");
  }

  else{

       for(unsigned int k=0; k < projectedList.size(); k++){

  	  ratios.push_back((TH1D*)projectedList[k]->Clone(Form("ratio_%d",k)));
	  ratios[k]->Divide(projectedBase);

  	  //draw the bottom panel
	  setOptions(ratios[k], true, true);
	  if(k==0) ratios[k]->Draw("PE");
	  else 
	     ratios[k]->Draw("samePE");
       }
   }

  pad2->SetGrid();
  
  pad2->Update();  
  pad2->RedrawAxis();  
  pad2->RedrawAxis("G");  

  pad1->Update();  
  pad1->RedrawAxis();  
  pad1->RedrawAxis("G");  

  gPad->RedrawAxis();
  gPad->RedrawAxis("G");


}

void histoCompare::drawLegend(TH1D *baseH, vector <TH1D*> list ){
	
	TLegend *leg = new TLegend(0.7,0.9,0.99,0.7);

	leg->AddEntry(baseH, names[0]);
	for(unsigned int i=0; i < list.size(); i++){ 
		leg->AddEntry(list[i], names[i+1]);
	}

	leg->Draw("same");


}

TString histoCompare::projectionInfo(){

   TString t = "Slice in ";
   if(projectionX) t.Append("Y between  ");
   else t.Append("X between  "); 

   t.Append(TString(Form("%1.1f - %1.1f ",projectionMin , projectionMax)));

   return t;

}




void histoCompare::setOptions(TH1D *h, bool dataLike, bool isBottom){

	h->GetXaxis()->SetTitle(titleX);
        h->GetYaxis()->SetTitle(titleY);
        if(dataLike) {
		h->SetMarkerColor(h->GetLineColor());
		h->SetMarkerStyle(20);
	}

        h->SetLineWidth(2);

	if(isBottom){
		h->GetYaxis()->SetRangeUser(0.5,1.5);
    		h->GetYaxis()->SetTitleSize(0.07);
      		h->GetYaxis()->SetLabelSize(0.07);
      		h->GetYaxis()->SetNdivisions(5, true);
        	h->GetXaxis()->SetTitleSize(0.07);
	  	h->GetXaxis()->SetLabelSize(0.07);
	}
}


void histoCompare:: project(){

     //clear projected histos
     projectedBase = NULL;
     projectedList.clear();

     // rebin
     base.Rebin2D(rebinX,rebinY);
     
     // finding the bins
     if(projectionX){
	   if(projectionMin !=UNDEFINED) binMin = base.GetYaxis()->FindBin(projectionMin);
	   if(projectionMax !=UNDEFINED) binMax = base.GetYaxis()->FindBin(projectionMax);
		}

       else{ 
	   if(projectionMin !=UNDEFINED) binMin = base.GetXaxis()->FindBin(projectionMin);
	   if(projectionMax !=UNDEFINED) binMax = base.GetXaxis()->FindBin(projectionMax);
	}


     // do the projection
     if(projectionX) 
	     projectedBase = base.ProjectionX("projectBaseX",binMin, binMax);
     else 
	     projectedBase = base.ProjectionY("projectBaseY",binMin, binMax);


     for(unsigned int i=0; i < compareList.size(); i++) {
	     compareList[i].Rebin2D(rebinX,rebinY); 
	     if(projectionX)
		     projectedList.push_back(compareList[i].ProjectionX(Form("projectCompX_%d",i),binMin, binMax));
	     else
		     projectedList.push_back(compareList[i].ProjectionY(Form("projectCompY_%d",i),binMin, binMax));
     }



}


void histoCompare::setNameofComponent(unsigned int i, TString n){

	if(i <= names.size() ) {
		names[i] = n;
	}
	else{
		cout <<"histoCompare::setNameofComponent - ERROR: no component for index " << i << endl;
	}
}





