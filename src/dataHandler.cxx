#include "dataHandler.h"

dataHandler::dataHandler(TString name) : errorHandler("dataHandler"), Name(name){

	DMdata = NULL;
	file = NULL;
	sumOfWeights=0;
	gs1s2w=0;
	
	
	s1 = 0.;
	s2 = 0.;
	weight = 1.;

	dataType = UNSPECIFIED_DATA; 

	FirstVarName   = "cs1";  //default var in data
	SecondVarName  = "cs2";  //default var name
}

dataHandler::dataHandler(TString name, TH2F *h2pdf, int N) : errorHandler("dataHandler"), Name(name){

      	DMdata = NULL;


	FirstVarName   = "cs1";  //default var in data
	SecondVarName  = "cs2";  //default var name
	
	s1 = 0.;
	s2 = 0.;
	weight = 1.;
	gs1s2w=new TGraph2D();
	sumOfWeights=0;
	generateDataSet(h2pdf,N);
}


dataHandler::dataHandler(TString name, TH2F *h2pdf) : errorHandler("dataHandler"), Name(name){

  DMdata = NULL;

  FirstVarName   = "cs1";  //default var in data
  SecondVarName  = "cs2";  //default var name

  sumOfWeights=0;
  s1 = 0.;
  s2 = 0.;
  weight = 1.;
  
  gs1s2w=new TGraph2D();
  dataType = ASIMOV_DATA; 

  generateAsimov(h2pdf);

  
}

dataHandler::dataHandler(TString name, TString fileName, TString dmTree) : errorHandler("dataHandler"), Name(name){

	file = TFile::Open(fileName);

	if(file == NULL)
		Error("dataHandler","file " + fileName + " does not exist. Quit.");

	if( file->FindKey(dmTree) == NULL) 
		Error("dataHandler","TTree " + dmTree+ " does not exist in file " 
			                       + fileName + ". Quit.");


	FirstVarName   = "cs1";  //default var in data
	SecondVarName  = "cs2";  //default var name


	DMdata = (TTree*) file->Get(dmTree);
	
	s1 = 0.;
	s2 = 0.;
	weight = 1.;
	sumOfWeights=0;
	DMdata->SetBranchAddress(FirstVarName,&s1);
	DMdata->SetBranchAddress(SecondVarName,&s2);
	
	dataType = DM_DATA;

	gs1s2w=new TGraph2D();
	for (int i=0; i< DMdata->GetEntries(); i++) {
	   weight = 1.;
	   DMdata->GetEntry(i);
	   gs1s2w->SetPoint(gs1s2w->GetN(),s1,s2,weight);
	   sumOfWeights+=weight;
	}


}



dataHandler::~dataHandler(){

	delete DMdata;
	delete gs1s2w;

}


Long64_t dataHandler::getEntries(){
  return gs1s2w->GetN(); }

double dataHandler::getSumOfWeights(){
  return sumOfWeights; }

void dataHandler::getEntry(Long64_t entry) {
  if(DMdata ==NULL)  Error("getEntry","No data is set.");
  if(entry > getEntries() )  Error("getEntry"," Entry number outside range");  
  s1=gs1s2w->GetX()[entry];
  s2=gs1s2w->GetY()[entry];
  weight=gs1s2w->GetZ()[entry];
	
}

void dataHandler::generateAsimov( TH2F *background ){
    delete DMdata;
	delete gs1s2w;
	double ts1 = 0., ts2 =0., tw = 0.;
	DMdata = new TTree(Name,"FakeData "+Name);
	DMdata->Branch(FirstVarName,&ts1);
	DMdata->Branch(SecondVarName,&ts2);
	DMdata->Branch("weight", &tw); 
	gs1s2w=new TGraph2D();
    sumOfWeights=0;

	for (int ix=1; ix<=background->GetNbinsX(); ix++) {
	  for (int iy=1; iy<=background->GetNbinsY(); iy++) {
	   ts1 = background->GetXaxis()->GetBinCenter(ix);
	   ts2 = background->GetYaxis()->GetBinCenter(iy);
	   tw=background->GetBinContent(ix,iy);
	   DMdata->Fill();
	   gs1s2w->SetPoint(gs1s2w->GetN(),ts1,ts2,tw);
	   sumOfWeights+=tw;
	  }
	}
}


double dataHandler::integrate(TH2F *histo, double s1_min, double s1_max, double s2_min, double s2_max){
	
	
	 TAxis *xaxis = histo->GetXaxis();
	 TAxis *yaxis = histo->GetYaxis();

	  int xmin = xaxis->FindBin(s1_min);
	  int xmax = xaxis->FindBin(s1_max);
	  int ymin = yaxis->FindBin(s2_min);
	  int ymax = yaxis->FindBin(s2_max);
	  
	  
	  double integral = 0.;

 	 //for each slice in Y we remove the excess in X, since TH2F::Integral() will always over estimate the integral
	  for(int y = ymin; y <= ymax ; y++){
		
		double temp_integral = 0.;
		temp_integral += histo->Integral(xmin,xmax, y, y);

		// removing x excess		
	        temp_integral -= histo->GetBinContent(xmin, y)*( s1_min - xaxis->GetBinLowEdge(xmin)) / xaxis->GetBinWidth(xmin);
		temp_integral -= histo->GetBinContent(xmax, y)*(xaxis->GetBinUpEdge(xmax)-s1_max)/ xaxis->GetBinWidth(xmax);

		//removing y excess
		if(y == ymin)  
			temp_integral -= temp_integral * ( s2_min - yaxis->GetBinLowEdge(ymin) ) / 
					yaxis->GetBinWidth(ymin);
		if(y == ymax) 
			temp_integral -= temp_integral * ( yaxis->GetBinUpEdge(ymax) - s2_max ) / 
					yaxis->GetBinWidth(ymax);

		integral += temp_integral;

	  }

	return integral;

}






double dataHandler::getValFromPdf( TH2F &histo ) {

	//s1 and s2 are the current values taken from the relative
	//branches of either DMdata or asimovData 
	int s1_bin = histo.GetXaxis()->FindBin(s1);
	int s2_bin = histo.GetYaxis()->FindBin(s2);

	return histo.GetBinContent(s1_bin,s2_bin);
}




void dataHandler::fillDataHisto(TH2F *hist){

	if(hist == NULL)
		Error("fillDataHisto", "you passed me a NULL pointer, quit.");


	for(Long64_t entry =0; entry < getEntries(); entry++){
		getEntry(entry);
		hist->Fill(s1,s2);
	}



}


void dataHandler::setFileAndTree(TString PathtoFile, TString nameTree){
	setFile(PathtoFile);
	setDataTree(nameTree);

}
void dataHandler::setFile(TString PathtoFile){

	file = TFile::Open(PathtoFile);

	if(file == NULL)
		Error("dataHandler","file " + PathtoFile + " does not exist. Quit.");

}


void dataHandler::setDataTree(TString nameTree){

	if( file->FindKey(nameTree) == NULL) 
		Error("dataHandler","TTree " + nameTree+ " does not exist in file " 
			                       + file->GetName() + ". Quit.");

	
	setDataTree( (TTree*) file->Get(nameTree) );
	
}

void dataHandler::setDataTree(TTree *tree){

	delete DMdata;
	delete gs1s2w;
	
	gs1s2w=new TGraph2D();
	DMdata = tree;
	
	DMdata->SetBranchAddress(FirstVarName,&s1);
	DMdata->SetBranchAddress(SecondVarName,&s2);
	
	dataType = DM_DATA;
	sumOfWeights=0;
	for (Long64_t i=0; i< DMdata->GetEntries(); i++) {
	   weight = 1.;
	   DMdata->GetEntry(i);
	   gs1s2w->SetPoint(gs1s2w->GetN(),s1,s2,weight);
	   sumOfWeights+=weight;
	}
}

void dataHandler::generateDataSet(TH2F *h2pdf, int N){
  delete DMdata;
  delete gs1s2w;
  sumOfWeights=0;
  gs1s2w=new TGraph2D();
  DMdata = new TNtuple(Name,"FakeData "+Name,FirstVarName+":"+SecondVarName+":weight"); 
  dataType = DM_SIMULATED_DATA;
  Name="Fake data set:";
  addToDataSet(h2pdf,N);
}


void dataHandler::addToDataSet(TH2F *h2pdf, int N){

  dataType = DM_SIMULATED_DATA;
  if(dataType != DM_SIMULATED_DATA)
    Error("addTodataSet","Cannot add fake data to this data set. Use generateDataSet first");
  double ts1,ts2,tw = 0.;
  DMdata->SetBranchAddress(FirstVarName,&ts1);
  DMdata->SetBranchAddress(SecondVarName,&ts2);
  Name+=Form("%s(%d),",h2pdf->GetName(),N);
  gRandom->SetSeed(0);
  for (int i=0; i<N; i++) {
    h2pdf->GetRandom2(ts1,ts2);
    tw=1;
    DMdata->Fill();
    gs1s2w->SetPoint(gs1s2w->GetN(),ts1,ts2,tw);
    sumOfWeights+=tw;
  }
  

}


void dataHandler::drawS1S2(TString opt) {

  if(DMdata != NULL){
    DMdata->Draw(FirstVarName+":"+SecondVarName,"","goff");
    TGraph *gr=new TGraph(DMdata->GetSelectedRows(),
			  DMdata->GetV1(), DMdata->GetV2());
    gr->SetTitle(Name+";"+FirstVarName+";"+SecondVarName);
    gr->Draw(opt);
  }


}

TGraph dataHandler::getS1S2() {
  TGraph gr;
 if (DMdata != NULL){
   gr = TGraph(gs1s2w->GetN(), gs1s2w->GetX(),gs1s2w->GetY());
   }

return gr;
}


void dataHandler::printSummary() {

  printf ("dataHandler:: summary:  name= %s,  N=%d \n", Name.Data(),gs1s2w->GetN());
  return; 

}



vector<int> dataHandler::getSimulatedInfo(unsigned int size){

	vector<int> out;
	for(unsigned int k=0; k <size; k++) out.push_back(0.);

	
	float type = 0.;

	DMdata->SetBranchAddress("type",&type);
	
	for (Long64_t i=0; i< DMdata->GetEntries(); i++) {

     	   DMdata->GetEntry(i);

	   for(unsigned int j=0; j< size; j++) {
	   	if(type == j) out[j] = out[j] + 1;
	   }
	}


	return out;
}



