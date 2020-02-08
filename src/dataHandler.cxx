#include "dataHandler.h"


void dataHandler::initialize(){
	DMdata = NULL;
	file = NULL;
	sumOfWeights = 0.;
	weights = NULL;

	x = 0.;
	y = 0.;
	z = 0.;
	weight = 1.;

	dataType = UNSPECIFIED_DATA; 

	x_name   = "x";  //default var in data
	y_name  = "y";  //default var name
	z_name  = "z";  //default var name
}

dataHandler::dataHandler(TString name) : errorHandler("dataHandler"), Name(name){

	DMdata = NULL;
	file = NULL;
	sumOfWeights = 0.;
	weights = NULL;

	x = 0.;
	y = 0.;
	z = 0.;
	weight = 1.;
	x_name   = "x";  //default var in data
	y_name  = "y";  //default var name
	z_name  = "z";  //default var name
	weights = new TNtuple("Weights", "Data weights", x_name+":"+y_name+":"+z_name+":weight");

	dataType = UNSPECIFIED_DATA; 



}

dataHandler::dataHandler(TString name, TH3F *h3pdf, int N) : errorHandler("dataHandler"), Name(name){

      	DMdata = NULL;


	x_name   = "x";  //default var in data
	y_name  = "y";  //default var name
	z_name  = "z";

	x = 0.;
	y = 0.;
	z = 0.;
	weight = 1.;
	weights = new TNtuple("Weights", "Data weights", x_name+":"+y_name+":"+z_name+":weight");
	sumOfWeights = 0;
	generateDataSet(h3pdf,N);
}


dataHandler::dataHandler(TString name, TH3F *h3pdf) : errorHandler("dataHandler"), Name(name){

  	DMdata = NULL;

	x_name   = "x";  //default var in data
	y_name  = "y";  //default var name
	z_name  = "z";
	sumOfWeights=0;
	x = 0.;
	y = 0.;
	weight = 1.;
  
	weights = new TNtuple("Weights", "Data weights", x_name+":"+y_name+":"+z_name+":weight");
	dataType = ASIMOV_DATA; 

	generateAsimov(h3pdf);

}

dataHandler::dataHandler(TString name, TString fileName, TString dmTree, 
						TString x_name, TString y_name, TString z_name) : errorHandler("dataHandler"), 
						Name(name), x_name(x_name), y_name(y_name), z_name(z_name){

	file = TFile::Open(fileName);

	if(file == NULL)
		Error("dataHandler","file " + fileName + " does not exist. Quit.");

	if( file->FindKey(dmTree) == NULL) 
		Error("dataHandler","TTree " + dmTree+ " does not exist in file " 
			                       + fileName + ". Quit.");

	// x_name  = x_name;  //default var in data
	// y_name  = y_name;  //default var name
	// z_name  = z_name;

	DMdata = (TTree*) file->Get(dmTree);
	
	x = 0.;
	y = 0.;
	z = 0.;
	weight = 1.;
	sumOfWeights=0;
	
	DMdata->SetBranchAddress(x_name,&x);
	DMdata->SetBranchAddress(y_name,&y);
	DMdata->SetBranchAddress(y_name,&z);
	dataType = DM_DATA;

	weights = new TNtuple("Weights", "Data weights",x_name+":"+y_name+":"+z_name+":weight");

	for (int i=0; i< DMdata->GetEntries(); i++) {
	   weight = 1.;
	   DMdata->GetEntry(i);
	   weights->Fill(x,y,z,weight);
	   sumOfWeights+=weight;
	}


}



dataHandler::~dataHandler(){

	delete DMdata;
	delete weights;

}


Long64_t dataHandler::getEntries(){
  return weights->GetEntries(); }

Float_t dataHandler::getSumOfWeights(){
  return sumOfWeights; }

void dataHandler::getEntry(Long64_t entry) {
  if(DMdata==NULL)  Error("getEntry","No data is set.");
  if(entry > getEntries() )  Error("getEntry"," Entry number outside range");  
  weights->GetEntry(entry);
//   x=GetX()[entry];
//   y=GetY()[entry];
//   z=GetZ()[entry];
//   weight=GetW()[entry];
}

void dataHandler::generateAsimov( TH3F *background ){
    delete DMdata;
	delete weights;
	Float_t tx = 0., ty =0., tz =0., tw = 0.;
	DMdata = new TTree(Name,"FakeData "+Name);
	DMdata->Branch(x_name,&tx);
	DMdata->Branch(y_name,&ty);
	DMdata->Branch(z_name,&tz);
	DMdata->Branch("weight", &tw); 
	weights = new TNtuple("Weights", "Data weights", x_name+":"+y_name+":"+z_name+":weight");
    sumOfWeights=0;

	for (int ix=1; ix<=background->GetNbinsX(); ix++) {
	  for (int iy=1; iy<=background->GetNbinsY(); iy++) {
		  for (int iz=1; iz<=background->GetNbinsZ(); iz++) {
	   tx = background->GetXaxis()->GetBinCenter(ix);
	   ty = background->GetYaxis()->GetBinCenter(iy);
	   tz = background->GetZaxis()->GetBinCenter(iz);
	   tw=background->GetBinContent(ix,iy,iz);
	   DMdata->Fill();
	   weights->Fill(tx,ty,tz,tw);
	   sumOfWeights+=tw;
		  }
	  }
	}
}


Float_t dataHandler::integrate(TH3F *hist, Float_t xmin, Float_t xmax, 
									Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax){
    	
	TAxis *xaxis = hist->GetXaxis();
	TAxis *yaxis = hist->GetYaxis();
	TAxis *zaxis = hist->GetZaxis();

    int bxmin = xaxis->FindBin(xmin);
    int bxmax = xaxis->FindBin(xmax);
    int bymin = yaxis->FindBin(ymin);
    int bymax = yaxis->FindBin(ymax);
    int bzmin = zaxis->FindBin(zmin);
    int bzmax = zaxis->FindBin(zmax);
	  
	Float_t integral = 0.;
    Float_t dx, dy, dz, dV = 0;

    for (int x=bxmin; x<=bxmax; x++){
        for (int y=bymin; y<=bymax; y++){
            for (int z=bzmin; z<=bzmax; z++){
                if (x==bxmin){
                    dx = abs(xaxis->GetBinUpEdge(x) - xmin)/xaxis->GetBinWidth(x);
                } else if (x==bxmax){
                    dx = abs(xmax - xaxis->GetBinLowEdge(x))/xaxis->GetBinWidth(x);
                } else {
                    dx = 1.; 
                }

                if (y==bymin){
                    dy = abs(yaxis->GetBinUpEdge(y) - ymin)/yaxis->GetBinWidth(y);
                } else if (y==bymax){
                    dy = abs(ymax - yaxis->GetBinLowEdge(y))/yaxis->GetBinWidth(y);
                } else {
                    dy = 1.; 
                }

                if (z==bzmin){
                    dz = abs(zaxis->GetBinUpEdge(z) - zmin)/zaxis->GetBinWidth(z);
                } else if (z==bzmax){
                    dz = abs(zmax - zaxis->GetBinLowEdge(z))/zaxis->GetBinWidth(z);
                } else {
                    dz = 1.; 
                }

                // dV = xaxis->GetBinWidth(x)*yaxis->GetBinWidth(y)*zaxis->GetBinWidth(z);
                integral += hist->GetBinContent(x,y,z)*dx*dy*dz;

             }

        }

    }

    return integral;

}


Float_t dataHandler::getValFromPdf( TH3F &histo ) {

	//x and y are the current values taken from the relative
	//branches of either DMdata or asimovData 
	int x_bin = histo.GetXaxis()->FindBin(x);
	int y_bin = histo.GetYaxis()->FindBin(y);
	int z_bin = histo.GetZaxis()->FindBin(z);
	return histo.GetBinContent(x_bin,y_bin,z_bin);
}




void dataHandler::fillDataHisto(TH3F *hist){

	if(hist == NULL)
		Error("fillDataHisto", "you passed me a NULL pointer, quit.");


	for(Long64_t entry =0; entry < getEntries(); entry++){
		getEntry(entry);
		hist->Fill(x,y,z);
	}

}

vector<double> dataHandler::getTrueParams(){

    // retrive the previously saved TList of parameters (done in ToyGenerator)
    TIter iterateMe(DMdata->GetUserInfo());

    TParameter<double> *parameter = NULL;
	vector <double> true_params;
	true_params.clear();

	// if not MC toys then does not have truth so return empty
	if(DMdata->GetUserInfo()->IsEmpty()) return true_params;

    // saving the parameters of the tree
    while ((parameter = (TParameter<double>*)iterateMe())) {
        true_params.push_back( parameter->GetVal() );
    }
	return true_params;
}

vector<string> dataHandler::getTrueParamsNames(){

    // retrive the previously saved TList of parameters (done in ToyGenerator)
    TIter iterateMe(DMdata->GetUserInfo());

    TParameter<double> *parameter = NULL;
	vector <string> true_params;
	true_params.clear();
	
	// if not MC toys then does not have truth so return empty
	if(DMdata->GetUserInfo()->IsEmpty()) return true_params;

    // saving the parameters of the tree
    while ((parameter = (TParameter<double>*)iterateMe())) {
        true_params.push_back( parameter->GetName() );
    }

	return true_params;
}

void dataHandler::setTreeIndex( int index ){
	if( TreePrefix == "" ) Error( "setTreeIndex", "Tree Prefix not set, use setPrefixTree().");
	
	TString newTree =  TreePrefix + "_" +TString::Itoa(index, 10);
	Debug("setTreeIndex", "setting new tree: " + newTree);
	setDataTree(newTree);
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

	//delete DMdata;  // no much reason to delete this, since adding from file or existing tree
	delete weights;

	// changing name to the data handler
	Name = TString("Data_") + tree->GetName();
	
	weights = new TNtuple("Weights", "Data weights", x_name+":"+y_name+":"+z_name+":weight");
	DMdata = tree;
	float tx,ty,tz = 0.;
	DMdata->SetBranchAddress(x_name,&tx);
	DMdata->SetBranchAddress(y_name,&ty);
	DMdata->SetBranchAddress(z_name,&tz);
	dataType = DM_DATA;
	sumOfWeights=0;
	weight = 1.;
	for (Long64_t i=0; i< DMdata->GetEntries(); i++) {
	   DMdata->GetEntry(i);
	   weights->Fill(tx,ty,tz,weight);
	   sumOfWeights+=weight;
	}
}

void dataHandler::generateDataSet(TH3F *h3pdf, int N){
  delete DMdata;
  delete weights;
  sumOfWeights=0;
  weights = new TNtuple("Weights", "Data weights", x_name+":"+y_name+":"+z_name+":weight");
  DMdata = new TNtuple(Name,"FakeData "+Name,x_name+":"+y_name+":"+z_name+":weight"); 
  dataType = DM_SIMULATED_DATA;
  Name="Fake data set:";
  addToDataSet(h3pdf, N);
}


void dataHandler::addToDataSet(TH3F *h3pdf, int N){

  dataType = DM_SIMULATED_DATA;
  if(dataType != DM_SIMULATED_DATA)
    Error("addTodataSet","Cannot add fake data to this data set. Use generateDataSet first");
  Float_t tx,ty,tz = 0.;
  Double_t x,y,z = 0.;
  Float_t tw = 1;
  DMdata->SetBranchAddress(x_name,  &tx);
  DMdata->SetBranchAddress(y_name, &ty);
  DMdata->SetBranchAddress(z_name,  &tz);
  Name+=Form("%s(%d),",h3pdf->GetName(),N);
  gRandom->SetSeed(0);
  
  for (int i=0; i<N; i++) {
    h3pdf->GetRandom3(x,y,z);
    DMdata->Fill();
	weights->Fill(x,y,z,tw);
    sumOfWeights+=tw;
  }
  
}


void dataHandler::drawxy(TString opt) {

  if(DMdata != NULL){
    DMdata->Draw(x_name+":"+y_name,"","goff");
    TGraph *gr=new TGraph(DMdata->GetSelectedRows(),
			  DMdata->GetV1(), DMdata->GetV2());
    gr->SetTitle(Name+";"+x_name+";"+y_name);
    gr->Draw(opt);
  }


}

TGraph dataHandler::getXY() {
  TGraph gr;
//  if (DMdata != NULL){
//    gr = TGraph(weights->GetN(), weights->GetX(),weights->GetY());
//    }

return gr;
}


void dataHandler::printSummary() {

  printf ("dataHandler:: summary:  name= %s,  N=%lld \n", Name.Data(),weights->GetEntries());
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

Float_t dataHandler::getX(long N) {
	 if (N>weights->GetEntries()) {
		printf ("ERROR %ld larger than Entries (%lld) \n",N,weights->GetEntries());
		return 0;
		} else {
			Float_t x;
			weights->SetBranchAddress(x_name,&x);
			weights->GetEntry(N); 
			return x;
		}
	}


Float_t dataHandler::getY(long N) {
	 if (N>weights->GetEntries()) {
		printf ("ERROR %ld larger than Entries (%lld) \n",N,weights->GetEntries());
		return 0;
		} else {
			Float_t y;
			weights->SetBranchAddress(y_name,&y);
			weights->GetEntry(N); 
			return y;
		}
	}


Float_t dataHandler::getZ(long N) {
	 if (N>weights->GetEntries()) {
		printf ("ERROR %ld larger than Entries (%lld) \n",N,weights->GetEntries());
		return 0;
		} else {
			Float_t z;
			weights->SetBranchAddress(z_name,&z);
			weights->GetEntry(N); 
			return z;
		}
	}


Float_t dataHandler::getW(long N) {
	 if (N>weights->GetEntries()) {
		printf ("ERROR %ld larger than Entries (%lld) \n",N,weights->GetEntries());
		return 0;
		} else {
			Float_t weight;
			weights->SetBranchAddress("weight",&weight);
			weights->GetEntry(N); 
			return weight;
		}
	}


