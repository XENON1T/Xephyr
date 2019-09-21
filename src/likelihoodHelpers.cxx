#include "likelihoodHelpers.h"

using json = nlohmann::json;

namespace likelihoodHelpers  {

pdfComponent* genModel(json model_json) {
    json model_def;

    if (model_json.is_string()) {
        std::ifstream model_file(model_json.get<std::string>());
        model_file >> model_def;
    }
    else {
        model_def = model_json;
    }
    
    TString hist_path = model_def["file_path"].get<std::string>();
    if (hist_path.Index("/") != 0) hist_path = xephyrDir+hist_path;

    pdfComponent *model = new pdfComponent(model_def["name"].get<std::string>(), model_def["histogram_name"].get<std::string>(), hist_path);
    if (model_def["shape_parameters"].is_null()){
        model->autoLoad();
    } else {
        for (unsigned int i=0; i<model_def["shape_parameters"].size(); i++) {
            json p_def = model_def["shape_parameters"][i];
            shapeSys *p = new shapeSys(p_def["name"].get<std::string>());
            p->setStep(p_def["step_size"].get<double>());                      // FIXME
            p->setMinimum(p_def["lower_limit"].get<double>());
            p->setMaximum(p_def["upper_limit"].get<double>());
            p->setType(parameterTypeMap[p_def["type"].get<std::string>()]);
            model->addShapeSys(p);
        }

    }
    

    for (unsigned int i=0; i<model_def["rate_parameters"].size(); i++) {
        json p_def = model_def["rate_parameters"][i];
        scaleSys *scale = new scaleSys(p_def["name"].get<std::string>(),p_def["default_value"].get<double>());
        scale->setMinimum(p_def["lower_limit"].get<double>());
        scale->setMaximum(p_def["upper_limit"].get<double>());
        scale->setType(parameterTypeMap[p_def["type"].get<std::string>()]);
        model->addScaleSys(scale);
    }

    if (not model_def["suffix"].is_null()) 
        model->suffix = (TString) model_def["suffix"].get<std::string>();
    
    if (not model_def["exp_events"].is_null())
        model->setEvents(model_def["exp_events"].get<uint>() ); 

    if (not model_def["scale_factor"].is_null())
        model->setScaleFactor(model_def["scale_factor"].get<double>() ); 
    return model;
    }


dataHandler* genDataHandler(json dh_json) {
    json dh_def;
    if (dh_json.is_string()) {
        std::ifstream dh_file(dh_json.get<std::string>());
        dh_file >> dh_def;
    }
    else 
        dh_def = dh_json;

    TString hist_path = dh_def["histogram_file"].get<std::string>();
    if (hist_path.Index("/") != 0) hist_path = xephyrDir+hist_path;
    dataHandler *dh = new dataHandler(dh_def["name"].get<std::string>(), hist_path, 
                                         dh_def["histogram_name"].get<std::string>());
    if (not dh_def["prefix_treename"].is_null()){
    dh->setPrefixTree(dh_def["prefix_treename"].get<std::string>());
        }
    return dh;
    
}

TH2F* genAdditionalSafeGuardComponent(json comp_json){
    json comp_def;
    if (comp_json.is_string()) {
        std::ifstream comp_file(comp_json.get<std::string>());
        comp_file >> comp_def;
    }
    else 
        comp_def = comp_json;
    TString emtpy = "";
    TString hist_path = emtpy+comp_def["histogram_file"].get<std::string>();
    if (hist_path.Index("/") != 0) hist_path = xephyrDir + comp_def["histogram_file"].get<std::string>();
    TFile *f1 = TFile::Open(hist_path);
    TH2F* comp = (TH2F*) f1->Get( emtpy+comp_def["histogram_name"].get<std::string>() );

    for (uint i=0; i<comp_def["extra_histograms"].size(); i++) {
        json extra_def = comp_def["extra_histograms"][i];
        TString emtpy = "";
        hist_path = emtpy+extra_def["histogram_file"].get<std::string>();
        if (hist_path.Index("/") != 0) hist_path = xephyrDir + extra_def["histogram_file"].get<std::string>();
        TFile *fextra = TFile::Open(hist_path);
        TH2F* extra_histo = (TH2F*) fextra->Get(emtpy+extra_def["histogram_name"].get<std::string>() );
        comp->Add(extra_histo, extra_def["multiplier"].get<double>() );
    }
    
    comp->Scale(comp_def["scale"].get<double>());
    return comp;
}

pdfLikelihood* genLikelihood(json pl_json) {
    json pl_def;
    if (pl_json.is_string()) {
        std::ifstream pl_file(pl_json.get<std::string>());
        pl_file >> pl_def;
    }
    else 
        pl_def = pl_json;
    
    pdfLikelihood *pl = new pdfLikelihood(pl_def["name"].get<std::string>(), pl_def["mass"].get<double>());
    if (not pl_def["index"].is_null()) pl->setExperiment(pl_def["index"].get<uint>());
    
    for (unsigned int i=0; i<pl_def["models"].size(); i++) {
        json model_def = pl_def["models"][i];
        pdfComponent* model = genModel(model_def);
        if (model_def["type"].get<std::string>() == "BACKGROUND"){
            pl->addBkgPdfComponent(model, model_def["safeguard"].get<bool>());
        } else if (model_def["type"].get<std::string>() == "SIGNAL")
        {
            pl->setSignalPdf(model);
        }
    }
    
    pl->setSignalDefaultNorm(pl_def["signal_default_norm"].get<double>());
    pl->setWithSafeGuard(pl_def["safeguard"].get<bool>());
    pl->setSafeGuardPosDef(pl_def["posdef_safeguard"].get<bool>());
    
    for (unsigned int i=0; i<pl_def["datasets"].size(); i++) {
        json dh_def = pl_def["datasets"][i];
        dataHandler *dh = genDataHandler(dh_def);
        if (dh_def["type"].get<std::string>() == "DATA"){
            pl->setDataHandler(dh);
        } else if (dh_def["type"].get<std::string>() == "CALIBRATION") {
            pl->setCalibrationData(dh);
        }
    }
    
    if (not pl_def["additional_safeguard_component"].is_null()){
        TH2F *comp = genAdditionalSafeGuardComponent(pl_def["additional_safeguard_component"]);
        pl->setAdditionalSafeGuardComponent(comp);
    }
    
    pl->initialize();

    return pl;
  
}


CombinedProfileLikelihood* genCombinedLikelihood(json cpl_json) {

    json cpl_def;
    if (cpl_json.is_string()){
        std::ifstream cpl_file(cpl_json.get<std::string>());
        cpl_file >> cpl_def;
    } else 
        cpl_def = cpl_json;

    CombinedProfileLikelihood *cpl = new CombinedProfileLikelihood(cpl_def["name"].get<std::string>());
    cout<<"Combining "<<cpl_def["likelihoods"].size()<<" likelihoods"<<endl;
    if (cpl_def["likelihoods"].size()<1) return cpl;
    TString corr_prefix = "CORR_";
    cout<<"Starting to generate individual likelihoods"<<endl;
    vector <pdfLikelihood*> pls;
    for (uint i=0; i<cpl_def["likelihoods"].size(); i++){
        json pl_def = cpl_def["likelihoods"][i];
        cout<<"Generating likelihood"<<i<<endl;
        pdfLikelihood *pl = genLikelihood(pl_def);
        pl->setExperiment(i + 1);
        pl->setPrintLevel(INFO);
        // pl->initialize();
        pls.push_back(pl);
    }
    cout<<"Done generating likelihoods"<<endl;
    cout<<"Correlating parameters"<<endl;

    for (uint i=0; i<cpl_def["combined_parameters"].size(); i++){
        json cp_def = cpl_def["combined_parameters"][i];
        CombinedParameter *cp = new CombinedParameter(cp_def["name"].get<std::string>());
        cout<<"Correlating "<<cp_def["name"].get<std::string>()<<endl;
        for (uint k=0; k<cp_def["parameters"].size(); k++) {
            json p_def = cp_def["parameters"][k];
            pdfLikelihood* pl = pls[p_def["likelihood"].get<uint>()];
            pdfComponent *model;
            if (cp_def["model_type"].get<std::string>()=="BACKGROUND") {
                model = pl->getBkgComponent(p_def["model_name"].get<std::string>());
            } else if(cp_def["model_type"].get<std::string>()=="SIGNAL") {
                model = pl->signal_component;
            } else {
                // cout<<"ERROR: UNRECOGNIZED MODEL TYPE "<<cp_def["model_type"].get<std::string>()<<endl;
                Error("genCombinedLikelihood", "UNRECOGNIZED MODEL TYPE");
            }

            if (cp_def["type"].get<std::string>()=="SHAPE") {
                cp->correlateParameter(model->getShapeSys(p_def["name"].get<std::string>()));
            } else if (cp_def["type"].get<std::string>()=="RATE") {
                cp->correlateParameter(model->getScaleSys(p_def["name"].get<std::string>()));
            } else {
                Error("genCombinedLikelihood", "UNRECOGNIZED PARAMETER TYPE");
                // cout<<"ERROR: UNRECOGNIZED PARAMETER TYPE "<<cp_def["sys_type"].get<std::string>()<<endl;
            }
            }
            cp->setMinimum(cp_def["lower_limit"].get<double>());
            cp->setMaximum(cp_def["upper_limit"].get<double>());
            if (not cp_def["parameter_type"].is_null()){
            cp->setType(parameterTypeMap[cp_def["parameter_type"].get<std::string>()]);
            }
            cpl->addParameter(cp, AUTO);
        
        }

    // combine!
    cout<<"Starting to combine individual likelihoods"<<endl;
    for(unsigned int i=0 ; i < pls.size(); i++) cpl->combine(pls[i]);

    cpl->initialize();
    cpl->setPrintLevel(INFO);
    cout<<"Done."<<endl;
    return cpl;

}


}