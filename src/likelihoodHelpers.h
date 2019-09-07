#include "XeLikelihoods.h"
#include "XePdfObjects.h"
#include "XeStat.h"
#include "XeUtils.h"
#include "ToyFitterExclusion.h"
#include "ToyGenerator.h"
#include "dataHandler.h"
#include "TSystem.h"
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

///---------------     PARAMETERS ------------         ///

#ifndef LIKELIHOOD_HELPERS
#define LIKELIHOOD_HELPERS


namespace likelihoodHelpers  {
    const double year         = 365. ;   // days
    const TString xephyrDir(gSystem->Getenv("XEPHYR_DIR"));

    TH2F* genAdditionalSafeGuardComponent(json comp_json);
    dataHandler* genDataHandler(json dh_json);
    pdfComponent* genModel(json model_json);
    pdfLikelihood* genLikelihood(json pl_json);
    CombinedProfileLikelihood* genCombinedLikelihood(json cpl_json);

    map <std::string, parameterType> parameterTypeMap = {
        {"PARAMETER_OF_INTEREST", PARAMETER_OF_INTEREST },
        {"NUISANCE_PARAMETER", NUISANCE_PARAMETER },
        {"FIXED_PARAMETER", FIXED_PARAMETER },
        {"FROZEN_PARAMETER", FROZEN_PARAMETER },
        {"FREE_PARAMETER", FREE_PARAMETER },
        {"N_PARAMETER_TYPES", N_PARAMETER_TYPES }
    };

}

#endif