#ifndef TOY_GENERATOR
#define TOY_GENERATOR


#include "TRandom.h"
#include "TRint.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "XeUtils.h"
#include "XeLikelihoods.h"
#include "TParameter.h"
#include "XeStat.h"
#include "TTree.h"
#include "TH2F.h"
#include <map>
#include <vector>
#include <stdio.h>

using namespace std;


/**
 * \class ToyGenerator 
 * \brief Helper class to generate toys given an input likelihood. 
 * It can generate toys for the calibration dataset and for the "science data" dataset.
*/
class ToyGenerator: public errorHandler {

    public:
        //! Constructor: @param sampleName: name of the file and tree prefix
        //! @param outDir: path to dir where you want to save toys
        ToyGenerator(TString sampleName, TString outDir);

        //! \brief pointer to the defined likelihood
        void setLikelihood(pdfLikelihood *like) {likeHood = like;};

        //! \brief the total number of events to which the calibtration pdfs will be rescaled.
        //
        //! The fraction of events given to each pdf component is defined in the likelihood. 
        //! This is just an overall scale.
        void setAverageCalibrationEvents(double evnt) {averageCalEvnt = evnt;} ;


        //! \brief the total number of events to which the science data pdfs will be rescaled.
        //
        //! This overrides the number of events set in the likelihood.
        //! The fraction of events given to each pdf component is defined in the likelihood. 
        //! This is just an overall scale. Typically you don't want to set this, but you must 
        //! set "setAverageCalibrationEvents" instead.
        void setAverageDataEvents(double evnt) {averageDataEvnt = evnt;} ;

        //! \brief generate N toys of the calibration dataset.
        //
        // generate calibration like toys, where averageCalEvnt is the poisson median
        // of the number of generated event per toy. N is the number of toy datasets.
        // The pdfs are taken only from pdfcomponent that are "safeguarded" in the 
        // likelihood, their relative rates are respected. Also the additional 
        // contribution from "AdditionalSafeGuardComponent" is included.
        void generateCalibration(int N, bool randomizeNP = false);

        //! \brief generate N toys of the 'science data' dataset with injected signal fraction 'mu'.
        void generateData(double mu, int N, bool randomizeNP = false);

        //! set the seed for the toy generation. YOU MUST CHANGE THIS for any run.
        void setSeed(int seed);
        
        //! set the Generation, this info will be available in the generated tree (so that u can hadd them), it is optional and non ncecessary.
        void setGeneration(int generation) { Gen = generation; };

        //! set the LikeliHood Type, this info will be available in the generated tree (so that u can hadd them), it is optional and non ncecessary.
        void setLikeType(int lktype) {likelihoodType = lktype; } ;

        //! \brief this will randomize the nominal value of the nuissance parameters, according to their distro.
        //
        //! The generated toys will not be generated now from nominal distro. YOU MUST change seed
        //! for each repetition of this.
        void randomizeNuissanceParameter();

        //! \brief sets a new tree name prefix
        void setTreeName(TString newname) { treeName = newname;};


    private:
        
        //! \biref saves the current values of NP to the UserInfo TList of the Tree
        void saveParameters(TTree *f);
        
        //! \brief return the sum of default integral of all bkg components.
        //! 
        //! NOTE: is always the 'default' integral for nominal t-value of NP, not the current. 
        //! scale factors are applied but not scale SYS.
        double getModelIntegral();

        //! \brief return the sum of default integral of all bkg safeguarded components.
        //!
        //! NOTE: is always the 'default' integral for nominal t-value of NP, not the current. 
        //! scale factors are applied but not scale SYS.
        double getModelIntegralSafeguarded();

        //! \brief returns a vector of sys interpolated TH2F of each bkg
        vector<TH2F> getTH2OfBkg();


        double    averageCalEvnt;
        double    averageDataEvnt;
        int       Gen;
        int       likelihoodType;
        TRandom3  rambo;
        TString   treeName;
        TString   dir;
        pdfLikelihood *likeHood;
};




#endif