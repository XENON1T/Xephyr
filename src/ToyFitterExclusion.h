#ifndef TOY_FITTER
#define TOY_FITTER

#include "XeLikelihoods.h"
#include "XeUtils.h"
#include "dataHandler.h"
#include "TFile.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TRint.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TParameter.h"
#include "XeStat.h"
#include "TTree.h"
#include "TH2F.h"
#include <map>
#include <vector>
#include <stdio.h>

using namespace std;

/**
 * \class ToyFitterExclusion
 * \brief Class to handle the fitting of toy dataset for limit production.
 * 
 *  
 * - Produces a Tree with all the post fit values.
 * - Produces distributions of f(mu=mu_true | H_mu) and the 90% quantiles for a set of mu values.
 * - Computes sensitivity + bands 
 * - Given a science dataset it compute the limit.
 */
class ToyFitterExclusion: public errorHandler {

  public:
    
    //! \brief constructor, pathToFile is the name of the file including path to fit, 
    //! by convenction this ia also the name of the file.
    ToyFitterExclusion(TString pathToFile);
    
    /**
     * \brief returns the Test Statistic q_tilde(mu) for a value of mu_test.
     * 
     * It is a reproduction of q_tilde test statistic from equation 16 of: https://arxiv.org/abs/1007.1727
     * this is suitable for upper limit in which the parameter of interest must be >0.
     * For it to work correctly is important to set up the likelihood beforehand.
     * @param mu: the signal strenght for conditional fit. 
     */
    double computeTS(double mu);

    /**
     * \brief given as input a file with many data toys, it writes a tree with fit outputs
     * 
     * Note that it expects a toy file generated with naming convention of ToyGenerator.
     * The output tree contains: conditional and unconditional nuisance parameter 
     * post fit values, and the test statistic. This is supposed to be used for the alternative 
     * hypotesis fits, to be able to generate aftewards the graph of 90th percentiles.
     * @parmam mu: the signal strenght of the conditional fit
     * @param nameTree: the name prefix of the tree inside the file
     */
    void fit(double mu, TString nameTree);

    //! \brief set the likelihood to fit
    void setTheLikelihood(pdfLikelihood *like) { likeHood = like; };

    //! \brief set the path to dir in which the fits are stored otherwise is current
    void setOutputDir(TString path){ OutDir = path; };



  private:


    //! \brief save the postfit parameters into params.
    void saveParameters(double *params);

    //! \brief fills the initial paramteres value of true parameter an their name
    void fillTrueParams( TTree *inputTree);

    //! \brief assign to each parameter a random measured t-value
    void measureParameters();

    void saveNames(string *names);

    pdfLikelihood *likeHood;
    TString       dirPath;
    TString       treeName;
    TString       OutDir;               //! output dir if not set is current dir
    
    double        likelihood_uncond;    //! value of unconditional fit
    double        likelihood_cond;      //! value of conditional fit
    double        true_params[50]  ;    //! t_value of truth generated param
    double        measured_params[50];  //! t_value of measured param
    double        uncond_params[50] ;   //! unconditional fit
    double        cond_params[50] ;     //! conditional fit
    vector<string> name_params;         //! names of parameters
    TRandom3      rambo;                //! random handler
};

#endif