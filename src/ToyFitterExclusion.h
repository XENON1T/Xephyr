#ifndef TOY_FITTER
#define TOY_FITTER

#include "XeLikelihoods.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
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
     * @parmam mu: the signal strenght of the conditional fit.
     * @param nameTree: the name prefix of the tree inside the file.
     * @param randomizeMeasure: if it is going to randomize the measure of the parameter for each toy.
     * @parma soptAt:  optional, the number of toy in file you want to fit.
     */
    void fit(double mu, TString nameTree,  bool randomizeMeasure =true, int stopAt=-999);

    //! \brief set the likelihood to fit
    void setTheLikelihood(pdfLikelihood *like) { likeHood = like; };

    //! \brief set the path to dir in which the fits are stored otherwise is current
    void setOutputDir(TString path){ OutDir = path; };

    //! \brief produces a TGraph with the 90% quantiles of the alternate Hypothesis.
    //
    //! it is used to compute limits. The TH1F of the test statistic f(q_mu | H_mu) are
    //! stored in a root file toghether with the TGraph of 90% quantiles.
    //! @param fileName: takes as input the hadded post fit output from fit() of this class.
    //! @param mu_list: is the list of true mu hypothesis present in the file tree.
    //! @param mu_size: is the size of the previous list.
    TGraphAsymmErrors computeTSDistros(TString fileName, double *mu_list, int mu_size);


    //! \brief produces a TGraph with the 90% quantiles of the alternate Hypothesis.
    //
    //! it is used to compute limits. The TH1F of the test statistic f(q_mu | H_mu) are
    //! stored in a root file toghether with the TGraph of 90% quantiles.
    //! @param fileName: takes as input the hadded post fit output from fit() of this class.
    //! @param mu_list: is the list of true mu hypothesis present in the file tree.
    //! @param mu_size: is the size of the previous list.
    TGraphAsymmErrors computeTSDistros(TTree *tree, double *mu_list, int mu_size);
    
    //! \biref given a input file with N null Hypo toy trees it produce an out tree containing post fits and limit.
    //!
    //! note that this can be used also with just a signle tree, for real data limit production.
    //! @param ninety_quantiles: graph of 90% quantiles as produced from computeTSDistros().
    //! @param inputTreeFile: file containing the input toy tree list or data.
    void spitTheLimit(TGraphAsymmErrors *ninety_quantiles, TFile *inputTreeFile);

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