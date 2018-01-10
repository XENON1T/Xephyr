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
#include "plotHelpers.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"

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
    
    //! \brief constructor, CollectionName is the name of the tree collection to fit, 
    //! by convention set by ToyGenerator this is also the name of the file.
    ToyFitterExclusion(TString CollectionName);
    
    /**
     * \brief returns the Test Statistic q_tilde(mu) for a value of mu_test.
     * 
     * It is a reproduction of q_tilde test statistic from equation 16 of: https://arxiv.org/abs/1007.1727
     * this is suitable for upper limit in which the parameter of interest must be >0.
     * For it to work correctly is important to set up the likelihood beforehand.
     * @param mu: the signal strenght for conditional fit. 
     */
    double computeTS(double mu) ;

    /**
     * \brief given as input a file with many data toys, it writes a tree with fit outputs
     * 
     * Note that it expects a toy file generated with naming convention of ToyGenerator.
     * The output tree contains: conditional and unconditional nuisance parameter 
     * post fit values, and the test statistic. This is supposed to be used for the alternative 
     * hypotesis fits, to be able to generate aftewards the graph of 90th percentiles.
     * @parmam mu: the signal strenght of the conditional fit.
     * @parma soptAt:  optional, the number of toy in file you want to fit.
     * NOTE: by deafulte this will randomize the NP measures during fit for each toy, setRandomizeMeasure(false) if you don't want.
     */
    void fit(double mu, int stopAt=-999);

    //! \brief set the likelihood to fit
    void setTheLikelihood(pdfLikelihood *like) { likeHood = like; };


    /**
     * \brief produces a TGraph with the 90% quantiles of the alternate Hypothesis.
     * 
     * it is used to compute limits. The TH1F of the test statistic f(q_mu | H_mu) are
     * stored in a root file toghether with the TGraph of 90% quantiles.
     * @param fileName: takes as input the hadded post fit output from fit() of this class.
     * @param mu_list: is the list of true mu hypothesis present in the file tree.
     * @param mu_size: is the size of the previous list.
     */
    TGraphAsymmErrors computeTSDistros(TString fileName, double *mu_list, int mu_size);


    /**
     * \brief produces a TGraph with the 90% quantiles of the alternate Hypothesis.
     * 
     * it is used to compute limits. The TH1F of the test statistic f(q_mu | H_mu) are
     * stored in a root file toghether with the TGraph of 90% quantiles.
     * @param tree: takes as input the tree of post fit output from fit() of this class.
     * @param mu_list: is the list of true mu hypothesis present in the file tree.
     * @param mu_size: is the size of the previous list.
     */
    TGraphAsymmErrors computeTSDistros(TTree *tree, double *mu_list, int mu_size);
    
    /**
     * \biref given a input file with N null Hypo toy trees it produce an out tree containing post fits and limit.
     * 
     * note that this can be used also with just a signle tree, for real data limit production.
     * @param ninety_quantiles: graph of 90% quantiles as produced from computeTSDistros().
     * @param TreeName: prefix name of the tree list you fit. To se the input file file use setPathToFile() method.
     * @parma soptAt:  optional, the number of toy in file you want to fit.
     * NOTE: by deafulte this will randomize the NP measures during fit for each toy, setRandomizeMeasure(false) if you don't want. 
     */
    void spitTheLimit(TGraphAsymmErrors *ninety_quantiles,  int sopAt = -999);

    //! \brief set path to input file (only dir no file name)
    void setInputDir(TString path) { dirPath = path; };
    
    //! \brief set the path to dir in which the fits are stored otherwise is current
    void setOutputDir(TString path){ OutDir = path; };

    // \brief set or unset the parameter randomization. Default true.
    void setRandomizeMeasure(bool doOrNot) { randomizeMeasure = doOrNot ;};

    // \brief in case you want to change the tree name on the fly
    void setTreeName(TString newName) { treeName = newName; };

    // \brief setting random seed for measured param generation
    // TRandom likes a integer seed.
    void setSeed(ULong_t seed) { rambo.SetSeed(seed); }; 

    // \brief you must set this if you want to use calibration in fit
    void setCalibrationTreeName(TString newName) { calTreeName = newName; };

    // \brief set output file suffix, useful for conditional fit of null 
    void setOutputSuffix(TString name) { Suffix = name; } ;
    
  private:

    //! \brief fancy method that loops over a list of tree in a file and run the p2method() function on each.
    //!
    //! @params stopAt: number of tree one wants to loop on
    void for_each_tree( double (ToyFitterExclusion::*p2method)(double), TTree *outTree, double mu, int stopAt = -999);

    //! \brief compute the limit starting from initial_mu via a loop on computeTS
    //! and using graph_of_quantiles
    double limitLoop(double initial_mu);

    //! \brief save the postfit parameters into params.
    void saveParameters(double *params);

    //! \brief fills the initial paramteres value of true parameter an their name
    void fillTrueParams( TTree *inputTree);

    //! \brief assign to each parameter a random measured t-value
    void measureParameters();

    //! \brief function implementation to be used in limit minuit minimization
    double eval_testStatMinuit( double mu );

    void saveNames(string *names);

    pdfLikelihood *likeHood;
    TString       dirPath;
    TString       treeName;
    TString       OutDir;               //! output dir if not set is current dir
    TGraphAsymmErrors *graph_of_quantiles;
    TString       calTreeName;
    TString       Suffix;
    int           IndexHolder;
    int           CurrentTreeIndex;
    
    double        likelihood_uncond;    //! value of unconditional fit
    double        likelihood_cond;      //! value of conditional fit
    double        true_params[50] = {0.} ;    //! t_value of truth generated param
    double        measured_params[50] = {0.};  //! t_value of measured param
    double        uncond_params[50] = {0.};   //! unconditional fit
    double        cond_params[50] = {0.};     //! conditional fit
    bool          limit_converged;      //! if limit finder converged or not
    double        testStat_limit;       //! value of test statistic at limit
    vector<string> name_params;         //! names of parameters
    TRandom3      rambo;                //! random handler
    double        mu_fit;               //! mu to be tested, conditional fit
    double        mu_limit;             //! value of the limit in terms of mu
    double        limit;                //! limit scaled in cross section
    double        testStat;             //! value of test statistic
    int           numberOfParams;       //! number of likelihood parameter (including POI)
    bool          randomizeMeasure;     //! to random or not the np central value
};

#endif
