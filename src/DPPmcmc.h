#include "Distribution.h"
#include "model.h"
using namespace Rcpp;

#ifndef DPPmcmc_H
#define DPPmcmc_H

RCPP_EXPOSED_CLASS(DPPmcmc)

class DPPmcmc {

public:
                       /*DPPmcmc(DoubleVector data_, NormalModel& normal_model,int num_auxiliary_tables_ ,double expected_k,int power_,
                               Function effectiveSizeFunction_,
                               Function pminFunction_,
                               Function textBarFunction_);*/
                       DPPmcmc(DoubleVector data_, Model& model,int num_auxiliary_tables_ ,double expected_k,int power_,Function effectiveSizeFunction_, Function pminFunction_);

                      ~DPPmcmc(){};

  // simulation functions
  IntegerVector       simulateChineseRestaurant(int num_elements_, double alpha_);
  void                run(int generations,bool auto_stop,int max_gen,double min_ess,bool random,int sample_freq);
  double              expectedNumberOfClusters(int num_elements_, double alpha);
  double              concentrationParameterFromK(int n, double k);
  //helper methods
  IntegerVector       integerSequence(int min,int max);
  std::vector<double> divideIntegerVectorByDouble(IntegerVector vector1,double denominator);
  std::vector<double> rep(double a,int num_reps);
  IntegerVector intRep(int a,int num_reps);

  // accessor functions

  int                 getNumElements() { return num_elements; }
  int                 getPower() { return power; }
  DoubleVector        getData() { return data; }
  int                 getNumAuxiliaryTables() { return num_auxiliary_tables; }
  int                 getNumCategories() { return num_categories; }
  int                 getNumParams() { return num_params; }
  std::vector<int>    getNumElementsPerTable() { return num_elements_per_table; }
  std::vector<int>    getNumCategoryTrace() { return return_num_cats_trace; }
  double              getConcentrationParameter() { return concentration_parameter; }
  double              getConcentrationParameterAlpha() { return concentration_parameter_alpha; }
  double              getConcentrationParameterBeta() { return concentration_parameter_beta; }
  bool                getEstimateConcentrationParameter() { return estimate_concentration_parameter; }

  IntegerVector       getAllocationVector() { return allocation_vector; }
  Rcpp::List          getParamVector() { return param_vector; }
  IntegerVector       getNumElementsInEachCategory(){return num_elements_in_each_category;}

  void                 setOutputPrefix(std::string);
  void                 setVerbose(bool);
  void                 setSampleNumClusters(bool);
  std::string                 getOutputPrefix(){ return outputPrefix; }
  // MCMC functions

  std::vector<double>       dummyFunction(std::vector<double> dummyInput);
  void writeOutputFiles(void);
  void makeOutputFiles(void);
  void postInitialization(void);

private:
  DoubleVector        data;
  int                 num_elements;
  int                 num_params;
  int                 power;
  int                 num_categories;
  std::vector<int>    num_elements_per_table;
  std::vector<int>    return_num_cats_trace;
  IntegerVector num_elements_in_each_category;
  int                 num_auxiliary_tables;

  double              concentration_parameter;
  double              concentration_parameter_alpha;
  double              concentration_parameter_beta;
  bool                estimate_concentration_parameter;
  bool                verbose;
  bool                sample_num_clusters;
  double              likelihood;
  double              prior;
  int                 generation;
  double min_ESS;

  Model&       model;
  Function effectiveSizeFunction;
  Function pminFunction;
  std::string         outputPrefix;
  //Function textBarFunction;

  IntegerVector       allocation_vector;
  Rcpp::List param_vector;
  std::vector<double> proposed_parameters;
  std::vector<double> placeholder;

  double              ln_prior_ratio; // for the parameter proposal
  double              ln_hastings_ratio; // for the parameter proposal
  //private methods

  void                allocationProposal(int element);
  void                concentrationParameterProposal(void);
  //private helper methodes
  std::vector<double> makeIntegerVectorStandardDoubleVector(IntegerVector vector1);
  std::vector<double> makeDoubleVectorStandardDoubleVector(DoubleVector vector1);
  std::vector<int> makeIntegerVectorStandardIntVector(IntegerVector vector1);
  std::vector<int> evaluateVectorGreaterThanInt(std::vector<int> vector1,int val);
  void                updateConcentrationParameter();
  void                write_text_to_log_file( const std::string &text );

};

#endif
