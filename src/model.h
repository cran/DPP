#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include "Distribution.h"
using namespace Rcpp;

#ifndef MODEL_H
#define MODEL_H

//RCPP_EXPOSED_CLASS(Model)
RCPP_EXPOSED_CLASS(NormalModel)
RCPP_EXPOSED_CLASS(GammaModel)


  ////////////////////////////////////
  // pure virtual distribution type //
  ////////////////////////////////////

  class Model {

  public:

    virtual                     ~Model(){};
                                 Model(){};

  //  virtual double              lnProb( std::vector<double> val ) = 0;
  //  virtual std::vector<double> sample( int n ) = 0;
    virtual List                getParameters() = 0;
   // virtual int                 getDimensions() = 0;
    virtual std::vector<double> likelihood_fn(DoubleVector data,IntegerVector allocation,Rcpp::List params, int power)=0;
    virtual std::vector<double> single_likelihood_fn(double data,IntegerVector allocation,Rcpp::List params, int power)=0;
    virtual Rcpp::List base_distn_sim(int num_categories)=0;
    virtual DoubleVector base_distn(Rcpp::List params)=0;
    virtual Rcpp::List proposal_distn(Rcpp::List params)=0;
    bool                 getEstimateConcentrationParameter() { return estimate_concentration_parameter; }
    double                 getConcentrationParameterAlpha() { return concentration_parameter_alpha; }
    double                 getProposalDisturbanceSd() { return proposal_disturbance_sd; }

   // virtual bool                isFixed() = 0;

  protected:
     bool estimate_concentration_parameter;
     double concentration_parameter_alpha;
     double proposal_disturbance_sd;

  };

//////////////////////////////
// the Normal Model         //
//////////////////////////////

class NormalModel : public Model {

public:
  ~NormalModel(){};
  NormalModel(double mean_prior_mean_,
              double mean_prior_mean_sd_,
              double sd_prior_shape_,
              double sd_prior_rate_,
              bool estimate_concentration_parameter_,
              double concentration_parameter_alpha_,
              double proposal_disturbance_sd_

                );

  NormalModel(const NormalModel& other);


  List                getParameters();


  std::vector<double> likelihood_fn(DoubleVector data,IntegerVector allocation,Rcpp::List params, int power);
  std::vector<double> single_likelihood_fn(double data,IntegerVector allocation,Rcpp::List params, int power);
  Rcpp::List base_distn_sim(int num_categories);
  DoubleVector base_distn(Rcpp::List params);
  Rcpp::List proposal_distn(Rcpp::List params);

private:
  const double        mean_prior_mean;
  const double        mean_prior_sd;
  const double        sd_prior_shape;
  const double        sd_prior_rate;


};

//////////////////////////////
// the Gamma Model         //
//////////////////////////////

class GammaModel : public Model {

public:
  ~GammaModel(){};
  GammaModel(double shape_prior_mean_,
              double shape_prior_sd_,
              double rate_prior_mean_,
              double rate_prior_sd_,
              bool estimate_concentration_parameter_,
              double concentration_parameter_alpha_,
              double proposal_disturbance_sd_

  );

  GammaModel(const GammaModel& other);


  List                getParameters();


  std::vector<double> likelihood_fn(DoubleVector data,IntegerVector allocation,Rcpp::List params, int power);
  std::vector<double> single_likelihood_fn(double data,IntegerVector allocation,Rcpp::List params, int power);
  Rcpp::List base_distn_sim(int num_categories);
  DoubleVector base_distn(Rcpp::List params);
  Rcpp::List proposal_distn(Rcpp::List params);

private:
  const double        shape_prior_mean;
  const double        shape_prior_sd;
  const double        rate_prior_mean;
  const double        rate_prior_sd;


};

#endif





