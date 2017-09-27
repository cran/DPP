// #include <Rcpp.h>
#include "model.h"
#include <limits>
//using namespace Rcpp;

#include <algorithm>
#include <functional>
#include "helper.h"


//////////////////////////////
// the normal mode          //
//////////////////////////////



NormalModel::NormalModel(double mean_prior_mean_,
            double mean_prior_sd_,
            double sd_prior_shape_,
            double sd_prior_rate_,
            bool estimate_concentration_parameter_,
            double concentration_parameter_alpha_,
            double proposal_disturbance_sd_
            ):mean_prior_mean(mean_prior_mean_),
                                  mean_prior_sd(mean_prior_sd_),
                                  sd_prior_shape(sd_prior_shape_),
                                  sd_prior_rate(sd_prior_rate_)//,
                                  {
  estimate_concentration_parameter=estimate_concentration_parameter_;
  concentration_parameter_alpha=concentration_parameter_alpha_;
  proposal_disturbance_sd=proposal_disturbance_sd_;
}

std::vector<double> NormalModel::likelihood_fn(DoubleVector data,IntegerVector allocation,Rcpp::List params, int power)  {

    Normal normal_dist=Normal(mean_prior_mean,mean_prior_sd*mean_prior_sd);//takes variance

    std::vector<double> params1=params(0);
    std::vector<double> params2=params(1);

    int dataSize=data.size();
    double mean;
    double sd;
    std::vector<double> ll(dataSize);

    for (int i=0;i<dataSize;i++){
      mean=params1[allocation[i]-1]; //-1 because vector array is 0 based
      sd=params2[allocation[i]-1];
      ll[i]=normal_dist.lnProb2(data[i],mean,sd)*power;
    }
  return ll;
}
std::vector<double> NormalModel::single_likelihood_fn(double data,IntegerVector allocation,Rcpp::List params, int power){

  Normal normal_dist=Normal(mean_prior_mean,mean_prior_sd*mean_prior_sd);//takes variance

  std::vector<double> params1=params(0);
  std::vector<double> params2=params(1);
 // int vectorSize=params1.size();
  int dataSize=allocation.size();
  double mean;
  double sd;
  std::vector<double> ll(dataSize);

  for (int i=0;i<dataSize;i++){
    mean=params1[allocation[i]-1]; //-1 because vector array is 0 based
    sd=params2[allocation[i]-1];
    ll[i]=normal_dist.lnProb2(data,mean,sd)*power;
  }
  return ll;
}

Rcpp::List NormalModel::base_distn_sim(int num_categories)
{

  Normal normal_dist=Normal(mean_prior_mean,mean_prior_sd*mean_prior_sd);
  DoubleVector means(num_categories);
  means=normal_dist.sample(num_categories);
  Gamma gamma_dist=Gamma(sd_prior_shape,sd_prior_rate);
  DoubleVector gamma_sds(num_categories);
  gamma_sds=gamma_dist.sample(num_categories);
  return Rcpp::List::create(Named("means")=means,Named("sds")=gamma_sds);

}

DoubleVector NormalModel::base_distn(Rcpp::List params)
{
  //initializing distributions

  Normal normal_dist=Normal(mean_prior_mean,mean_prior_sd*mean_prior_sd);//takes variance
  Gamma gamma_dist=Gamma(sd_prior_shape,sd_prior_rate);

  DoubleVector params1=params(0);
  DoubleVector params2=params(1);

  DoubleVector lp(params1.size());
   lp=normal_dist.lnDNorm(params1);
    lp=lp+gamma_dist.lnDGamma(params2);
 return lp;
}


Rcpp::List NormalModel::proposal_distn(Rcpp::List params){
  std::vector<double> params1=params(0);
  std::vector<double> params2=params(1);
  int params_length=params.length();

  std::vector<double> new_params1=params1;
  std::vector<double> new_params2=params2;

  Normal normal_dist=Normal(0,1);

 // choose a parameter at random
  int j;
  int k = normal_dist.sample_int(params_length);

  if(k==1) {
     j=normal_dist.sample_int(params1.size());
     new_params1[j-1]=params1[j-1]+normal_dist.rnorm(0,proposal_disturbance_sd);
  } else { //if  k==2

    j=normal_dist.sample_int(params2.size());
    new_params2[j-1]=params2[j-1]+normal_dist.rnorm(0,proposal_disturbance_sd);
    // if the new variance is negative,
    // reflect it back into the positive
    new_params2[j-1]=fabs(new_params2[j-1]);
  }

  return Rcpp::List::create(Named("means")=new_params1,Named("sds")=new_params2);
}


List NormalModel::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("mean_prior_mean")=mean_prior_mean,
                                   Rcpp::Named("mean_prior_sd")=mean_prior_sd,
                                   Rcpp::Named("sd_prior_shape")=sd_prior_shape,
                                   Rcpp::Named("sd_prior_rate")=sd_prior_rate,
                                   Rcpp::Named("estimate_concentration_parameter")=estimate_concentration_parameter,
                                   Rcpp::Named("concentration_parameter_alpha")=concentration_parameter_alpha,
                                   Rcpp::Named("proposal_disturbance_sd")=proposal_disturbance_sd);

  return params;

}

//////////////////////////////
// the gamma mode          //
//////////////////////////////



GammaModel::GammaModel(  double shape_prior_mean_,
                         double shape_prior_sd_,
                         double rate_prior_mean_,
                         double rate_prior_sd_,
                         bool estimate_concentration_parameter_,
                         double concentration_parameter_alpha_,
                         double proposal_disturbance_sd_
):shape_prior_mean(shape_prior_mean_),
shape_prior_sd(shape_prior_sd_),
rate_prior_mean(rate_prior_mean_),
rate_prior_sd(rate_prior_sd_)//,
{
  estimate_concentration_parameter=estimate_concentration_parameter_;
  concentration_parameter_alpha=concentration_parameter_alpha_;
  proposal_disturbance_sd=proposal_disturbance_sd_;
}

std::vector<double> GammaModel::likelihood_fn(DoubleVector data,IntegerVector allocation,Rcpp::List params, int power)  {

  Gamma gamma_dist=Gamma(shape_prior_mean,shape_prior_sd);//this values won't be used, but gamma dist
                                                           //needs to be initialized

  std::vector<double> params1=params(0); //shape
  std::vector<double> params2=params(1); //rate

  int dataSize=data.size();
  double shape;
  double rate;
  std::vector<double> ll(dataSize);

  for (int i=0;i<dataSize;i++){
    shape=params1[allocation[i]-1]; //-1 because vector array is 0 based
    rate=params2[allocation[i]-1];
    ll[i]=gamma_dist.lnProb2(data[i],shape,rate)*power;
  }
  return ll;
}
std::vector<double> GammaModel::single_likelihood_fn(double data,IntegerVector allocation,Rcpp::List params, int power){

  Gamma gamma_dist=Gamma(shape_prior_mean,shape_prior_sd);//this values won't be used, but gamma dist
  //needs to be initialized

  std::vector<double> params1=params(0);
  std::vector<double> params2=params(1);
  // int vectorSize=params1.size();
  int dataSize=allocation.size();
  double shape;
  double rate;
  std::vector<double> ll(dataSize);

  for (int i=0;i<dataSize;i++){
    shape=params1[allocation[i]-1]; //-1 because vector array is 0 based
    rate=params2[allocation[i]-1];
    ll[i]=gamma_dist.lnProb2(data,shape,rate)*power;
  }
  return ll;
}

Rcpp::List GammaModel::base_distn_sim(int num_categories)
{

  Normal shape_normal_dist=Normal(shape_prior_mean,shape_prior_sd*shape_prior_sd); //takes variance



  DoubleVector shapes(num_categories);
  shapes=shape_normal_dist.sample(num_categories);
  for(int i=0;i<shapes.size();i++){shapes[i]=fabs(shapes[i]);} //making sure this don't go negative

  Normal rate_normal_dist=Normal(rate_prior_mean,rate_prior_sd *rate_prior_sd ); //takes variance
  DoubleVector rates(num_categories);
  rates=rate_normal_dist.sample(num_categories);
  for(int i=0;i<rates.size();i++){rates[i]=fabs(rates[i]);} //making sure this don't go negative
  return Rcpp::List::create(Named("shapes")=shapes,Named("rates")=rates);

}

DoubleVector GammaModel::base_distn(Rcpp::List params)
{
  //initializing distributions

  Normal shapes_normal_dist=Normal(shape_prior_mean,shape_prior_sd*shape_prior_sd);//takes variance
  Normal rates_normal_dist=Normal(rate_prior_mean,rate_prior_sd*rate_prior_sd);//takes variance
  //Gamma gamma_dist=Gamma(sd_prior_shape,sd_prior_rate);

  DoubleVector params1=params(0);
  DoubleVector params2=params(1);

  DoubleVector lp(params1.size());
  lp=shapes_normal_dist.lnDNorm(params1);
  lp=lp+rates_normal_dist.lnDNorm(params2);
  return lp;
}


Rcpp::List GammaModel::proposal_distn(Rcpp::List params){
  std::vector<double> params1=params(0);
  std::vector<double> params2=params(1);
  int params_length=params.length();

  std::vector<double> new_params1=params1;
  std::vector<double> new_params2=params2;

  Normal normal_dist=Normal(0,1);

  // choose a parameter at random
  int j;
  int k = normal_dist.sample_int(params_length);

  if(k==1) {
    j=normal_dist.sample_int(params1.size());
    new_params1[j-1]=fabs(params1[j-1]+normal_dist.rnorm(0,proposal_disturbance_sd)); //fabs 'cause this should > 0
  } else { //if  k==2

    j=normal_dist.sample_int(params2.size());
    new_params2[j-1]=fabs(params2[j-1]+normal_dist.rnorm(0,proposal_disturbance_sd)); //fabs 'cause this should > 0

  }

  return Rcpp::List::create(Named("shapes")=new_params1,Named("rates")=new_params2);
}


List GammaModel::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("shape_prior_mean")=shape_prior_mean,
                                   Rcpp::Named("shape_prior_sd")=shape_prior_sd,
                                   Rcpp::Named("rate_prior_mean")=rate_prior_mean,
                                   Rcpp::Named("rate_prior_sd")=rate_prior_sd,
                                   Rcpp::Named("estimate_concentration_parameter")=estimate_concentration_parameter,
                                   Rcpp::Named("concentration_parameter_alpha")=concentration_parameter_alpha,
                                   Rcpp::Named("proposal_disturbance_sd")=proposal_disturbance_sd);
  return params;

}
/////////////////
// Rcpp Module //
/////////////////

RCPP_MODULE(Models) {

  using namespace Rcpp;

 class_<Model>( "Model" )


    .method("getParameters", &Model::getParameters)
    .method("likelihood_fn",&Model::likelihood_fn)
    .method("single_likelihood_fn",&Model::single_likelihood_fn)
    .method("base_distn_sim",&Model::base_distn_sim)
    .method("base_distn",&Model::base_distn)
    .method("proposal_distn",&Model::proposal_distn)
    .method("getConcentrationParameterAlpha",&NormalModel::getConcentrationParameterAlpha)
    .method("getEstimateConcentrationParameter",&NormalModel::getEstimateConcentrationParameter)

    ;

  class_<NormalModel>( "NormalModel" )
    .derives<Model>( "Model" )
    .constructor<double, double, double,double,bool,double,double/*,double*/>()

      ;

  class_<GammaModel>( "GammaModel" )
    .derives<Model>( "Model" )
    .constructor<double, double, double,double,bool,double,double/*,double*/>()

  ;

}



