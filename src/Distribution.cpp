//#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include "Distribution.h"
#include <limits>
using namespace Rcpp;

//////////////////////////////
// the uniform distribution //
//////////////////////////////

Uniform::Uniform(double min_, double max_) : min(min_), max(max_), fixed(false) {

}

Uniform::Uniform(double min_, double max_, bool fixed_) : min(min_), max(max_), fixed(fixed_) {

}

Uniform::Uniform(const Uniform& other) : min(other.min), max(other.max), fixed(other.fixed) {

}

double Uniform::lnProb( std::vector<double> val ) {

  double lnp = 0.0;

  if (fixed) {
    lnp = ::Rf_dunif(val.at(0),min,max,1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dunif(val.at(i),min,max,1);
    }
  }

  return lnp;

}

std::vector<double> Uniform::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_runif(min,max);
    std::fill(res.begin(),res.end(),f);
  } else {
    for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_runif(min,max);
    }
  }
  return res;

}

List Uniform::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("min")=min,
                                   Rcpp::Named("max")=max);

  return params;

}

/////////////////////////////
// the normal distribution //
/////////////////////////////

Normal::Normal(double m_, double v_) : mean(m_), variance(v_), fixed(false) {

}

Normal::Normal(double m_, double v_, bool fixed_) : mean(m_), variance(v_), fixed(fixed_) {

}

Normal::Normal(const Normal& other) : mean(other.mean), variance(other.variance), fixed(other.fixed) {

}
double Normal::rnorm(double mean, double sd)
{
  return R::rnorm(mean,sd);
}
double Normal::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  if (fixed) {
    lnp = ::Rf_dnorm4(val.at(0),mean,sqrt(variance),1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dnorm4(val.at(i),mean,sqrt(variance),1);
    }
  }
  return lnp;

}
double Normal::lnProb2(double val,double mean, double sd ){
  return Rf_dnorm4(val,mean,sd,1);
}
DoubleVector Normal::lnDNorm(DoubleVector val ){
  int size = val.size();
  DoubleVector lnp_vector(size);



    for(int i = 0; i < size; ++i) {
      lnp_vector[i] = ::Rf_dnorm4(val.at(i),mean,sqrt(variance),1);

    }
  return lnp_vector;

}

int Normal::sample_int(int max)
{ //mimics R sample.int sampling one value between 1 and max

  double returnVal;
  RNGScope rngScope;
  returnVal=R::runif(1,max+0.99999);
  return std::floor(returnVal);

}

int Normal::sample_int_prob(std::vector<double> probs)
{
  std::vector<double> limits(probs.size());
  int chosen_int=probs.size();
  RNGScope rngScope;
  limits[0]=probs[0];


  for(int i=1;i<probs.size();i++){
       limits[i]=limits[i-1]+probs[i];
  }

  double sampledVal=  ::Rf_runif(0,limits[limits.size()-1]);
  //double sampledVal=  ::Rf_runif(0,1);


  for(int i=0;i<limits.size();i++){
    if(sampledVal<limits[i]){
      chosen_int=i+1;
      break;
    }
    //limits[i]=limits[i-1]+probs[i];
  }

  return chosen_int;

}


std::vector<double> Normal::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_rnorm(mean,sqrt(variance));
    std::fill(res.begin(),res.end(),f);
  } else {
   for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_rnorm(mean,sqrt(variance));
    }
  }
  return res;

}

List Normal::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("mean")=mean,
                                   Rcpp::Named("variance")=variance);

  return params;

}

////////////////////////////////
// the lognormal distribution //
////////////////////////////////

Lognormal::Lognormal(double m_, double s_) : logmu(m_), logsigma(s_), fixed(false) {

}

Lognormal::Lognormal(double m_, double s_, bool fixed_) : logmu(m_), logsigma(s_), fixed(fixed_) {

}

Lognormal::Lognormal(const Lognormal& other) : logmu(other.logmu), logsigma(other.logsigma), fixed(other.fixed) {

}

double Lognormal::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  if (fixed) {
    lnp = ::Rf_dlnorm(val.at(0),logmu,logsigma,1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dlnorm(val.at(i),logmu,logsigma,1);
    }
  }
  return lnp;

}




std::vector<double> Lognormal::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_rlnorm(logmu,logsigma);
    std::fill(res.begin(),res.end(),f);
  } else {
    for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_rlnorm(logmu,logsigma);
    }
  }
  return res;

}

List Lognormal::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("logmu")=logmu,
                                   Rcpp::Named("logsigma")=logsigma);

  return params;

}

////////////////////////////
// the gamma distribution //
////////////////////////////

Gamma::Gamma(double s_, double r_) : shape(s_), rate(r_), fixed(false) {

}

Gamma::Gamma(double s_, double r_, bool fixed_) : shape(s_), rate(r_), fixed(fixed_) {

}

Gamma::Gamma(const Gamma& other) : shape(other.shape), rate(other.rate), fixed(other.fixed) {

}

double Gamma::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  if (fixed) {
    lnp = ::Rf_dgamma(val.at(0),shape,1.0/rate,1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dgamma(val.at(i),shape,1.0/rate,1);
    }
  }
  return lnp;

}
double Gamma::lnProb2(double val,double shape, double rate ){
  //return Rf_dnorm4(val,mean,sd,1);
  return Rf_dgamma(val,shape,1.0/rate,1);
  //Rf_dgamma(x, shp, scl, lg)
}
std::vector<double> Gamma::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_rgamma(shape,1.0/rate);
    std::fill(res.begin(),res.end(),f);
  } else {
    for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_rgamma(shape,1.0/rate);
    }
  }
  return res;

}

DoubleVector Gamma::lnDGamma( DoubleVector val )
{
  RNGScope scope;
  int size = val.size();
  DoubleVector lnp_vector(size);

    for(int i = 0; i < size; ++i) {
      lnp_vector[i]= ::Rf_dgamma(val.at(i),shape,1.0/rate,1);

    }
  return lnp_vector;

}

List Gamma::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("shape")=shape,
                                   Rcpp::Named("rate")=rate);

  return params;

}

///////////////////////////
// the beta distribution //
///////////////////////////

Beta::Beta(double a_, double b_) : alpha(a_), beta(b_), fixed(false) {

}

Beta::Beta(double a_, double b_, bool fixed_) : alpha(a_), beta(b_), fixed(fixed_) {

}

Beta::Beta(const Beta& other) : alpha(other.alpha), beta(other.beta), fixed(other.fixed) {

}

double Beta::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  if (fixed) {
    lnp = ::Rf_dbeta(val.at(0),alpha,beta,1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dbeta(val.at(i),alpha,beta,1);
    }
  }
  return lnp;

}

std::vector<double> Beta::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_rbeta(alpha,beta);
    std::fill(res.begin(),res.end(),f);
  } else {
    for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_rbeta(alpha,beta);
    }
  }
  return res;

}

List Beta::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                                   Rcpp::Named("beta")=beta);

  return params;

}

////////////////////////////////
// the geometric distribution //
////////////////////////////////

Geometric::Geometric(double p_, int offset_) : p(p_), offset(offset_), fixed(false) {

}

Geometric::Geometric(double p_, int offset_, bool fixed_) : p(p_), offset(offset_), fixed(fixed_) {

}

Geometric::Geometric(const Geometric& other) : p(other.p), offset(other.offset), fixed(other.fixed) {

}

double Geometric::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  if (fixed) {
    lnp = ::Rf_dgeom(val.at(0) - offset,p,1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dgeom(val.at(i) - offset,p,1);
    }
  }
  return lnp;

}

std::vector<double> Geometric::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_rgeom(p) + offset;
    std::fill(res.begin(),res.end(),f);
  } else {
    for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_rgeom(p) + offset;
    }
  }
  return res;

}

List Geometric::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("p")=p);

  return params;

}

//////////////////////////////
// the poisson distribution //
//////////////////////////////

Poisson::Poisson(double lambda_, int offset_) : lambda(lambda_), offset(offset_), fixed(false) {

}

Poisson::Poisson(double lambda_, int offset_, bool fixed_) : lambda(lambda_), offset(offset_), fixed(fixed_) {

}

Poisson::Poisson(const Poisson& other) : lambda(other.lambda), offset(other.offset), fixed(other.fixed) {

}

double Poisson::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  if (fixed) {
    lnp = ::Rf_dpois(val.at(0) - offset,lambda,1);
  } else {
    int size = val.size();
    for(int i = 0; i < size; ++i) {
      lnp += ::Rf_dpois(val.at(i) - offset,lambda,1);
    }
  }
  return lnp;

}

std::vector<double> Poisson::sample( int n ){

  RNGScope scope;
  std::vector<double> res(n,0.0);
  if (fixed) {
    double f = ::Rf_rpois(lambda) + offset;
    std::fill(res.begin(),res.end(),f);
  } else {
    for(int i = 0; i < n; ++i) {
      res.at(i) = ::Rf_rpois(lambda) + offset;
    }
  }
  return res;

}

List Poisson::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("lambda")=lambda);

  return params;

}

/////////////////////////////////
// the degenerate distribution //
/////////////////////////////////

Degenerate::Degenerate(double point_) : point(point_) {

}

Degenerate::Degenerate(const Degenerate& other) : point(other.point) {

}

double Degenerate::lnProb( std::vector<double> val ){

  double lnp = 0.0;
  for(int i = 0; i < val.size(); ++i) {
    if(val.at(i) != point) {
      lnp = std::numeric_limits<double>::infinity();
      break;
    }
  }
  return lnp;

}

std::vector<double> Degenerate::sample( int n ){

  std::vector<double> res(n,point);
  return res;

}

List Degenerate::getParameters(){

  List params = Rcpp::List::create(Rcpp::Named("point")=point);

  return params;

}

/////////////////
// Rcpp Module //
/////////////////

RCPP_MODULE(Distributions) {

  using namespace Rcpp;
/*
  class_<Distribution>( "Distribution" )

    .method("lnProb", &Distribution::lnProb)
    .method("sample", &Distribution::sample)
    .method("getParameters", &Distribution::getParameters)
    .method("getDimensions", &Distribution::getDimensions)

  ;

  class_<Uniform>( "Uniform" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, double>()
    .constructor<double, double, bool>()
  ;

  class_<Normal>( "Normal" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, double>()
    .constructor<double, double, bool>()
    .method("lnDNorm",&Normal::lnDNorm)
    .method("lnProb2",&Normal::lnProb2)
    .method("sample_int",&Normal::sample_int)
    .method("sample_int_prob",&Normal::sample_int_prob)


  ;

  class_<Lognormal>( "Lognormal" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, double>()
    .constructor<double, double, bool>()
  ;

  class_<Gamma>( "Gamma" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, double>()
    .constructor<double, double, bool>()
    .method("lnDGamma",&Gamma::lnDGamma)
  ;

  class_<Beta>( "Beta" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, double>()
    .constructor<double, double, bool>()
  ;

  class_<Geometric>( "Geometric" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, int>()
    .constructor<double, int, bool>()
  ;

  class_<Poisson>( "Poisson" )
    .derives<Distribution>( "Distribution" )
    .constructor<double, int>()
    .constructor<double, int, bool>()
  ;

  class_<Degenerate>( "Degenerate" )
    .derives<Distribution>( "Distribution" )
    .constructor<double>()
  ;
*/
}



