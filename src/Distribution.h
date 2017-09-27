#include <Rcpp.h>
//#include <RcppArmadillo.h>
using namespace Rcpp;

#ifndef DIST_H
#define DIST_H

RCPP_EXPOSED_CLASS(Distribution)
RCPP_EXPOSED_CLASS(Uniform)
RCPP_EXPOSED_CLASS(Normal)
RCPP_EXPOSED_CLASS(Gamma)
RCPP_EXPOSED_CLASS(Beta)
RCPP_EXPOSED_CLASS(Geometric)
RCPP_EXPOSED_CLASS(Poisson)
RCPP_EXPOSED_CLASS(Degenerate)

////////////////////////////////////
// pure virtual distribution type //
////////////////////////////////////

class Distribution {

public:

  virtual                     ~Distribution(){};
                              Distribution(){};

  virtual double              lnProb( std::vector<double> val ) = 0;
  virtual std::vector<double> sample( int n ) = 0;
  virtual List                getParameters() = 0;
  virtual int                 getDimensions() = 0;
  virtual bool                isFixed() = 0;

private:

};

//////////////////////////////
// the uniform distribution //
//////////////////////////////

class Uniform : public Distribution {

public:
                      ~Uniform(){};
                      Uniform(double min_, double max_);
                      Uniform(double min_, double max_, bool fixed_);
                      Uniform(const Uniform& other);

  double              lnProb( std::vector<double> val );
  std::vector<double> sample(int n);
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        min;
  const double        max;
  const bool          fixed;

};

/////////////////////////////
// the normal distribution //
/////////////////////////////

class Normal : public Distribution {

public:

                      ~Normal(){};
                      Normal(double m_, double v_);
                      Normal(double m_, double v_, bool fixed_);
                      Normal(const Normal& other);

  double              lnProb( std::vector<double> val );
  double              lnProb2(double val,double mean, double sd );
  double              rnorm(double mean, double sd );
  int                 sample_int(int max);
  int                 sample_int_prob(std::vector<double> probs);
  DoubleVector lnDNorm(DoubleVector val );

  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        mean;
  const double        variance;
  const bool          fixed;

};

////////////////////////////////
// the lognormal distribution //
////////////////////////////////

class Lognormal : public Distribution {

public:

  ~Lognormal(){};
  Lognormal(double m_, double s_);
  Lognormal(double m_, double s_, bool fixed_);
  Lognormal(const Lognormal& other);

  double              lnProb( std::vector<double> val );
  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        logmu;
  const double        logsigma;
  const bool          fixed;

};

////////////////////////////
// the gamma distribution //
////////////////////////////

class Gamma : public Distribution {

public:

                      ~Gamma(){};
                      Gamma(double s_, double r_);
                      Gamma(double s_, double r_, bool fixed_);
                      Gamma(const Gamma& other);

  double              lnProb( std::vector<double> val );
  double              lnProb2(double val,double shape, double rate );
  DoubleVector lnDGamma( DoubleVector val );
  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        shape;
  const double        rate;
  const bool          fixed;

};

///////////////////////////
// the beta distribution //
///////////////////////////

class Beta : public Distribution {

public:

                      ~Beta(){};
                      Beta(double a_, double b_);
                      Beta(double a_, double b_, bool fixed_);
                      Beta(const Beta& other);

  double              lnProb( std::vector<double> val );
  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        alpha;
  const double        beta;
  const bool          fixed;

};

////////////////////////////////
// the geometric distribution //
////////////////////////////////

class Geometric : public Distribution {

public:

                      ~Geometric(){};
                      Geometric(double p_, int offset_);
                      Geometric(double p_, int offset_, bool fixed_);
                      Geometric(const Geometric& other);

  double              lnProb( std::vector<double> val );
  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        p;
  const int           offset;
  const bool          fixed;

};

//////////////////////////////
// the poisson distribution //
//////////////////////////////

class Poisson : public Distribution {

public:

                      ~Poisson(){};
                      Poisson(double lambda_, int offset_);
                      Poisson(double lambda_, int offset_, bool fixed_);
                      Poisson(const Poisson& other);

  double              lnProb( std::vector<double> val );
  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return fixed; }

private:

  const double        lambda;
  const int           offset;
  const bool          fixed;

};

/////////////////////////////////
// the degenerate distribution //
/////////////////////////////////

class Degenerate : public Distribution {

public:

  ~Degenerate(){};
  Degenerate(double points_);
  Degenerate(const Degenerate& other);

  double              lnProb( std::vector<double> val );
  std::vector<double> sample( int n );
  List                getParameters();
  int                 getDimensions() { return 1; }
  bool                isFixed() { return true; }

private:

  const double        point;

};


#endif





