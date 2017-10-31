\name{DPP-package}
\alias{DPP-package}
\alias{DPP}
\alias{DPPmcmc}
\alias{Model}
\alias{expectedNumberOfClusters}
\alias{Rcpp_DPPmcmc}
\alias{Rcpp_Model}
\docType{package}
\title{
\packageTitle{DPP}
}
\description{
\packageDescription{DPP}
}
\details{
 DPP implements a Bayesian method to get a posterior
 probability for k normal distributions form a vector of n numeric values.
We implemented an MCMC method as described in Escobar (1994). Using a Dirichlet process prior we describe stochastically our prior assumptions about the dimensionality of the data without specifying a fix number of clusters k, allowing us to infer the number of normal distributions or categories from a potentially infinite number of categories.
DPP is implemented in C++ and made available to used within the R statistical environment using Rcpp (Eddelbuettel and Francois, 2011).

}

\author{
 Luis M. Avila, Michael R. May, Jeff Ross-Ibarra

Maintainer: Luis M. Avila <lmavila@gmail.com>
}
\references{
Ferguson, Thomas S. A Bayesian analysis of some nonparametric problems. The annals of statistics (1973): 209-230.

Antoniak, Charles E. Mixtures of Dirichlet processes with applications to Bayesian nonparametric problems. The Annals of Statistics (1974): 1152-1174.

Escobar, Michael D. Estimating Normal Means With a Dirichlet Process Prior. Journal of the American Statistical Association, 89(425), 1994.

Neal, Radford M. Markov chain sampling methods for Dirichlet process mixture models. Journal of Computational and Graphical Statistics 9.2 (2000): 249-265.

Eddelbuettel, Dirk and Romain Francois. Rcpp: Seamless R and C++ Integration. Journal Of Statistical Software, 40(8):1-18, 2011.
}

\keyword{ package }

\examples{
normal.model<-new(NormalModel,
                  mean_prior_mean=0.5,
                  mean_prior_sd=0.1,
                  sd_prior_shape=3,
                  sd_prior_rate=20,
                  estimate_concentration_parameter=TRUE,
                  concentration_parameter_alpha=10,
                  proposal_disturbance_sd=0.1)

 #simulating three normal distributions
 y <- c(rnorm(100,mean=0.2,sd=0.05), rnorm(100,0.7,0.05), rnorm(100,1.3,0.1))
 hist(y,breaks=30)

 #setwd("~/yourwd") #mcmc log files will be saved here
 my_dpp_analysis <- dppMCMC_C(data=y,
                              output = "output_prefix_",
                              model=normal.model,
                              num_auxiliary_tables=4,
                              expected_k=1.5,
                              power=1)
 #running the mcmc  , generations will be ignored because auto_stop=true
 \dontrun{
 my_dpp_analysis$run(generations=1000,auto_stop=TRUE,max_gen = 10000,min_ess = 500)

 #we get rid of the first 25\% of the output (burn-in)
 hist(my_dpp_analysis$getNumCategoryTrace(0.25))
 my_dpp_analysis$getNumCategoryProbabilities(0.25)
 }
}


