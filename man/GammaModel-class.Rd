\name{GammaModel-class}
\alias{GammaModel}
\docType{class}
\title{Class \code{"GammaModel"}}
\description{
 Objects of the GammaModel class are initialized with prior parameters to be used by the MCMC algorithm in dppMCMC_C class.
 An object of the class GammaModel will be passed as an argument upon creation
 of the dppMCMC_C object that will run the MCMC code.
}


\examples{
#creating an object of the class NormalModel
gamma.model<-new(GammaModel,
                    shape_prior_mean=4,
                    shape_prior_sd=1,
                    rate_prior_mean=1.5,
                    rate_prior_sd=0.54,
                    estimate_concentration_parameter=TRUE,
                    concentration_parameter_alpha=10,
                    proposal_disturbance_sd=0.1)
gamma.model$getParameters()
}
\keyword{classes}

\section{Methods}{
  \preformatted{new(GammaModel,
                    shape_prior_mean=4,
                    shape_prior_sd=1,
                    rate_prior_mean=1.5,
                    rate_prior_sd=0.54,
                    estimate_concentration_parameter=TRUE,
                    concentration_parameter_alpha=10,
                    proposal_disturbance_sd=0.1)
       instantiates a GammaModel object

getParameters()
  returns a list with the parameters and arguments supplied upon object initialization }

 }


