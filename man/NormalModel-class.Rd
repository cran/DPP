\name{NormalModel-class}
\alias{NormalModel}
\docType{class}
\title{Class \code{"NormalModel"}}
\description{
 Objects of the NormalModel class are initialized with prior parameters to be used by the MCMC algorithm in dppMCMC_C class.
 An object of the class NormalMode will be passed as an argument upon creation
 of the dppMCMC_C object that will run the MCMC code.
}


\examples{
#creating an object of the class NormalModel
normal.model<-new(NormalModel,
                  mean_prior_mean=0.5,
                  mean_prior_sd=0.1,
                  sd_prior_shape=3,
                  sd_prior_rate=20,
                  estimate_concentration_parameter=TRUE,
                  concentration_parameter_alpha=10,
                  proposal_disturbance_sd=0.1)
normal.model$getParameters()
}
\keyword{classes}

\section{Methods}{
  \preformatted{new(NormalModel,
    mean_prior_mean=0.5,
    mean_prior_sd=0.1,
    sd_prior_shape=3,
    sd_prior_rate=20,
    estimate_concentration_parameter=TRUE,
    concentration_parameter_alpha=10,
    proposal_disturbance_sd=0.1)
       instantiates a NormalModel object

getParameters()
  returns a list with the parameters and arguments supplied upon object initialization }

 }


