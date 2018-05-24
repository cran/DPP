#'Class dppMCMC_C
#'Class \code{dppMCMC_C} A Reference Class that provides DPP functionality
#'
#'@name dppMCMC_C
#'@rdname dppMCMC_C
#'@exportClass dppMCMC_C
#'
#'@import coda methods
#'@title A Reference Class that provides DPP functionality
#'
#'
#' @description This class implements the main functionality of this package.
#' The consturctor receives a numeric vector (Y) and priors. A model should be provided specifying  the distributions
#' to be used for inference (e.g. NormalModel for Normal distributions or GammaModel for Gamma distributions).
#' Then an MCMC algorithm will be used to infer a number of distributions (k) that fit the data.
#' The prior for the number of distributions is specified by the concentration_parameter_alpha and expected_k.
#' Once the data and priors are specified the method run is used to start the inference.
#'
#' @field  dpp_mcmc_object a DPPmcmc object
#' @examples
#' normal.model<-new(NormalModel,
#'                  mean_prior_mean=0.5,
#'                  mean_prior_sd=0.1,
#'                  sd_prior_shape=3,
#'                  sd_prior_rate=20,
#'                  estimate_concentration_parameter=TRUE,
#'                  concentration_parameter_alpha=10,
#'                   proposal_disturbance_sd=0.1)
#'
#' #simulating three normal distributions
#' y <- c(rnorm(100,mean=0.2,sd=0.05), rnorm(100,0.7,0.05), rnorm(100,1.3,0.1))
#' hist(y,breaks=30)
#'
#' #setwd("~/yourwd") #mcmc log files will be saved here
#' \dontrun{
#' my_dpp_analysis <- dppMCMC_C(data=y,
#'                              output = "output_prefix_",
#'                              model=normal.model,
#'                              num_auxiliary_tables=4,
#'                              expected_k=1.5,
#'                              power=1)
#' #running the mcmc  , generations will be ignored because auto_stop=true
#' my_dpp_analysis$run(generations=1000,auto_stop=TRUE,max_gen = 10000,min_ess = 500)
#'
#' #we get rid of the first 25% of the output (burn-in)
#' hist(my_dpp_analysis$getNumCategoryTrace(0.25))
#'
#' my_dpp_analysis$getNumCategoryProbabilities(0.25)
#'}
#'


# generate documentation with roxygen2::roxygenise() or devtools::document(pkg=".")

dppMCMC_C <- setRefClass(

  Class = "dppMCMC_C",

  field = c(
    "dpp_mcmc_object"

  ), # end fields

  methods = list(
       ###############
       # Constructor #
       ###############

    initialize = function(data,
                          output,
                          model,
                          num_auxiliary_tables = 4,
                          expected_k = 2,
                          power = 1,
                          verbose=TRUE,
                          sample.num.clusters=TRUE) {

      "the class constructor, initializes DPPmcmc object with data and parameters"  # a docstring for documentation
       #initialize DPPmcmc object with data and parameters
       dpp_mcmc_object <<- new(DPPmcmc,data, model,num_auxiliary_tables,expected_k,power,effectiveSize,pmin)
       dpp_mcmc_object$setOutputPrefix(output)
       dpp_mcmc_object$setVerbose(verbose)
       dpp_mcmc_object$setSampleNumClusters(sample.num.clusters)
       dpp_mcmc_object$postInitialization()
       dpp_mcmc_object$makeOutputFiles()

    }, # end initialize

    #####################
    ###### Methods ######
    #####################



    run = function(generations,
                   sample_freq = generations/1000,
                   log_file,
                   allocation_file,
                   param_file,
                   append = TRUE,
                   random = FALSE,
                   auto_stop = FALSE,
                   min_ess = 500,
                   max_gen = 1e5) {
       "starts the MCMC run" # a docstring for documntation
       dpp_mcmc_object$run(generations,auto_stop,max_gen,min_ess,random,sample_freq)
    },

    getNumCategoryTrace = function(burnin_cutoff=0.25) {
      "returns the trace vector for the inferred number of categories in the data"
      catTrace<-dpp_mcmc_object$getNumCategoryTrace()
      catTrace<-catTrace[!is.na(catTrace)]
      #getting rid of 25% first results Burnin or whaever the cutoff is set for
      catTrace <-catTrace[c(ceiling(length(catTrace) * burnin_cutoff):length(catTrace))]
      return (catTrace)

    },

    getNumCategoryProbabilities = function(burnin_cutoff=0.25) {
      "returns the probabilities vector for inferred number of categories"
      catTrace<-getNumCategoryTrace(burnin_cutoff)
      category_probabilities <- table(catTrace)/length(catTrace)
      probs<-c()
      for(catIndex in 1:max(catTrace)){
        probs<-c(probs,category_probabilities[as.character(catIndex)])
      }
      probs[is.na(probs)]<-0
      names(probs)<-c(1:max(catTrace))
      return (probs)
    }


  ) # end methods

)
