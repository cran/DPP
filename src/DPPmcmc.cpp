#include <Rcpp.h>
#include <limits.h>
#include "DPPmcmc.h"
#include "helper.h"
#include <fstream>
#include <sstream>
#include "math.h"
using namespace Rcpp;

DPPmcmc::DPPmcmc(
    DoubleVector data_,
    Model& model_,
    int num_auxiliary_tables_,
    double expected_k_,
    int power_,
    Function effectiveSizeFunction_,
    Function pminFunction_

  ) :    data(data_),
         model(model_),
         effectiveSizeFunction(effectiveSizeFunction_),
         pminFunction(pminFunction_)
         //textBarFunction(textBarFunction_)

{
   RNGScope rngScope;
   power=power_;
   data=data_;
   num_auxiliary_tables=num_auxiliary_tables_;
   num_elements=data_.size();

   //Compute the alpha parameter
   //with prior expected number of
  //categories equal to expected_k

   estimate_concentration_parameter=model.getEstimateConcentrationParameter();
     if (estimate_concentration_parameter) {
       concentration_parameter_alpha= model.getConcentrationParameterAlpha();
       double concentration_parameter_mean= concentrationParameterFromK(num_elements,expected_k_);
       concentration_parameter_beta= concentration_parameter_alpha / concentration_parameter_mean;
       concentration_parameter=::Rf_rgamma(concentration_parameter_alpha,1/concentration_parameter_beta);
     } else {
       concentration_parameter=concentrationParameterFromK(num_elements,expected_k_);
     }

     //the following are initialized as TRUE but can be updated after object init.
     verbose=TRUE;
     sample_num_clusters=TRUE;

     //some code here was moved to the postInitialization method
     //and will be called after some values are sent from the R wrapper class.
}


void DPPmcmc::postInitialization(void) {
  // randomly generate starting values


  RNGScope rngScope;
  if(sample_num_clusters){
    allocation_vector=simulateChineseRestaurant(num_elements,concentration_parameter);
  } else {  //initialize with fix vector of 1
    allocation_vector=intRep(1,num_elements);
  }



  num_categories=max(allocation_vector);
  num_elements_in_each_category=Rcpp::table(allocation_vector);

  param_vector=model.base_distn_sim(num_categories);

/*
  //removing stochasticity
  DoubleVector means(1);
  means[0]=0.7;
  DoubleVector gamma_sds(1);
  gamma_sds[0]=0.3;

  Rcpp::List new_params=Rcpp::List::create(Named("means")=means,Named("sds")=gamma_sds);
  param_vector=new_params;
  //end removing stochasticity
*/

  /*
   //  TODO For debugging only, can be removed later
   std::ostringstream strs;
  strs << "allocation_vector:\t";
  strs << allocation_vector;
  strs<<"\n";
  strs << "max allocation_vector:\t";
  strs << max(allocation_vector);
  strs<<"\n";
  strs << "param_vector:\n";

  DoubleVector param1=param_vector(0);
  DoubleVector param2=param_vector(1);

  strs << "means:\t";

    for (int i=0;i<param1.size();i++){
      strs <<"\t"<<param1[i];
    }
  strs << "\nsds:\t";

    for (int i=0;i<param2.size();i++){
      strs <<"\t"<<param2[i];
    }


  strs<<"\n";
  write_text_to_log_file(strs.str());
*/

  num_params=param_vector.length();

  std::vector<double> likelihood_vector=model.likelihood_fn(data,allocation_vector,param_vector,power);
  likelihood=sumVector(likelihood_vector);

  std::vector<double> base_dist=makeDoubleVectorStandardDoubleVector(model.base_distn(param_vector));
  prior=sumVector(base_dist);

  generation=0;
  min_ESS=0;
}




double DPPmcmc::expectedNumberOfClusters(int num_elements_, double alpha_) {

  double e = 0.0;
  for(int i = 0; i < num_elements_; ++i) {
    e += alpha_ / (alpha_ + static_cast<double>(i));
  }
  /*  expectedNumberOfClusters = function(n,a) sum(a / (a + 1:n - 1))*/
  return e;

}


double DPPmcmc::concentrationParameterFromK(int num_elements_, double expected_num_clusters) {

  double alpha = 1.0;
  double e;
  double min_alpha = 0.0;
  //double max_alpha = std::numeric_limits<double>::max();
  double max_alpha = 1000;

  while ( true ) {
    e=expectedNumberOfClusters(num_elements_, alpha);
    if ( std::fabs(e- expected_num_clusters) < 1e-6 ) break;

    if ( e > expected_num_clusters ) {
      max_alpha = alpha;
      alpha = min_alpha + (alpha - min_alpha) / 2.0;
    } else {
      min_alpha = alpha;
      alpha = min_alpha + (max_alpha - alpha) / 2.0;
    }

  }

  return alpha;

}

IntegerVector DPPmcmc::simulateChineseRestaurant(int num_elements_, double alpha_) {

  IntegerVector res(num_elements_);
  RNGScope rngScope;
  int current_table = 0;

  for (int i = 0; i < num_elements_; ++i) {
    double new_cat_prob = alpha_ / (alpha_ + static_cast<double>(i));
    double u = runif(1,0,1)[0];
    if ( u < new_cat_prob ) {
      res[i] = current_table++;
    } else {
      int old_category_with_element = static_cast<int>(R::runif(0,i));
      res[i] = res[old_category_with_element];
    }
  }
  return (res+1);
}



void DPPmcmc::setOutputPrefix(std::string outputPrefix_) {
  outputPrefix=outputPrefix_;
}

void DPPmcmc::setVerbose(bool verbose_) {
  verbose=verbose_;
}

void DPPmcmc::setSampleNumClusters(bool sample_num_clusters_) {
  sample_num_clusters=sample_num_clusters_;

}


void DPPmcmc::concentrationParameterProposal() {

  // Following Escobar and West 1995:

  // First, sample the value eta.
  // eta ~ beta(concentration_parameter + 1, num_elements)

  RNGScope scope;
  double eta = ::Rf_rbeta(concentration_parameter + 1, num_elements);

  // Next, sample a new value of alpha.
  // alpha ~ pi * Gamma(shape = alpha + num_tables, rate = beta - log(eta)) &
  //         (1 - pi) * Gamma(shape= alpha + num_tables + 1, rate = beta - log(eta))

  double pi = (concentration_parameter_alpha + num_categories - 1) / ( num_elements * (concentration_parameter_beta - std::log(eta)) );
  double rate = concentration_parameter_beta - std::log(eta);

  double u = runif(1,0,1)[0];
  if ( (u / (1 - u)) < pi ) {
    concentration_parameter = ::Rf_rgamma(concentration_parameter_alpha + num_categories, 1 / rate);
  } else {
    concentration_parameter = ::Rf_rgamma(concentration_parameter_alpha + num_categories - 1 , 1 / rate);
  }

}




void  DPPmcmc::run(int generations,bool auto_stop,int max_gen,double min_ess,bool random,int sample_freq){
    //the following definition of rngScope is necesary
    //for propor behavioyr of the calls ro R::unif
    RNGScope rngScope;
    bool append=true;
// Set up the output file
    if (!append | auto_stop) {
      generation=0;
      makeOutputFiles();
    }

  int num_logged=0;
  return_num_cats_trace=std::vector<int>(1);

  std::vector<double> likelihood_trace(1000,NA_REAL);
  std::vector<int> num_cats_trace(1000,NA_INTEGER);
  std::vector<double> alpha_trace(1000,NA_REAL);

  int init_generations=generation;

  sample_freq=std::max(1,sample_freq);

  int i=0;
  Normal normal_distribution(0,1);
  int k=0;

  while(TRUE) {

    i = i + 1;
    generation=generation + 1;
    /*
    //  TODO For debugging only, can be removed later
    std::ostringstream strs;
    strs << "before allocation proposal:\t";
    strs << i;
    strs << "likelihood:\t";
    strs << likelihood;
    strs << "allocation:\t";
    strs << allocation_vector;
    strs << "\n";
    //strs<<"\t prior_theta:\t";
    // strs<<prior_theta;
    // strs<<"\t param_vector:\t";
    //  strs<<param_theta;
    write_text_to_log_file(strs.str());
    //end of TODO
*/

    //do allocation proposal


    if(sample_num_clusters){
      if ( random ) {
       k=normal_distribution.sample_int(num_elements);
       allocationProposal(k);
      } else {
      for(k=1;k<=num_elements;k++) allocationProposal(k);
      }
    }



    //update the likelihood
    std::vector<double> likelihood_vector=model.likelihood_fn(data,allocation_vector,param_vector,power);
    likelihood=sumVector(likelihood_vector);

    std::vector<double> base_dist= makeDoubleVectorStandardDoubleVector(model.base_distn(param_vector));
    prior=sumVector(base_dist);

/*
     //  TODO For debugging only, can be removed later
     std::ostringstream strs2;
    strs2 << "iteration:\t";
    strs2 << i;
    strs2 << "likelihood:\t";
    strs2 << likelihood;
    strs2 << "\nallocation:\t";
    strs2 << allocation_vector;
    strs2 << "\n";
    //strs<<"\t prior_theta:\t";
    // strs<<prior_theta;
    // strs<<"\t param_vector:\t";
    //  strs<<param_theta;
    write_text_to_log_file(strs2.str());
    //end of TODO
*/

    Uniform uniform_dist=Uniform(0,1);


    for(k=1;k<=num_categories;k++){

      Rcpp::List param_theta=model.proposal_distn(param_vector);
      double likelihood_theta=0;
      std::vector<double> likelihood_theta_vector=model.likelihood_fn(data,allocation_vector,param_theta,power);
      likelihood_theta=sumVector(likelihood_theta_vector);

      double prior_theta=0;
      std::vector<double> prior_theta_vector= makeDoubleVectorStandardDoubleVector(model.base_distn(param_theta));
      prior_theta=sumVector(prior_theta_vector);

      double R=std::exp(likelihood_theta-likelihood+prior_theta-prior);
      //double u=uniform_dist.sample(1)[0];
      double u=R::runif(0.0,1.0);

/*
         //  TODO For debugging only, can be removed later
       std::ostringstream strs;
      strs << "likelihood:\t";
      strs << likelihood;
      strs << "param_theta_0:\t";
      strs << makeDoubleVectorStandardDoubleVector(param_theta(0))[0];
      strs << "param_theta_1:\t";
      strs << makeDoubleVectorStandardDoubleVector(param_theta(1))[0];
      strs<<"\t prior_theta:\t";
      strs<<prior_theta;
      strs<<"\t u:\t";
      strs<<u;
      strs<<"\t R:\t";
      strs<<R;
      strs<<"\t likelihood_theta:\t";
      strs<<likelihood_theta;
      strs<<"\t prior:\t";
      strs<<prior;
      strs<<"\nsuperSuma:\t";
      strs<<(likelihood_theta-likelihood+prior_theta-prior);
      strs<<i<<"\t"<<R<<"\t"<<u<<"\t"<<makeDoubleVectorStandardDoubleVector(param_theta(0))[0];


      // strs<<"\t param_vector:\t";
      //  strs<<param_theta;
      write_text_to_log_file(strs.str());
      //end TODO
*/
      if(u < R) {
        likelihood=likelihood_theta;
        prior=prior_theta;
        param_vector=param_theta;
  /*    //  TODO For debugging only, can be removed later
         std::ostringstream strs;
        strs << "likelihood:\t";
        strs << likelihood;
        strs<<"\t prior_theta:\t";
        strs<<prior_theta;
       // strs<<"\t param_vector:\t";
      //  strs<<param_theta;
        write_text_to_log_file(strs.str());
        //end TODO
        */

      }

   }


    if (estimate_concentration_parameter ) {
     //  do concentration parameter proposal
       concentrationParameterProposal();
    }


   if ((i % sample_freq) == 0) writeOutputFiles();


    // check for MCMC stop
    if (auto_stop) {
      // auto-stopping rule
      // keep track of the things to sample
      if (i % sample_freq == 0) {
        num_logged=num_logged + 1;
        likelihood_trace[num_logged]=likelihood;
        num_cats_trace[num_logged]=num_categories;
        alpha_trace[num_logged]= concentration_parameter;
      }

      // extend the trace if necessary
      if ( (num_logged + 1) % 1000 == 0 ) {
        std::vector<double> emptyDoubleVector(1000,NA_REAL);
        std::vector<int> emptyIntegerVector(1000,NA_INTEGER);
        for(int i=0;i<1000;i++){
         emptyDoubleVector[i]=NA_REAL;//nan,not a numer
        }
        likelihood_trace=concatenateVectors(likelihood_trace,emptyDoubleVector);
        num_cats_trace=concatenateVectors(num_cats_trace,emptyIntegerVector);
        alpha_trace=concatenateVectors(alpha_trace,emptyDoubleVector);
      }

      // compute the minimum ESS using a 25% burnin
      if (num_logged >= 10) {
        //DoubleVector likelihood_ESS_vector=effectiveSizeFunction(elementsInRange((int)floor(num_logged / 4),num_logged-1,likelihood_trace));
       // double likelihood_ESS=(double)likelihood_ESS_vector[0];

        DoubleVector num_cats_ESS_vector=effectiveSizeFunction(elementsInRange((int)std::floor(num_logged / 4.0),num_logged-1,num_cats_trace));
        double num_cats_ESS=num_cats_ESS_vector[0];

        NumericVector pminResult=pminFunction(num_cats_ESS);
        min_ESS = (double)pminResult[0];

        double  alpha_ESS=0;
        if (estimate_concentration_parameter) {
         NumericVector alpha_ESS_vector=effectiveSizeFunction(elementsInRange((int)std::floor(num_logged / 4.0),num_logged-1,alpha_trace));
          alpha_ESS= alpha_ESS_vector[0];
          NumericVector pminResult=pminFunction(alpha_ESS,num_cats_ESS);
          min_ESS = (double)pminResult[0];
        }


       //  TODO For debugging only, can be removed later
   /*    std::ostringstream strs;
       strs << "num_logged:\t";
       strs << num_logged;
       strs << "min_ESS:\t";
       strs << min_ESS;
       strs<<"\t min_ess:\t";
       strs<<min_ess;
       // strs<<"\t param_vector:\t";
       //  strs<<param_theta;
       write_text_to_log_file(strs.str());*/

        if ((min_ESS > min_ess) || (i >= max_gen))  {
          return_num_cats_trace=num_cats_trace;
          break;
        }

     }

    } else { //if not auto_stop
      // traditional stopping rule
      if (generation == (init_generations + generations)) break;
    }

  } //from while(true)

}

void DPPmcmc::makeOutputFiles(void){
  std::string mcmcLogName=outputPrefix+"mcmc.log";
  std::ofstream of(mcmcLogName.c_str());

  if(of.is_open())
  {
    of << "generation\tlikelihood\tnum_categories\tconc\tmin_ESS"<< std::endl;
    of.flush();
    of.close();//std::cout<<"wrote the file successfully!"<<std::endl;
  }
  else
  {
    Rcerr<<"Failed to open file : "<<mcmcLogName<<std::endl; // return -1;
   }

  if(verbose){

    /* allocation log file */
    std::string allocationOutputName=outputPrefix+"allocation.log";
    std::ofstream of2(allocationOutputName.c_str());

    if(of2.is_open()){
      of2 << "generation\t";
      for (int i=1;i<=num_elements;i++) {of2 <<"x_"<<i<<"\t";}
      of2<<std::endl;
      of2.flush();of2.flush();
    } else {
      Rcerr<<"Failed to open file : "<<allocationOutputName<<std::endl; // return -1;
    }

    /* parameters log files */
    std::string param1OutputName=outputPrefix+"param_1.log";
    std::ofstream of3(param1OutputName.c_str());
    if(of3.is_open()) {
      of3 << "generation\tparameters"<<std::endl;
      of3.flush();of3.flush();
    }  else {
      Rcerr<<"Failed to open file : "<<param1OutputName<<std::endl; // return -1;
    }
    std::string param2OutputName=outputPrefix+"param_2.log";
    std::ofstream of4(param2OutputName.c_str());
    if(of4.is_open()) {
      of4 << "generation\tparameters"<<std::endl;
      of4.flush();of4.flush();
    } else {
      Rcerr<<"Failed to open file : "<<param2OutputName<<std::endl; // return -1;
    }

    /* allocation params log files */
    std::string param1AllocOutputName=outputPrefix+"allocation_param_1.log";
    std::ofstream of5(param1AllocOutputName.c_str());
    if(of5.is_open()) {
      of5 << "generation\t";
      for (int i=1;i<=num_elements;i++) {of5 <<"x_"<<i<<"\t";}
      of5<<std::endl;
      of5.flush();of5.flush();
    }  else  {
      Rcerr<<"Failed to open file : "<<param1AllocOutputName<<std::endl; // return -1;
    }

    std::string param2AllocOutputName=outputPrefix+"allocation_param_2.log";
    std::ofstream of6(param2AllocOutputName.c_str());
    if(of6.is_open()) {
      of6 << "generation\t";
      for (int i=1;i<=num_elements;i++) {of6 <<"x_"<<i<<"\t";}
      of6<<std::endl;
      of6.flush();of6.flush();
    }  else {
      Rcerr<<"Failed to open file : "<<param2AllocOutputName<<std::endl; // return -1;
    }
  } //end of if (verbose)
  writeOutputFiles();

}

void DPPmcmc::writeOutputFiles(void){

    /* updating mcmc log*/
    std::ofstream of;
    of.open( (outputPrefix+"mcmc.log").c_str(), std::ios_base::out | std::ios_base::app );
    if(of.is_open())
    {
      //of << generation<<"\t"<<likelihood<<"\t"<<num_categories<<"\t"<<num_elements_in_each_category<<"\t"<<concentration_parameter<<std::endl;
      of << generation<<"\t"<<likelihood<<"\t"<<num_categories<<"\t"<<concentration_parameter<<"\t"<<min_ESS<<std::endl;
      of.flush();
      of.close();  //TODO Luis, perhaps I could close the files at the end to gain speed
    }

    if (verbose){

      /* updating allocation log file */
      std::ofstream of2;
      of2.open( (outputPrefix+"allocation.log").c_str(), std::ios_base::out | std::ios_base::app );
      if(of2.is_open())
      {
        of2 <<generation<<"\t";
        of2 <<allocation_vector<<std::endl;
        of2.flush();of2.flush();
      }

      /* updating parameter ouput files */
      DoubleVector param1=param_vector(0);
      DoubleVector param2=param_vector(1);
      std::ofstream of3;
      of3.open( (outputPrefix+"param_1.log").c_str(), std::ios_base::out | std::ios_base::app );
      if(of3.is_open()){
        of3 <<generation;
        for (int i=0;i<param1.size();i++){
          of3 <<"\t"<<param1[i];
        }
        of3 <<std::endl;
        of3.flush();of3.flush();
      }

      std::ofstream of4;
      of4.open( (outputPrefix+"param_2.log").c_str(), std::ios_base::out | std::ios_base::app );
      if(of4.is_open()){
        of4 <<generation;
        for (int i=0;i<param2.size();i++){
          of4 <<"\t"<<param2[i];
        }
        of4   <<std::endl;
        of4.flush();of4.flush();
      }

      /* updating allocation parameter ouput files */

      DoubleVector allocatedParam1(data.size());
      DoubleVector allocatedParam2(data.size());

      for (int i=0;i<data.size();i++){
        allocatedParam1[i]=param1[allocation_vector[i]-1];
        allocatedParam2[i]=param2[allocation_vector[i]-1];
      }

      std::ofstream of5;
      of5.open( (outputPrefix+"allocation_param_1.log").c_str(), std::ios_base::out | std::ios_base::app );
      if(of5.is_open()){
        of5 <<generation<<"\t"<<allocatedParam1<<std::endl;
        of5.flush();of5.flush();
      }

      std::ofstream of6;
      of6.open( (outputPrefix+"allocation_param_2.log").c_str(), std::ios_base::out | std::ios_base::app );
      if(of6.is_open())  {
        of6 <<generation<<"\t"<<allocatedParam2<<std::endl;
        of6.flush();of6.flush();
      }

    }//end of if (verbose)
}



void DPPmcmc::allocationProposal(int element){
/*  TODO For debugging only, can be removed later
  std::ostringstream strs;
  strs << "num_elements_in_each_category:\t";
  strs << num_elements_in_each_category;
  strs<<"\t num_categories:\t";
  strs<<num_categories;
  write_text_to_log_file(strs.str());
*/

  // get the current category,
  // decrement the number of elements there
  int current_category=allocation_vector[element-1];
  num_elements_in_each_category[current_category-1]=num_elements_in_each_category[current_category-1]-1;




  // make the rate parameters
  Rcpp::List new_params =model.base_distn_sim(num_auxiliary_tables);
  //compute the likelihood for the
  //categories that exist.







  std::vector<double> old_category_likelihoods=model.single_likelihood_fn(data[element-1],
                                                                          integerSequence(1,num_categories),
                                                                          param_vector,
                                                                          power);



  std::vector<double> new_category_likelihoods=model.single_likelihood_fn(data[element-1],
                                                                          integerSequence(1,num_auxiliary_tables),
                                                                          new_params,
                                                                          power);
 /*
  //  TODO For debugging only, can be removed later
  std::ostringstream strs;
  strs << "\nnum_elements_in_each_category:\t";
  strs << num_elements_in_each_category;
  strs<<"\t old_category_likelihoods:\t";
  for(int i=0;i<old_category_likelihoods.size();i++) strs<<old_category_likelihoods[i]<<"\t";
   strs<<" ## end old category likelihoods \n";
   strs<<"\t new_category_likelihoods:\t";
   for(int i=0;i<new_category_likelihoods.size();i++) strs<<new_category_likelihoods[i]<<"\t";
   strs<<" ## end  mew category likelihoods \n";
  write_text_to_log_file(strs.str());
   for (int z=0;z<new_category_likelihoods.size();z++){
     if(std::isnan(new_category_likelihoods[z])) {
       std::ostringstream strs3;
       strs3 << "\n new cat.like. dio nan. data:\t";
       strs3 << data[element-1];
       strs3 << "\tsequence [1:4]";
       strs3 << "\new params:\t";
       std::vector<double> new_params1=new_params(0);
       std::vector<double> new_params2=new_params(1);
       for(int i=0;i<new_params1.size();i++) strs3<<new_params1[i]<<"\t";

        strs3 << "\t";
        for(int i=0;i<new_params2.size();i++) strs3<<new_params2[i]<<"\t";

       strs3 << "\n";
       write_text_to_log_file(strs3.str());
     }
   }

  // TODO end
*/




  //category_likelihoods= old_category_likelihoods.insert(old_category_likelihoods.end(),new_category_likelihoods.begin(),new_category_likelihoods.end());
   std::vector<double> category_likelihoods=concatenateVectors(old_category_likelihoods,new_category_likelihoods);

   /*
   //for debugging
   std::ostringstream strs;
   strs<<"\t Rcpp_category_likelihoods:\t";
   for(int i=0;i<category_likelihoods.size();i++) strs<<category_likelihoods[i]<<"\t";
   strs<<" ## end old category likelihoods \n";
   write_text_to_log_file(strs.str());
   //end for debugging
  */
  //compute the category priors
   std::vector<double> existing_priors=divideIntegerVectorByDouble(num_elements_in_each_category, (concentration_parameter + num_elements - 1));
     std::vector<double> augmented_priors=
   divideVectorByDouble(rep(concentration_parameter,num_auxiliary_tables), (num_auxiliary_tables * (concentration_parameter + num_elements - 1)));
   std::vector<double> priors=logVector(concatenateVectors(existing_priors,augmented_priors));

  //choose the new rate category
  std::vector<double> category_probs=expVector(category_likelihoods+priors);

  /*
  //  TODO For debugging only, can be removed later
   std::ostringstream strs2;
  strs2 << "num_elements_in_each_category:\t";
  //strs2 << num_elements_in_each_category;
  for(int i=0;i<num_elements_in_each_category.size();i++) strs2<<num_elements_in_each_category[i]<<"\t";
  strs2<<"\t category_likelihoods:\t";
  for(int i=0;i<category_likelihoods.size();i++) strs2<<category_likelihoods[i]<<"\t";
  //strs<<"\t category_likelihoo\n";
  strs2<<"\t category_probs:\t";
  for(int i=0;i<category_probs.size();i++) strs2<<category_probs[i]<<"\t";
  //strs<<"\n";
  write_text_to_log_file(strs2.str());

  // TODO end
*/


  //new_category <- sample.int(num_categories + num_auxilliary_tables,size=1,prob=category_probs)
  Normal normal_dist=Normal(0,1); //this parameters are irrelevant, consider moving the
                                  //methode sample_int_prob somewhere else
  int new_category = normal_dist.sample_int_prob(category_probs);


/*
  //for debuging TODO remove
  DoubleVector probs(category_probs.size());

  for (int i=0;i<probs.size();i++){
    probs[i]=category_probs[i];

  }

  DoubleVector oldCats(old_category_likelihoods.size());

  for (int i=0;i<oldCats.size();i++){
    oldCats[i]=old_category_likelihoods[i];

  }

  //  TODO For debugging only, can be removed later
  std::ostringstream strs;
  strs << "old cat likes:\t";
  strs << oldCats<<"\n";
  strs << "category_probs:\t";
  strs <<probs<<"\n";
  write_text_to_log_file(strs.str());
*/


  // update the parameters
  if(new_category != current_category) {
   // write_text_to_log_file("new cat  != current cat");
    if(new_category > num_categories) {
     // write_text_to_log_file("new cat  > num_cat");
     // writeOutputFiles();
      //add the new parameters to the existing parameter vector
      //hardcoded for two parameters of the normal_model
      //TODO make more general

        std::vector<double> old_params1=param_vector(0);
        std::vector<double> old_params2=param_vector(1);
        std::vector<double> new_params1=new_params(0);
        std::vector<double> new_params2=new_params(1);

        std::vector<double> really_new_params1(1);
        std::vector<double> really_new_params2(1);
        really_new_params1[0]=new_params1[new_category - num_categories-1];
        really_new_params2[0]=new_params2[new_category - num_categories-1];
        std::vector<double> return_params1;
        std::vector<double> return_params2;
       // write_text_to_log_file("new cat  > num_cat par2");
        return_params1=concatenateVectors(old_params1,really_new_params1);
        return_params2=concatenateVectors(old_params2,really_new_params2);

        param_vector=Rcpp::List::create(Named("means")=return_params1,Named("sds")=return_params2);

        allocation_vector[element-1]=num_categories+1;
        num_elements_in_each_category.push_back(1);
        //num_elements_in_each_category[num_categories]=1;
        num_categories=num_categories+1;
        //write_text_to_log_file("new cat  > num_cat  part 3");
        //write_text_to_log_file("a new category was added");

    } else {
      //write_text_to_log_file("new cat  < num_cat");
       allocation_vector[element-1] = new_category;
       num_elements_in_each_category[new_category-1] = num_elements_in_each_category[new_category-1] + 1;
    }
  }else{ //write_text_to_log_file("new cat  = current cat");
    num_elements_in_each_category[current_category-1]= num_elements_in_each_category[current_category-1] + 1;
  }



  if (anyEqual(makeIntegerVectorStandardDoubleVector(num_elements_in_each_category),0.0)){

    //num
      //TODO  comment
     //result is zero based for now, potential problems

     std::vector<double> these_categories_are_empty=whichAreEqual(makeIntegerVectorStandardDoubleVector(num_elements_in_each_category),0.0);
     //these_categories_are_empty <- which(num_elements_in_each_category == 0)
    // write_text_to_log_file("made to stage #1");
     //writeOutputFiles();
     for(int i=0;i<these_categories_are_empty.size();i++) {
      // write_text_to_log_file("made to stage #2");
       //TODO Luis comment
       //again, hardcoded for two parametersof the normal distribution, mean and sd

       int c=these_categories_are_empty[i];

       /*
       for(c in these_categories_are_empty) {
         for(i in 1:num_params) {
           param_vector[[i]] <<- param_vector[[i]][-c]
         }
         allocation_vector <<- allocation_vector - (allocation_vector > c)
           num_categories <<- num_categories - 1
         num_elements_in_each_category <<- num_elements_in_each_category[-c]
       }*/

       std::vector<double> old_params1=param_vector(0);
       std::vector<double> old_params2=param_vector(1);

       old_params1=removeElementAtPosition(old_params1,c);
       old_params2=removeElementAtPosition(old_params2,c);

       //write_text_to_log_file("made to stage #6");

       param_vector=Rcpp::List::create(Named("means")=old_params1,Named("sds")=old_params2);

        allocation_vector=makeIntegerVectorStandardIntVector(allocation_vector)-
        evaluateVectorGreaterThanInt(makeIntegerVectorStandardIntVector(allocation_vector),c+1);


        num_categories=num_categories-1;
        num_elements_in_each_category=removeElementAtPosition(makeIntegerVectorStandardIntVector(num_elements_in_each_category),c);
      // num_elements_in_each_category.erase(num_elements_in_each_category.begin()+c);
      //write_text_to_log_file("category removed");

     }
  }

  /*
  //  TODO For debugging only, can be removed later
  std::ostringstream strs8;
  strs8 << "\nnum_elements_in_each_category:\t";
  strs8 << num_elements_in_each_category<<"\n";
 strs8 <<"\n";
  write_text_to_log_file(strs8.str());
  // end TODO
  */
}


std::vector<int> DPPmcmc::makeIntegerVectorStandardIntVector(IntegerVector vector1){
  std::vector<int> returnVector(vector1.size());
  for(int i=0;i<vector1.size();i++){
    returnVector[i]=vector1[i];
  }
  return returnVector;
}

std::vector<double> DPPmcmc::makeIntegerVectorStandardDoubleVector(IntegerVector vector1){
  std::vector<double> returnVector(vector1.size());
  for(int i=0;i<vector1.size();i++){
    returnVector[i]=vector1[i];
  }
  return returnVector;
}


std::vector<double> DPPmcmc::makeDoubleVectorStandardDoubleVector(DoubleVector vector1){
  std::vector<double> returnVector(vector1.size());
  for(int i=0;i<vector1.size();i++){
    returnVector[i]=vector1[i];
  }
  return returnVector;
}



std::vector<double> DPPmcmc::dummyFunction(std::vector<double> dummyInput){
/* this fuction is used for testing only TODO  */
/*res[0]=2;
res[1]=1;
res[2]=1;
res[3]=1;
res[4]=1;*/

 //std::vector<double> causa(100);
 //for(int i=0;i<100;i++) { causa[i]=::Rf_rgamma(dummyInput[0],dummyInput[1]); }
// causa= makeIntegerVectorStandardDoubleVector(integerSequence(1100,5000));

  Normal myNormal=Normal(0,0.5);
  Uniform uniform_dist=Uniform(0,1);
  //myNormal.sample_int(2);
  //std::vector<double> causa(100);
  //for(int i=0;i<100;i++) {causa[i]=myNormal.rnorm(dummyInput[0],dummyInput[1]); }

   std::vector<double> causa(100);
  // for(int i=0;i<100;i++) {causa[i]=myNormal.sample_int_prob(dummyInput); }
  for(int i=0;i<100;i++) {causa[i]=myNormal.sample_int(200); }
  //for(int i=0;i<100;i++) {causa[i]= uniform_dist.sample(1)[0]; }
  //for(int i=0;i<dummyInput.size();i++) {causa[i]= exp(dummyInput[i]); }

 /*
IntegerVector res(causa.size());
for (int i=0;i<causa.size();i++)
{
   res[i]=causa[i];
  res[i]=causa[i];
}*/

 //Normal normal_dist=Normal(0,1); //this parameters are irrelevant, consider moving the
 //methode sample_int_prob somewhere else

 // int new_category = normal_dist.sample_int_prob(dummyInput);

//return Rcpp::table(res);
//return  concatenateVectors(dummyInput,causa);
//  return rep((double)new_category,3);
// get the current category,
// decrement the number of elements there

// get the current category,
// decrement the number of elements there

  // NumericVector testor=effectiveSizeFunction(dummyInput);
  //NumericVector testor=pminFunction(dummyInput,3.14);
//  double testor2=testor[0]+0;
  // return rep((double)(5%2),2);
  return causa;

}

IntegerVector DPPmcmc::integerSequence(int min,int max){

  IntegerVector my_sequence(max-min+1);
  int counter=0;
  for (int i=min;i<=max;i++){
    my_sequence[counter]=i;
    counter++;
  }
  return my_sequence;
}





std::vector<int> DPPmcmc::evaluateVectorGreaterThanInt(std::vector<int> vector1,int val){
  std::vector<int> returnVector(vector1.size());
  for(int i=0;i<vector1.size();i++){
    int eval=0;
    if (vector1[i]>val) eval=1;
    returnVector[i]=eval;
  }
  return returnVector;
}


std::vector<double> DPPmcmc::divideIntegerVectorByDouble(IntegerVector vector1,double denominator){
  std::vector<double> returnVector(vector1.size());
  for(int i=0;i<vector1.size();i++){
    returnVector[i]=(vector1[i]/denominator);
  }
  return returnVector;

}


std::vector<double> DPPmcmc::rep(double a,int num_reps)
{
  std::vector<double> result(num_reps);
  for(int i=0;i<num_reps;i++)
  {
     result[i]=a;
  }
  return result;
}


IntegerVector DPPmcmc::intRep(int a,int num_reps)
{

  IntegerVector result(num_reps);
  for(int i=0;i<num_reps;i++)
  {
    result[i]=a;
  }
  return result;
}


void DPPmcmc::write_text_to_log_file( const std::string &text )
{
  std::ofstream log_file(
      "c_log_file.txt", std::ios_base::out | std::ios_base::app );
  log_file << text << std::endl;
}

/////////////////
// Rcpp Module //
/////////////////

RCPP_MODULE(DPPmcmc) {

  using namespace Rcpp;

  class_<DPPmcmc>( "DPPmcmc" )
    //.constructor<DoubleVector, NormalModel&,int, double, int,Function,Function,Function>()
      .constructor<DoubleVector, NormalModel&,int, double, int,Function,Function>()
    .method("simulateChineseRestaurant", &DPPmcmc::simulateChineseRestaurant)

    .method("expectedNumberOfClusters", &DPPmcmc::expectedNumberOfClusters)
    .method("concentrationParameterFromK", &DPPmcmc::concentrationParameterFromK)
    .method("getNumElements", &DPPmcmc::getNumElements)

    .method("getNumCategories", &DPPmcmc::getNumCategories)
    .method("getData", &DPPmcmc::getData)
    .method("getPower", &DPPmcmc::getPower)
    .method("getNumAuxiliaryTables", &DPPmcmc::getNumAuxiliaryTables)
    .method("getNumElementsPerTable", &DPPmcmc::getNumElementsPerTable)
    .method("getNumElementsInEachCategory", &DPPmcmc::getNumElementsInEachCategory)
    .method("getNumCategoryTrace", &DPPmcmc::getNumCategoryTrace)


    .method("getConcentrationParameter", &DPPmcmc::getConcentrationParameter)
    .method("getConcentrationParameterAlpha", &DPPmcmc::getConcentrationParameterAlpha)
    .method("getConcentrationParameterBeta", &DPPmcmc::getConcentrationParameterBeta)
    .method("getEstimateConcentrationParameter", &DPPmcmc::getEstimateConcentrationParameter)
    .method("getAllocationVector", &DPPmcmc::getAllocationVector)
    .method("getNumParams", &DPPmcmc::getNumParams,"place holder for getNumParams docs")

    .method("dummyFunction", &DPPmcmc::dummyFunction)
    .method("getParamVector", &DPPmcmc::getParamVector,"place holder for getParamVector docs")
    .method("integerSequence", &DPPmcmc::integerSequence)
   // .method("allocationProposal", &DPPmcmc::allocationProposal)
  //  .method("concentrationParameterProposal", &DPPmcmc::concentrationParameterProposal)
    .method("writeOutputFiles", &DPPmcmc::writeOutputFiles)
    .method("postInitialization", &DPPmcmc::postInitialization)
    .method("setOutputPrefix", &DPPmcmc::setOutputPrefix)
    .method("setVerbose", &DPPmcmc::setVerbose)
    .method("setSampleNumClusters", &DPPmcmc::setSampleNumClusters)
    .method("getOutputPrefix", &DPPmcmc::getOutputPrefix)
    .method("makeOutputFiles", &DPPmcmc::makeOutputFiles)
    .method("run", &DPPmcmc::run)

  ;

}
