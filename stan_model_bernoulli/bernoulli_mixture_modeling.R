pacman::p_load(tidyverse, bayesplot, BayesLCA, parallel, cmdstanr, data.table)

num_cores <- detectCores()
options(mc.cores = num_cores)



# read in data

outcome_variables <-  c("ABINSPAY", "ABMEDGOV1", "ABHELP1", "ABHELP2", "ABHELP3", "ABHELP4", 
                        "ABMORAL", "ABSTATE1", "ABSTATE2")

survey_data <- read_csv("model_data.csv") 
items_only_data <- survey_data %>% select(outcome_variables)


#######################################
############### modeling ##############

# 1. BayesLCA package - Gibbs sampler

# function to iterate over multiple class models
blca_fun <- function(df, num_classes, iterations, burn_in){
  blca.gibbs(df, num_classes, iter = iterations, burn.in = burn_in)
}

# iterate over x number of classes using the gibbs sampler
## model for endorsing items
blca_model_results_list <- mclapply(seq(1,4), function(classes){
  blca_fun(items_only_data, classes, 6000, 4000)},
  mc.cores = num_cores) 

# print BIC goodness of fit estimates for each model & inspect
## lower BIC the better
## 1-class model is our baseline model to compare
map(blca_model_results_list, function(model){
  list(
    model$BICM,
    model$classprob,
    model$itemprob)
})


# print a few plots on one of the models
plot(blca_model_results_list[[2]], which=1) # item probabilities
plot(blca_model_results_list[[2]], which=5) # trace plots check for class switching



#################
# 2. Stan model #

# install this for the first time to use cmdstanr
# install_cmdstan(cores = 12)

# compile model
mod <- cmdstan_model("bernoulli_mixture_model.stan")

stan_fit <-
  mod$sample(
    data = list(y = items_only_data, #response matrix
                J = nrow(items_only_data), #number of units/persons
                I = ncol(items_only_data), #number of items
                C = 2), # of classes to model
    chains = num_cores, 
    parallel_chains = num_cores,
    iter_warmup = 20000,
    iter_sampling = 2000,
    refresh = 1000
  )

# print model results
stan_fit$summary(c("alpha", "p"))

# plot trace plots
mcmc_trace(stan_fit$draws("alpha"))


# function for mode
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


# save out class predictions
class_prediction <- stan_fit$draws(format = "df") %>% 
  select(starts_with("pred_class_dis")) %>% 
  map_dfr(., function(x) {(calculate_mode(x))}
          ) %>% 
  transpose() %>% 
  rename(class_prediction=V1)

# merge class prediction to survey data
survey_data <- survey_data %>% add_column(class_prediction) 

# save combined data
write_rds(survey_data,"survey_data_class_predictions.rds")

# also save the model fit
stan_fit$save_object(file = "stan_fit.rds")

