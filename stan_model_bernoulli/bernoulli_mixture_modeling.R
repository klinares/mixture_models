pacman::p_load(bayesplot, BayesLCA, parallel, cmdstanr, poLCA, rstan, label.switching, data.table, tidyverse)

num_cores <- detectCores()
options(mc.cores = num_cores)



# read in data

outcome_variables <-  c("ABINSPAY", "ABMEDGOV1", "ABHELP1", "ABHELP2", "ABHELP3", "ABHELP4", 
                        "ABMORAL", "ABSTATE1", "ABSTATE2")

survey_data <- read_csv("/home/linarek/mixture_models/stan_model_bernoulli/model_data.csv") 


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
  df <- survey_data %>% select(outcome_variables)
  blca_fun(df, classes, 5000, 1000)},
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


# using polca
polca_model <- poLCA(
  cbind(ABINSPAY, ABMEDGOV1, ABHELP1, ABHELP2, ABHELP3, ABHELP4, ABMORAL, ABSTATE1, ABSTATE2)~1,
      nclass=3, 
      maxiter=1000,
  # recode variables, polca does not take "0"s
  data = survey_data %>% select(outcome_variables) %>% 
    mutate_at(outcome_variables, 
              list(~recode(., '0'=1, '1'=2))) 
)

#################
# 2. Stan model #

# install this for the first time to use cmdstanr
# install_cmdstan(cores = 12)

# compile model
mod <- cmdstan_model("/home/linarek/mixture_models/stan_model_bernoulli/bernoulli_mixture_model_test.stan")

model_data_list <-  list(y = items_only_data, #response matrix
                    J = nrow(items_only_data), #number of units/persons
                    I = ncol(items_only_data), #number of items
                    C = 2) # of classes to model

stan_fit <-
  mod$sample(
    data = model_data_list,
    chains = num_cores, 
    parallel_chains = num_cores,
    iter_warmup = 35,
    iter_sampling = 10,
    refresh = 100  )


# also save the model fit
stan_fit$save_object(file = "/home/linarek/Documents/stan models/stan_fit.rds")

# print model results
stan_fit$summary(c("nu", "probs"))

# plot trace plots
mcmc_trace(stan_fit$draws("nu"))

# plot area plot
mcmc_areas(stan_fit$draws("nu"))







# fix label switching

# variational bayes
fit_vb <- mod$variational(
  data = model_data_list, 
  iter = 15000,
  elbo_samples = 1000,
  algorithm = c("fullrank"),
  output_samples = 10000,
  tol_rel_obj = 0.00001)

# in rstan
# stan_vb <- stan_model("/home/linarek/mixture_models/stan_model_bernoulli/bernoulli_mixture_model.stan")

# vb_fit <- vb(
#     stan_vb,
#     data = data,
#     iter = 15000,
#     elbo_samples = 1000,
#     algorithm = c("fullrank"),
#     # imporance_resampling = T,
#     output_samples = 10000,
#     tol_rel_obj = 0.00001
#   )

# vb estimates
# print(vb_fit, c("alpha", "p"))




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
  data.table::transpose() %>% 
  rename(class_prediction=V1)

# merge class prediction to survey data
survey_data <- survey_data %>% add_column(class_prediction) 

# save combined data
write_rds(survey_data,"/home/linarek/mixture_models/model_results/survey_data_class_predictions.rds")

