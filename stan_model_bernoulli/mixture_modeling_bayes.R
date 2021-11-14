library(tidyverse)
library(sjlabelled)
library(bayesplot)
library(rstan)
library(BayesLCA)
library(parallel)
library(cmdstanr)
num_cores <- detectCores()
options(mc.cores = num_cores)


df <- read_spss("GSS2018.sav")


outcome_variables <- c("ABANY", "ABHELP1", "ABHELP2", "ABHELP3", "ABHELP4", "ABHLTH", "ABPOOR")

demogrpahic_variables <- c("SEX", "RACE", "CLASS", "PARTYID")

sample <- df %>% select(outcome_variables, demogrpahic_variables, ID) 

# print questions
sample %>% get_label()

# print value codes
sample %>% get_label()

# recode "no" "yes" to = 0, missing as "1"
model_data_bernoulli <- sample %>% 
  mutate_at(outcome_variables, 
            list(~recode(., '1'=0, '2'=0))) %>% 
  mutate_at(outcome_variables, 
            list( ~replace_na(., "1"))) %>% 
  drop_na() %>% 
  select(outcome_variables) %>% 
  mutate(across(where(is.character), as.numeric))

# print counts & proportion for each item
map(outcome_variables, function(x){
  model_data_bernoulli %>% count_(x) %>% 
    mutate(proportion = n/sum(n))
})

# save data with demographics
write_csv(sample %>% 
  mutate_at(outcome_variables, 
            list(~recode(., '1'=0, '2'=0))) %>% 
  mutate_at(outcome_variables, 
            list( ~replace_na(., "1"))) %>% 
  drop_na(), "model_data.csv")



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
  blca_fun(model_data_bernoulli, classes, 4000, 3000)},
  mc.cores = num_cores) 

# print BIC goodness of fit estimates for each model
map(blca_model_results_list, function(model){
  list(
  model$BICM,
  model$classprob,
  model$itemprob)
  })


# print a few plots on one of the models
plot(blca_model_results_list[[3]], which=1) # item probabilities
plot(blca_model_results_list[[3]], which=5) # trace plots check for class switching



#################
# 2. Stan model #

# install this for the first time to use cmdstanr
#install_cmdstan(cores = 12)

# compile model
mod <- cmdstan_model("bernoulli_mixture_model.stan")

stan_fit <-
  mod$sample(
    data = list(y = model_data_bernoulli, #response matrix
                J = nrow(model_data_bernoulli), #number of units/persons
                I = ncol(model_data_bernoulli), #number of items
                C = 2), # of classes to model
    chains = 12, 
    parallel_chains = 12,
    iter_warmup = 4000,
    iter_sampling = 4000,
    refresh = 500
  )


# print posterior class & item probabilities 
stan_fit$summary(c("alpha", "p")) %>% print(n=30)

draws_arr <- stan_fit$draws() # or format="array"
str(draws_arr)

