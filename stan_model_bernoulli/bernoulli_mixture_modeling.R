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
      nclass=2, 
      maxiter=1000,
  # recode variables, polca does not take "0"s
  data = survey_data %>% select(outcome_variables) %>% 
    mutate_at(outcome_variables, 
              list(~recode(., '0'=1, '1'=2))) 
)

#################
# 2. Stan model #
model_data_list <-  list(y = survey_data %>% select(outcome_variables), #response matrix
                         J = nrow(survey_data %>% select(outcome_variables)), #number of units/persons
                         I = ncol(survey_data %>% select(outcome_variables)), #number of items
                         C = 2) # of classes to model

###################### cmdstanr ###############################
# install this for the first time to use cmdstanr
# install_cmdstan(cores = 12)

# compile model
# mod <- cmdstan_model("/home/linarek/mixture_models/stan_model_bernoulli/bernoulli_mixture_model.stan")
 

# stan_fit <-
#   mod$sample(
#     data = model_data_list,
#     chains = num_cores, 
#     parallel_chains = num_cores,
#     iter_warmup = 2000,
#     iter_sampling = 200,
#     refresh = 500  )

# print model results
#stan_fit$summary(c("alpha", "p"))

# plot trace plots
#mcmc_trace(stan_fit$draws("alpha"))

# plot area plot
#mcmc_areas(stan_fit$draws("alpha"))

# save fit
#r_stan_fit$save_object(file = "/home/linarek/Documents/stan models/stan_fit.rds")


##############################################

model <- stan_model("/home/linarek/mixture_models/stan_model_bernoulli/bernoulli_mixture_model.stan")

# using rstan
r_stan_fit <-
  stan(
    file=model,
    data = model_data_list,
    chains = num_cores, 
    iter = 3000,
    warmup = 2600,
    refresh = 500  )


#  save the model fit
saveRDS(r_stan_fit, "/home/linarek/Documents/stan models/stan_fit.rds")

# print model results
summary(r_stan_fit, pars = c("alpha", "p"), probs = c(0.05, 0.95))$summary %>% as_tibble()

# plot trace plots
mcmc_trace(as.matrix(r_stan_fit,"alpha"))

# plot area plot
mcmc_areas(as.matrix(r_stan_fit,"alpha"))



# fix label switching

# Option 1

# variational bayes with cmdstanr
# fit_vb <- mod$variational(
#   data = model_data_list, 
#   iter = 25000,
#   elbo_samples = 1000,
#   algorithm = c("fullrank"),
#   output_samples = 10000,
#   tol_rel_obj = 0.00001)

# fit the model
 vb_fit <- vb(
     model,
     data = model_data_list,
     iter = 35000,
     elbo_samples = 1000,
     algorithm = c("fullrank"),
     output_samples = 10000,
     tol_rel_obj = 0.00001
   )
 
# vb estimates
print(vb_fit, c("alpha", "p"))



# option 2

# extract stan fit as the required format of the input
pars <- r_stan_fit %>% names %>% `[`(1:20) # important to include parameter estimates needed
r_stan_fit@model_pars

post_par <- rstan::extract(r_stan_fit,
                           c("alpha", "p", "pred_class", "pred_class_dis", "lp__"),
                           permuted = TRUE)

# simulated allocation vectors
post_class <- post_par$pred_class_dis
# classification probabilities
post_class_p <- post_par$pred_class


m = 4800 # of draws
K = 2 # of classes
J = 5 # of component-wise parameters

# initialize mcmc arrays
mcmc <- array(data = NA, dim = c(m = m, K = K, J = J))

# assign posterior draws to the array
mcmc[, , 1] <- post_par$alpha
for (i in 1:(J - 1)) {
  mcmc[, , i + 1] <- post_par$p[, , i]
}

# set of selected relabeling algorithm
set <-
  c("PRA",
    "ECR",
    "ECR-ITERATIVE-1",
    "AIC",
    "ECR-ITERATIVE-2",
    "STEPHENS",
    "DATA-BASED")

# find the MAP draw as a pivot
mapindex = which.max(post_par$lp__)

# switch labels
ls_lcm <- label.switching(
    method = set,
    zpivot = post_class[mapindex,],
    z = post_class,
    K = K,
    prapivot = mcmc[mapindex, ,],
    constraint = 1,
    mcmc = mcmc,
    p = post_class_p,
    data = model_data_list$y
  ) 

# permuted posterior based on ECR method
mcmc_permuted <- permute.mcmc(mcmc, ls_lcm$permutations$ECR)

# change dimension for each parameter defined as in the Stan code
mcmc_permuted <-
  array(
    data = mcmc_permuted$output,
    dim = c(4800, 1, 20),
    dimnames = list(NULL, NULL, pars)
  )


# reassess the model convergence after switch the labels
fit_permuted <-
  monitor(mcmc_permuted, warmup = 0,  digits_summary = 3)



sim_summary <- as.data.frame(fit_permuted)

estimated_values <- sim_summary[pars %>% sort(), c("mean", "2.5%", "97.5%")]


# area plot
mcmc_permuted %>% 
  as_tibble() %>% 
  select(starts_with('1.alpha')) %>% 
  mcmc_areas()



# Combine predictions with data set

# function for mode
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


# save out class predictions
class_prediction <- as.matrix(r_stan_fit)%>%
  as_tibble %>% 
  select(starts_with("pred_class_dis")) %>% 
  map_dfr(., function(x) {(calculate_mode(x))}
          ) %>% 
  data.table::transpose() %>% 
  rename(class_prediction=V1)


# merge class prediction to survey data
survey_data <- survey_data %>% add_column(class_prediction) 

# save combined data
write_rds(survey_data,"/home/linarek/mixture_models/model_results/survey_data_class_predictions.rds")

