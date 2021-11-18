data {
  int<lower=1> I; // # of items
  int<lower=1> J; // # of respondents
  int<lower=1> C; // # of classes
  int y[J,I]; // response  matrix, [row, columns]
}

parameters {
  simplex[C] alpha; // vector of probabilities of being in one group
  
  real <lower = 0, upper = 1> p[C, I]; // [C, I] matrix of endoresment probabilities for I items & C classes
}

transformed parameters{
  vector[I] p_prod; // product of endorsing item i across the C classes
  
  for(i in 1:I){
    p_prod[i] = prod(p[, i]); 
  }
  
}


  // when no prior is included, a uniform prior distribution is the default
  // compute logs of the contributions to the marginal probabilities from each class, same in lmix

model {
  real lmix[C];
  
  for (j in 1:J){
    for (c in 1: C){
      lmix[c] = log(alpha[c]) + bernoulli_lpmf(y[j, ] | p[c,]);
    }
    
    // log_sum_exp computes the logarithm of the sum of the exponentiated elements of lmix
    target += log_sum_exp(lmix); 
  }
}

  // prediction of latent class membership
  // approximation of the expectation of our posterior prediction of class membership
generated quantities {
  int<lower = 1> pred_class_dis[J];     // posterior prediction for respondent j in latent class c
  simplex[C] pred_class[J];     // posterior probabilities of respondent j in latent class c
  real lmix[C];
  
  
  for (j in 1:J){
    for (c in 1: C){
      //  The log Bernoulli probability mass of y given chance of success theta
      lmix[c] = log(alpha[c]) + bernoulli_lpmf(y[j, ] | p[c,]);
    }               
    for (c in 1: C){
      pred_class[j][c] = exp((lmix[c])-log_sum_exp(lmix));
    }
    // categorical_rng = Generate a categorical variate with N-simplex distribution parameter theta
    pred_class_dis[j] = categorical_rng(pred_class[j]); 
  }
}

