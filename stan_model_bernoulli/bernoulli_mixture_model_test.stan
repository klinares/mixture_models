data {
  int<lower=1> I;                                     // number of items
  int<lower=1> J;                                     // number of respondents
  matrix[J, I] y;                                     // score for obs n
  int<lower=1> C;           // # of attribute profiles (latent classes) 
}
parameters {
  simplex[C] nu;
  // average intercept and main effect
  real mean_intercept;
  // item level deviations from average effects
  vector[I] dev_intercept;
  vector<lower=rep_array(-fabs(dev_intercept), C), upper=rep_array(fabs(dev_intercept), C)>[I] class_interaction[C];
  ordered[C] class_effect_non;
  
}

transformed parameters {
  matrix[I, C] log_prob;    // Probability of correct response for each class on each item
  // Probability of correct response for each class on each item
  for (c in 1:C) 
     log_prob[, c] = log_inv_logit(mean_intercept + dev_intercept + class_interaction[c] + class_effect_non[c]);
}
model{
  // Priors
  mean_intercept ~ normal(0, 5);
  
  for (c in 1:C)
    class_interaction[c] ~ std_normal();
  
  class_effect_non ~ normal(0, 3);
  dev_intercept ~ normal(0, 3);
 
  {
    for (j in 1:J) {
      row_vector[C] ps = log(nu') + y[j] * log_prob + (1 - y[j]) *  log1m_exp(log_prob);
      target += log_sum_exp(ps);
    }
  }
}


generated quantities {
    //int<lower = 1> pred_class_dis[J];     // posterior prediction for respondent j in latent class c
    simplex[C] pred_class[J];     // posterior probabilities of respondent j in latent class c
    matrix[I, C] probs;
    matrix[J, C] ps;

  
    for (c in 1:C){
      probs[, c] = inv_logit(mean_intercept + dev_intercept + class_interaction[c] + class_effect_non[c]);
       for (j in 1:J){
         pred_class[j] = exp((ps[, c])-log_sum_exp(ps));
       }
  }
}

