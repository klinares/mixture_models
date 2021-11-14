### Model Results:
- Results can be found on [rpubs](https://rpubs.com/klinares/834991). The idea of the LCA model was to distinguish "classes" of individuals that can explain the missingness in the data. This type of model is better suited for when modeling a particular response (i.e., "Don't know", "Approval", etc.) rather than item nonresponse because other factors contribute to missingness in the data.
- 9 abortion attitudes were used as outcome variables in this bayes LCA model.
  - Each outcome variable was coded as NA = 1, other values  = 0.
  - A bernoulli distribution was applied to each outcome variable within [stan](https://github.com/klinares/mixture_models/blob/main/stan_model_bernoulli/bernoulli_mixture_model.stan).
  - bayesLCA an R package was used to explore how many classes to model. cmdstanr was used in the final model, script found [here](https://github.com/klinares/mixture_models/blob/main/stan_model_bernoulli/bernoulli_mixture_modeling.R). 
- The modeling idea was meant to identify patterns of missing data in the abortion attitudes data. Although a two class solution was selected, rhats did not show evidence of congergence. Additionally, traceplots did not show evidence of class switching, which is a common issue in LCA models. Item probabilities were small suggesting that perhaps there is not pattern of missingness in the data that can be attributed to a general "class" of respondents.
- From the posterior distribution of each draw, a categorical class prediction was computed at each draw between {1, 2} for a 2-class solution.
- The mode was taken for each set of draws and a class membership prediction was created for each survey respondent.
- Class membership predictions were added to the origin survey data to use within a group_by function and compute proportions of demographics to examine whether there were differences in the demographic composition of those in classes that were likely to have more missing data.
